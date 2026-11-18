"""YaMAS output as an upload: turn the files YaMAS writes into our shared counts format.

YaMAS (the lab's pipeline from raw reads to tables) writes, for 16S data, two files:

- ``otu.csv``: one row per OTU/ASV (its ID in the first column) and one column per sample,
  holding read counts;
- ``taxonomy.csv``: one row per OTU ID, with the lineage in a ``Taxon`` column (QIIME 2 style;
  other columns, e.g. ``Confidence``, are ignored).

``counts_from_16s`` joins them into the shared format (sample IDs in the first column, one column
per OTU, and a last ``taxonomy`` row with each OTU's lineage), standardising every lineage with the
lab's ``normalize_taxonomy_name``. The **read counts are kept as they are** (whole numbers): no
division by the sample total and no percentages, so the data check, the quality figures and MIPMLP
all see raw counts. OTUs that end up with the same lineage stay separate columns: MIPMLP groups them
itself (its sub-PCA uses the variation among them), and richness counts OTUs.

The pair is recognised by its shape (``find``), before the uploads are sorted into counts, labels and
metabolites. Each YaMAS format is one entry of ``FORMATS``; shotgun output will be added as another
entry once it is built against a real YaMAS shotgun file.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable

import pandas as pd

from microplatform.core.contract import ID_COLUMN, TAXONOMY_ROW
from microplatform.core.errors import UserFacingError

TAXON_COLUMN = "Taxon"  # the lineage column of taxonomy.csv (QIIME 2's name), matched ignoring case
UNASSIGNED = "Unassigned"


class YamasError(UserFacingError):
    """The YaMAS files cannot be combined as they are."""


# --- the lab's taxonomy standardisation ----------------------------------------------------------
# Copied unchanged from the lab's bridge script, docs/preprocess_16S.py (normalize_taxonomy_name and the
# three constants it uses), so names are standardised exactly as the lab does it.

RANKS = ("k", "p", "c", "o", "f", "g", "s")
D_RANKS = {f"D_{i}": rank for i, rank in enumerate(RANKS)}
MISSING_TAXON_VALUES = {"", "nan", "none", "unassigned", "unknown"}


def normalize_taxonomy_name(taxon):
    """Return one canonical taxonomy string for the formats used by YAMAS.

    Rank prefixes are converted to ``k__/p__/...`` and numeric confidence
    annotations such as ``Bacteria(0.71)`` are removed. Empty ranked levels
    are retained so that a missing level cannot shift all subsequent ranks.
    """
    value = str(taxon).strip()
    if value.lower() in MISSING_TAXON_VALUES:
        return "Unassigned"

    parts = re.split(r"\s*[;,|]\s*", value)
    ranked = {}
    plain = []

    for part in parts:
        part = part.strip()
        if not part or part.lower() in MISSING_TAXON_VALUES:
            continue

        match = re.match(
            r"^(?:(D_[0-6])__|([dkpcofgs])(?:__|:))\s*(.*)$",
            part,
            flags=re.IGNORECASE,
        )

        if match:
            rank = D_RANKS.get(
                (match.group(1) or "").upper(),
                (match.group(2) or "").lower(),
            )
            name = match.group(3).strip()
        else:
            rank = None
            name = part

        # Some classifiers append a per-rank confidence score to the name.
        # Restrict removal to a trailing number in [0, 1], so meaningful
        # parenthesized text in a taxon name is preserved.
        name = re.sub(
            r"\s*\((?:0(?:\.\d+)?|1(?:\.0+)?)\)\s*$",
            "",
            name,
        ).strip()

        if rank == "d":
            rank = "k"

        if rank:
            ranked[rank] = name
        elif name:
            plain.append(name)

    if ranked:
        deepest_rank = max(RANKS.index(rank) for rank in ranked)
        return ";".join(
            f"{rank}__{ranked.get(rank, '')}"
            for rank in RANKS[:deepest_rank + 1]
        )

    return ";".join(plain) if plain else "Unassigned"


# --- 16S: otu.csv + taxonomy.csv -----------------------------------------------------------------

@dataclass(frozen=True)
class Combined:
    """YaMAS files turned into one counts table in the shared format."""

    counts: pd.DataFrame  # ID column, one column per OTU, a last 'taxonomy' row
    layout: str  # e.g. "yamas-16s" (the run records it; the data check says how the files were read)
    files: tuple[str, ...]  # the uploaded file names, in the order the format names them
    details: dict = field(default_factory=dict)  # numbers for the data check: OTUs, samples, unassigned


def taxon_column(table: pd.DataFrame) -> str | None:
    """The ``Taxon`` column of a taxonomy table (any case), or None."""
    return next((c for c in table.columns[1:] if str(c).strip().lower() == TAXON_COLUMN.lower()), None)


def _find_16s(tables: list[tuple[str, pd.DataFrame]]) -> tuple[int, int] | None:
    """Positions of (otu table, taxonomy table) among the uploads, or None if there is no taxonomy table.

    The taxonomy table has a ``Taxon`` column; the OTU table is the other upload whose first column holds
    the taxonomy's OTU IDs (most of its rows; ``counts_from_16s`` then says which ones are missing).
    """
    taxonomies = [i for i, (_, t) in enumerate(tables) if taxon_column(t) is not None]
    if not taxonomies:
        return None
    if len(taxonomies) > 1:
        names = " and ".join(f"'{tables[i][0]}'" for i in taxonomies)
        raise YamasError(f"{names} both look like YaMAS taxonomy files (a '{TAXON_COLUMN}' column).",
                         "Upload one otu.csv and its taxonomy.csv at a time.")
    tax_at = taxonomies[0]
    otu_ids = set(_ids(tables[tax_at][1]))
    overlap = {i: len(set(_ids(t)) & otu_ids) / max(len(t), 1)
               for i, (_, t) in enumerate(tables) if i != tax_at}
    best = max(overlap, key=overlap.get, default=None)
    if best is None or overlap[best] < 0.5:
        raise YamasError(
            f"'{tables[tax_at][0]}' looks like a YaMAS taxonomy file (a '{TAXON_COLUMN}' column), but no "
            "uploaded table lists its OTUs.",
            "Upload it together with YaMAS's otu.csv: one row per OTU (the same IDs in the first column) and "
            "one column per sample.")
    return best, tax_at


def counts_from_16s(otu: pd.DataFrame, taxonomy: pd.DataFrame, names: tuple[str, str] = ("otu.csv",
                    "taxonomy.csv")) -> tuple[pd.DataFrame, dict]:
    """The shared-format counts table from YaMAS's 16S files, and numbers for the data check.

    ``otu``: OTU IDs in the first column, then one column per sample (read counts). ``taxonomy``:
    OTU IDs in the first column and a ``Taxon`` column. The counts are kept as read counts.
    """
    otu_name, tax_name = names
    column = taxon_column(taxonomy)
    if column is None:
        raise YamasError(f"'{tax_name}' has no '{TAXON_COLUMN}' column.",
                         f"Upload YaMAS's taxonomy.csv: OTU IDs first, then a '{TAXON_COLUMN}' column.")
    if otu.shape[1] < 2:
        raise YamasError(f"'{otu_name}' has no sample columns.",
                         "Upload YaMAS's otu.csv: OTU IDs in the first column, then one column per sample.")
    ids = _ids(otu)
    duplicated = sorted({i for i in ids if ids.count(i) > 1})
    if duplicated:
        raise YamasError(f"Some OTU IDs appear more than once in '{otu_name}': {', '.join(duplicated[:3])}.",
                         "Each OTU should have one row, as YaMAS writes it.")
    lineages = pd.Series(taxonomy[column].to_numpy(), index=_ids(taxonomy))
    lineages = lineages[~lineages.index.duplicated()]
    missing = [i for i in ids if i not in lineages.index]
    if missing:
        raise YamasError(
            f"{len(missing):,} OTU{'s' if len(missing) != 1 else ''} in '{otu_name}' "
            f"{'are' if len(missing) != 1 else 'is'} not in '{tax_name}' (for example {', '.join(missing[:3])}).",
            "Upload the otu.csv and taxonomy.csv from the same YaMAS run.")

    raw = otu.iloc[:, 1:]
    values = raw.apply(pd.to_numeric, errors="coerce")
    bad = values.isna() & raw.notna()
    if bad.to_numpy().any():
        row, col = next(zip(*bad.to_numpy().nonzero()))
        raise YamasError(
            f"'{otu_name}' has entries that are not numbers (for example, OTU '{ids[row]}' in sample "
            f"'{raw.columns[col]}' contains '{raw.iat[row, col]}').",
            "Every cell except the OTU IDs must be a read count.")
    values = values.fillna(0)
    array = values.to_numpy(dtype=float)
    if (array < 0).any() or (array % 1 != 0).any():
        raise YamasError(
            f"'{otu_name}' does not hold read counts (some values are negative or not whole numbers).",
            "Upload YaMAS's otu.csv as it was written (whole-number read counts), not a normalised table.")

    samples = [str(s) for s in values.columns]
    duplicated = sorted({s for s in samples if samples.count(s) > 1})
    if duplicated:
        raise YamasError(f"Some sample IDs appear more than once in '{otu_name}': {', '.join(duplicated[:3])}.",
                         "Give every sample column its own ID.")
    standard = [normalize_taxonomy_name(lineages[i]) for i in ids]
    counts = pd.DataFrame(array.astype("int64").T, columns=ids)
    table = pd.concat([pd.DataFrame({ID_COLUMN: samples}), counts], axis=1)
    table = pd.concat([table, pd.DataFrame([[TAXONOMY_ROW, *standard]], columns=table.columns)],
                      ignore_index=True)
    details = {"otus": len(ids), "samples": len(samples),
               "unassigned": sum(s == UNASSIGNED for s in standard),
               "lineages": len(set(standard)), **replicate_suffixes(samples)}
    return table, details


def replicate_suffixes(samples: list[str]) -> dict:
    """How many sample IDs end in ``_1`` or ``_2``, and how many names appear with both (only reported, for
    the data check: the lab's bridge script drops ``_2`` and strips ``_1``, which is not done here until the
    meaning of the suffixes in YaMAS output is settled). Empty if there are none."""
    ones = {s[:-2] for s in samples if s.endswith("_1")}
    twos = {s[:-2] for s in samples if s.endswith("_2")}
    if not (ones or twos):
        return {}
    return {"suffix_1": len(ones), "suffix_2": len(twos), "suffix_both": len(ones & twos)}


@dataclass(frozen=True)
class YamasFormat:
    """One kind of YaMAS output: how to find its files among the uploads, and how to combine them."""

    layout: str
    description: str  # for the data check, e.g. "YaMAS 16S output (otu.csv + taxonomy.csv)"
    find: Callable[[list[tuple[str, pd.DataFrame]]], tuple[int, ...] | None]
    combine: Callable[..., tuple[pd.DataFrame, dict]]


FORMATS: tuple[YamasFormat, ...] = (
    YamasFormat("yamas-16s", "YaMAS 16S output (otu.csv + taxonomy.csv)", _find_16s, counts_from_16s),
    # Shotgun: to be added here, built against a real YaMAS shotgun file.
)
LAYOUTS = {f.layout: f.description for f in FORMATS}


def combine(tables: list[tuple[str, pd.DataFrame]]) -> tuple[Combined, list[int]] | None:
    """YaMAS files among the uploads ``(name, table)`` combined into one counts table, and the positions
    of the uploads used; None if there are none."""
    for fmt in FORMATS:
        found = fmt.find(tables)
        if found is None:
            continue
        names = tuple(tables[i][0] for i in found)
        counts, details = fmt.combine(*(tables[i][1] for i in found), names=names)
        return Combined(counts, fmt.layout, names, details), list(found)
    return None


def _ids(table: pd.DataFrame) -> list[str]:
    return [str(v).strip() for v in table.iloc[:, 0]]
