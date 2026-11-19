"""MetaPhlAn 4 output (shotgun) as an upload: keep the species, in the shared counts format.

MetaPhlAn lists every clade at every rank: ``k__Bacteria``, ``k__Bacteria|p__Firmicutes``, … down to
the species (``…|s__Ruminococcus_torques``) and its SGB (``…|t__SGB4608``). Each rank adds up to the
same total, so keeping all of them would count the same reads up to eight times. Only the **species**
are kept (a lineage whose last rank is ``s__``); the SGB rows split a species further and the higher
ranks are sums of the species. Lineages are standardised with the lab's ``normalize_taxonomy_name``
(``k__…;p__…;…;s__…``), and an ``UNCLASSIFIED`` row, if any, is left out (its share is reported).

The values are recognised as one of two forms:

- **read counts**: MetaPhlAn's estimated reads per clade (``-t rel_ab_w_read_stats``, the column
  ``estimated_number_of_reads_from_the_clade``, or a merged table of them): whole numbers. They are kept,
  so the data check sees the sequencing depth and MIPMLP normalises them as it does 16S counts;
- **relative abundance** (MetaPhlAn's default, %): passed on as it is. The data check then says plainly
  what that limits: no sequencing depth, and MIPMLP's log (log10 of value + 0.1) brings species below
  about 0.1% close to the value for "absent".

Layouts accepted (all recognised by MetaPhlAn's ``#mpa_…`` line or a ``clade_name`` header):

- clades as rows, one column per sample (``merge_metaphlan_tables.py``), optionally with
  ``clade_taxid``/``NCBI_tax_id`` columns (ignored);
- samples as rows, one column per clade (a transposed merged table: ``clade_name`` heads the sample IDs);
- one sample's profile (``#clade_name``, ``clade_taxid``, ``relative_abundance``, … and, with read stats,
  ``estimated_number_of_reads_from_the_clade``, which is then used). The sample is named by MetaPhlAn's
  ``#SampleID`` line, else by the file name.
"""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

from microplatform.core.contract import ID_COLUMN, RANK_SEPARATOR, TAXONOMY_ROW
from microplatform.core.errors import UserFacingError
from microplatform.core.yamas import normalize_taxonomy_name, replicate_suffixes

LAYOUT = "metaphlan"
DESCRIPTION = "MetaPhlAn 4 output (species level)"
CLADE = "clade_name"
READS = "estimated_number_of_reads_from_the_clade"
RELATIVE = "relative_abundance"
IGNORED = {"clade_taxid", "ncbi_tax_id", "coverage", "additional_species"}  # per-clade columns, not samples
UNCLASSIFIED = "unclassified"
LOW_PERCENT = 0.1  # MIPMLP adds 0.1 before its log: % values below this end up close to "absent"
_RANKS = re.compile(r"\s*[|,;]\s*")


class MetaphlanError(UserFacingError):
    """A MetaPhlAn file that cannot be used as it is."""


def is_metaphlan(content: bytes) -> bool:
    """MetaPhlAn output: an ``#mpa_…`` line, or a first column headed ``clade_name`` (``#clade_name``)."""
    for line in _lines(content)[:10]:
        if line.startswith("#mpa_"):
            return True
        first = _split(line)[0].strip().lstrip("#").strip().lower() if line.strip() else ""
        if first == CLADE:
            return True
    return False


def to_shared_format(name: str, content: bytes) -> tuple[pd.DataFrame, dict] | None:
    """The species as a shared-format counts table (ID column, one column per species, a last 'taxonomy'
    row), and numbers for the data check; None if the file is not MetaPhlAn output."""
    if not is_metaphlan(content):
        return None
    lines = _lines(content)
    sample_line = next((ln for ln in lines if ln.lower().startswith("#sampleid")), None)
    header_at = next((i for i, ln in enumerate(lines) if _split(ln)[0].strip().lstrip("#").strip().lower() == CLADE),
                     None)
    if header_at is None:
        raise MetaphlanError(f"'{name}' looks like MetaPhlAn output, but has no '{CLADE}' header line.",
                             "Upload the table as MetaPhlAn (or merge_metaphlan_tables.py) wrote it.")
    body = [lines[header_at].lstrip("#")] + [ln for ln in lines[header_at + 1:] if ln.strip() and not ln.startswith("#")]
    rows = list(csv.reader(body, delimiter=_delimiter(body[0])))
    header = [h.strip() for h in rows[0]]
    if len(rows) < 2 or len(header) < 2:
        raise MetaphlanError(f"'{name}' has no clades.", "Upload the MetaPhlAn table with its rows of clades.")

    if sum(_is_lineage(h) for h in header[1:]) >= 0.9 * (len(header) - 1):
        # samples as rows, clades as columns
        samples = [r[0].strip() for r in rows[1:]]
        values = pd.DataFrame([r[1:] for r in rows[1:]], columns=header[1:])
        clades = header[1:]
        orientation = "samples as rows"
    else:
        # clades as rows: the sample columns, or one profile's read (or relative abundance) column
        clades = [r[0].strip() for r in rows[1:]]
        table = pd.DataFrame([r[1:] for r in rows[1:]], columns=header[1:])
        lower = {c.lower(): c for c in table.columns}
        if READS in lower or RELATIVE in lower:
            column = lower.get(READS) or lower[RELATIVE]
            sample = _sample_name(sample_line) or Path(name).stem
            table = table[[column]].set_axis([sample], axis=1)
            orientation = "one sample's profile"
        else:
            table = table[[c for c in table.columns if c.lower() not in IGNORED]]
            orientation = "clades as rows"
        samples = [str(c) for c in table.columns]
        values = table.T
        values.columns = clades
    numbers = values.apply(lambda c: pd.to_numeric(c, errors="coerce"))
    bad = numbers.isna() & values.notna() & (values.astype(str).apply(lambda c: c.str.strip()) != "")
    if bad.to_numpy().any():
        row, col = next(zip(*bad.to_numpy().nonzero()))
        raise MetaphlanError(
            f"'{name}' has entries that are not numbers (for example '{values.iat[row, col]}' for sample "
            f"'{samples[row]}').", "Every value must be a relative abundance or an estimated read count.")
    numbers = numbers.fillna(0)
    numbers.index = samples
    if len(set(samples)) != len(samples):
        twice = sorted({s for s in samples if samples.count(s) > 1})
        raise MetaphlanError(f"Some sample IDs appear more than once in '{name}': {', '.join(twice[:3])}.",
                             "Give every sample its own ID.")

    ranks = [_ranks(c) for c in clades]
    species = [i for i, r in enumerate(ranks) if r and r[-1].lower().startswith("s__")]
    unclassified = [i for i, c in enumerate(clades) if c.strip().lower() == UNCLASSIFIED]
    if not species:
        raise MetaphlanError(
            f"'{name}' has no species-level clades (lineages ending in s__…).",
            "Upload MetaPhlAn's full output (all ranks, as it writes it), not a table cut to a higher rank.")
    if (numbers.to_numpy() < 0).any():
        raise MetaphlanError(f"'{name}' has negative values.",
                             "Upload MetaPhlAn's output as it was written, not a transformed table.")

    kept = numbers.iloc[:, species]
    array = kept.to_numpy(dtype=float)
    totals = array.sum(axis=1)
    whole = bool((array % 1 == 0).all())
    near_100 = bool(len(totals) and (np.abs(totals[totals > 0] - 100) < 1).all())
    values_kind = "counts" if whole and not near_100 else "relative"
    if values_kind == "relative" and (array > 100 + 1e-6).any():
        raise MetaphlanError(
            f"'{name}' holds values above 100 that are not whole numbers, so they are neither percentages "
            "nor read counts.", "Upload MetaPhlAn's relative abundances (%) or its estimated read counts.")
    if whole:
        array = array.astype("int64")

    lineages = [normalize_taxonomy_name(RANK_SEPARATOR.join(ranks[i])) for i in species]
    shared = pd.DataFrame(array, columns=lineages)
    shared.insert(0, ID_COLUMN, samples)
    shared = pd.concat([shared, pd.DataFrame([[TAXONOMY_ROW, *lineages]], columns=shared.columns)],
                       ignore_index=True)
    all_values = numbers.to_numpy(dtype=float)
    unclassified_share = 0.0
    if unclassified:
        # the unclassified row and the species are shares of the same total (MetaPhlAn's --unclassified_estimation)
        u = all_values[:, unclassified].sum(axis=1)
        whole_total = u + totals
        unclassified_share = float(np.mean(np.divide(u, whole_total, out=np.zeros_like(u), where=whole_total > 0)) * 100)
    kingdoms = pd.Series([lin.split(RANK_SEPARATOR)[0].removeprefix("k__") for lin in lineages]).value_counts()
    details = {
        "samples": len(samples), "clades": len(clades), "species": len(species),
        "values": values_kind,
        "orientation": orientation,
        "kingdoms": {str(k): int(v) for k, v in kingdoms.items()},
        "unclassified_percent": round(unclassified_share, 2),
        "database": next((ln.lstrip("#").strip() for ln in lines[:5] if ln.startswith("#mpa_")), ""),
        **replicate_suffixes(samples),
    }
    if values_kind == "relative":
        present = array[array > 0]
        details["below_low_percent"] = round(float((present < LOW_PERCENT).mean()) * 100, 1) if present.size else 0.0
    return shared, details


def _lines(content: bytes) -> list[str]:
    return content.decode("utf-8-sig", errors="replace").splitlines()


def _delimiter(line: str) -> str:
    return "\t" if "\t" in line else ","


def _split(line: str) -> list[str]:
    return next(csv.reader([line], delimiter=_delimiter(line)), [""]) or [""]


def _ranks(clade: str) -> list[str]:
    return [p for p in _RANKS.split(str(clade).strip()) if p]


def _is_lineage(name: str) -> bool:
    ranks = _ranks(name)
    return bool(ranks) and all(re.match(r"^[a-z]__", r) for r in ranks)


def _sample_name(line: str | None) -> str | None:
    if not line:
        return None
    parts = [p.strip() for p in re.split(r"\t|,", line.lstrip("#")) if p.strip()]
    return parts[1] if len(parts) > 1 and parts[1].lower() != "metaphlan_analysis" else None
