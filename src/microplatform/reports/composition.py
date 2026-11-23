"""Phylum composition, in the lab's figure style (the fig_1B scripts in references/composition/).

From the user's own abundance table (counts or already relative, never the log-transformed
table the tools use):

1. **Phylum table** (``phylum_percent``): each column's phylum is read from its lineage
   (``k__Bacteria;p__Firmicutes;…``; ``;``, ``|`` and ``,`` all separate ranks). A table that
   lists every rank (MetaPhlAn: ``k__…``, ``k__…|p__…``, ``k__…|p__…|c__…``) uses only its
   phylum rows, so nothing is counted twice. Phyla are summed, in % of each sample's whole total;
   columns with no phylum (``Bacteria``, ``Archaea``, ``Unassigned``, ``Unclassified``, and the names
   ``NO_NAME`` lists) are the segment "Not named at phylum", as at the finer ranks, so every bar adds
   up to 100%. Samples with a total of zero are left out. Values must be finite
   and not negative. For display, one phylum gets one bar (unless ``merge_gtdb_variants=False``;
   ``MERGE_GTDB_VARIANTS``): GTDB variants (``Bacillota_A_368345``, ``Bacillota_C``) are merged into their
   base phylum (``Bacillota``), and an old and a new name of the same phylum (``KNOWN_PHYLA``: Firmicutes
   and Bacillota, as MetaPhlAn's database mixes them) into the one the table uses more. The analysis
   itself (MIPMLP and the tools) keeps the real names. A phylum-rank ``…_unclassified`` (MetaPhlAn's
   ``p__Bacteria_unclassified``) is not a phylum: it is counted in "Other", so every bar still adds up to 100%.
2. **Top phyla** (``top_phyla``): the ``TOP_PHYLA`` phyla with the highest mean over all
   samples, most abundant first; the rest are summed as "Other", drawn last in grey.
3. **Colours** (``phylum_colors``): one colour per phylum, the same in every panel, figure and
   dataset. Common gut phyla have a fixed slot (``KNOWN_PHYLA``; old and new names share it),
   so Firmicutes looks the same in every run; other phyla take the next free colour.
   ``palette="lab"`` (the default) is the lab's look, tab20 + Set3 + Pastel1;
   ``palette="colorblind"`` is a colour-blind-safer set of 12 hues.
4. **Figures**: ``composition_figure`` (left: one stacked bar per sample, sorted by Firmicutes
   if present, else by the top phylum; right: the mean of all samples) and, only when asked
   and the samples form two groups, ``split_violin_figure`` (per phylum: group A on the left
   half, group B on the right, medians marked).
"""

from __future__ import annotations

import re

import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import to_hex

from microplatform.core.errors import UserFacingError
from microplatform.reports.figures import OTHER, Figure

TOP_PHYLA = 12
OTHER_GREY = "#BDBDBD"
MEAN = "Mean"  # the right panel's one bar
DROPPED = ("Bacteria", "Archaea", "Unassigned", "Unclassified")  # no phylum: left out (case-insensitive)
SORT_BY = ("Firmicutes", "Bacillota")  # the same phylum, old and new name
# GTDB splits some phyla into variants named by a capital-letter suffix and/or a numeric ID
# (Bacillota_A_368345, Bacillota_C, Desulfobacterota_G_459543, Methanobacteriota_A_1229). For display they
# are merged into the base phylum, so a phylum has one bar and its usual colour. Set to False to keep them.
MERGE_GTDB_VARIANTS = True
_GTDB_VARIANT = re.compile(r"^(?P<base>.+?)(?:_[A-Z])?(?:_\d+)?$")

def _not_grey(color: str) -> bool:
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return max(r, g, b) - min(r, g, b) > 16


PALETTES = {
    # The lab's look: matplotlib's tab20, then Set3, then Pastel1 (as in the fig_1B scripts), without
    # their greys, which would read as "Other".
    "lab": tuple(c for name in ("tab20", "Set3", "Pastel1") for c in map(to_hex, colormaps[name].colors)
                 if _not_grey(c)),
    # Colour-blind safer, 12 hues: Okabe-Ito (without black) + five of Paul Tol's muted/dark hues.
    # No 12-colour set is safe for everyone; neighbours in a stack differ in lightness too.
    "colorblind": ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442",
                   "#332288", "#882255", "#999933", "#117733", "#661100"),
}
PALETTE_DEFAULT = "lab"
# Fixed colour slots, in the usual order of abundance in gut samples, so that a phylum keeps its
# colour across datasets. Each entry lists one phylum's names (older name first, then current ones).
KNOWN_PHYLA = (
    ("Firmicutes", "Bacillota"), ("Bacteroidetes", "Bacteroidota"), ("Proteobacteria", "Pseudomonadota"),
    ("Actinobacteria", "Actinobacteriota", "Actinomycetota"), ("Verrucomicrobia", "Verrucomicrobiota"),
    ("Euryarchaeota", "Methanobacteriota"), ("Tenericutes", "Mycoplasmatota"), ("Fusobacteria", "Fusobacteriota"),
    ("Cyanobacteria", "Cyanobacteriota"), ("Spirochaetes", "Spirochaetota"), ("Synergistetes", "Synergistota"),
    ("Lentisphaerae", "Lentisphaerota"), ("Desulfobacterota",), ("Campylobacterota",),
    ("Patescibacteria", "Saccharibacteria", "TM7"), ("Chloroflexi", "Chloroflexota"),
    ("Elusimicrobia", "Elusimicrobiota"), ("Deferribacteres", "Deferribacterota"),
    ("Planctomycetes", "Planctomycetota"), ("Acidobacteria", "Acidobacteriota"),
)
_SLOT = {name.lower(): i for i, names in enumerate(KNOWN_PHYLA) for name in names}
_RANKS = re.compile(r"\s*[;|,]\s*")
_PREFIX = re.compile(r"^[a-z]__")

__all__ = ["phylum_percent", "gtdb_base", "merged_note", "top_phyla", "phylum_colors", "composition_figure",
           "split_violin_figure", "rank_percent", "top_taxa", "fine_colors", "fine_composition_figure",
           "taxon_choices", "taxon_abundance"]


class CompositionError(UserFacingError):
    """The table cannot be shown as relative abundance (e.g. negative, log-transformed values)."""


def gtdb_base(phylum: str) -> str:
    """'Bacillota_A_368345' -> 'Bacillota' (GTDB variant suffixes removed); other names unchanged."""
    match = _GTDB_VARIANT.match(str(phylum))
    return match["base"] if match else str(phylum)


def phylum_percent(table: pd.DataFrame,
                   merge_gtdb_variants: bool = MERGE_GTDB_VARIANTS) -> tuple[pd.DataFrame, dict]:
    """Samples × phyla in % of each sample's whole total, as at the finer ranks: abundance with no phylum is
    the column "Not named at phylum" (only if some sample has any), so every row adds up to 100.

    ``table`` is samples × taxa, columns named by lineage, values counts or relative abundances.
    The second value says how many samples had a total of zero (left out), the mean share of a sample with
    no phylum (``unassigned_percent``, = ``not_named_percent``), and which GTDB variants and old/new names
    were merged (``merged``; empty if none).
    """
    phylum, totals, merged = _phylum_sums(table, merge_gtdb_variants)
    valid = totals > 0
    if not (phylum.loc[valid] > 0).to_numpy().any():
        raise CompositionError("No sample has any abundance at phylum level.",
                               "Check that the column names hold the taxonomy (k__…;p__…).")
    percent = phylum.loc[valid].div(totals[valid], axis=0) * 100
    gap = (100 - percent.sum(axis=1)).clip(lower=0)
    gap = gap.where(gap > NOISE, 0.0)  # rounding left over is not abundance
    if (gap > 0).any():
        percent[not_named("phylum")] = gap
    share = float(gap.mean())
    return percent, {"samples_dropped": int((~valid).sum()), "unassigned_percent": share,
                     "not_named_percent": share, "merged": merged}


def _phylum_sums(table: pd.DataFrame, merge_gtdb_variants: bool) -> tuple[pd.DataFrame, pd.Series, dict]:
    """Samples × phyla in the table's own values (a phylum-rank ``…_unclassified`` summed as "Other"), each
    sample's whole total, and the merges (``phylum_percent``). The phylum figures and the one-taxon view
    both use it, so a phylum has the same name and the same share in both."""
    values = table.apply(pd.to_numeric, errors="coerce")
    array = values.to_numpy(dtype=float)
    if not (pd.notna(array).all() and abs(array).max(initial=0) < float("inf")) or (array < 0).any():
        raise CompositionError(
            "Composition needs finite, non-negative abundances (counts or relative abundances).",
            "Upload the raw counts or relative abundances, not a log-transformed table.",
        )
    ranks = [_ranks(c) for c in values.columns]
    if _all_ranks(ranks):  # MetaPhlAn style: keep the phylum rows only
        keep = [i for i, r in enumerate(ranks) if r and r[-1].startswith("p__")]
    else:
        keep = list(range(len(ranks)))
    names = [_phylum(ranks[i]) for i in keep]
    variants: dict[str, str] = {}
    if merge_gtdb_variants:
        variants = {n: gtdb_base(n) for n in dict.fromkeys(names) if gtdb_base(n) != n}
        names = [variants.get(n, n) for n in names]
    phylum = values.iloc[:, keep].T.groupby(names, sort=False).sum().T
    synonyms: dict[str, str] = {}
    if merge_gtdb_variants:
        by_slot: dict[int, list[str]] = {}
        for name in phylum.columns:
            if (slot := _SLOT.get(str(name).lower())) is not None:
                by_slot.setdefault(slot, []).append(name)
        for same in by_slot.values():
            if len(same) > 1:
                main = max(same, key=lambda n: phylum[n].sum())
                synonyms |= {n: main for n in same if n != main}
        if synonyms:
            phylum = phylum.T.groupby([synonyms.get(c, c) for c in phylum.columns], sort=False).sum().T
    totals = values.sum(axis=1) if not _all_ranks(ranks) else phylum.sum(axis=1)
    no_phylum = {d.lower() for d in DROPPED} | NO_NAME  # "Not named at phylum", as at the finer ranks
    phylum = phylum.drop(columns=[c for c in phylum.columns if str(c).strip().lower() in no_phylum])
    unclassified = [c for c in phylum.columns if str(c).lower().endswith("_unclassified")]
    if unclassified:  # abundance known to be bacterial but with no phylum: shown in "Other"
        phylum = phylum.T.groupby([OTHER if c in unclassified else c for c in phylum.columns], sort=False).sum().T
    return phylum, totals, {"variants": variants, "synonyms": synonyms}


def merged_note(merged: dict[str, dict[str, str]]) -> str:
    """One short sentence for the captions, with one example of each kind of merge ('' if none)."""
    kinds = []
    for key, what in (("variants", "GTDB variants"), ("synonyms", "old and new names")):
        pairs = merged.get(key) or {}
        if pairs:
            variant, base = next(iter(pairs.items()))
            kinds.append(f"{what} (e.g. {variant} → {base})")
    if not kinds:
        return ""
    return f" Merged for display: {' and '.join(kinds)}; the analysis keeps the original names."


def top_phyla(percent: pd.DataFrame, top_n: int = TOP_PHYLA) -> pd.DataFrame:
    """The ``top_n`` phyla by mean over all samples (most abundant first), then "Other" (the rest), then
    "Not named at phylum" (if the table has any)."""
    gap = not_named("phylum")
    order = list(percent.drop(columns=[gap], errors="ignore").mean()
                 .sort_values(ascending=False, kind="stable").index)
    named = [p for p in order if p != OTHER][:top_n]
    out = percent[named].copy()
    rest = [p for p in percent.columns if p not in named and p != gap]
    if rest:
        out[OTHER] = percent[rest].sum(axis=1)
    if gap in percent and (percent[gap] > 0).any():
        out[gap] = percent[gap]
    return out


def phylum_colors(phyla, palette: str = PALETTE_DEFAULT) -> dict[str, str]:
    """One colour per phylum (in the given order); "Other" is always grey.

    Known phyla keep their fixed slot; other phyla take the next free colour (after the fixed
    slots in the lab palette, so they never borrow a known phylum's colour).
    """
    if palette not in PALETTES:
        raise ValueError(f"palette must be one of {sorted(PALETTES)}")
    colors = PALETTES[palette]
    named = [p for p in dict.fromkeys(phyla) if p != OTHER and not _is_not_named(p)]
    unnamed = {p: NOT_NAMED_GREY for p in dict.fromkeys(phyla) if _is_not_named(p)}
    slots: dict[str, int] = {}
    for p in named:
        slot = _SLOT.get(str(p).lower())
        if slot is not None and slot < len(colors) and slot not in slots.values():
            slots[p] = slot
    start = len(KNOWN_PHYLA) if len(colors) > len(KNOWN_PHYLA) + len(named) else 0
    free = [i for i in range(start, len(colors)) if i not in slots.values()] + \
           [i for i in range(start) if i not in slots.values()]
    for p in named:
        if p not in slots:
            slots[p] = free.pop(0) if free else len(slots) % len(colors)
    return {p: colors[slots[p]] for p in named} | {OTHER: OTHER_GREY} | unnamed


def sort_samples(top: pd.DataFrame) -> pd.DataFrame:
    """Samples by Firmicutes (most first) if present, else by the most abundant phylum."""
    lead = next((p for p in top.columns if p in SORT_BY), top.columns[0])
    return top.sort_values(lead, ascending=False, kind="stable")


def composition_figure(top: pd.DataFrame, colors: dict[str, str], dataset: str | None = None,
                       left_out: str = "", descriptive: str = "") -> Figure:
    """Left: one stacked bar per sample; right: the mean of all samples. ``top`` from ``top_phyla``."""
    ordered = sort_samples(top)
    lead = next((p for p in ordered.columns if p in SORT_BY), ordered.columns[0])
    mean = ordered.mean().to_frame(MEAN).T
    rows = pd.concat([ordered, mean]).stack(future_stack=True).round(4).rename("percent")
    data = rows.rename_axis(["sample", "phylum"]).reset_index()
    data.insert(0, "panel", ["Average" if s == MEAN else "Samples" for s in data["sample"]])
    named = len([p for p in top.columns if p != OTHER and not _is_not_named(p)])
    rest = f"; the rest are summed as {OTHER} (grey)" if OTHER in top.columns else ""
    return Figure(
        name="composition_bars",
        title="Composition of each sample (phylum)",
        caption=f"Left: each bar is one sample, split into the share of each phylum (adding up to 100%: the "
                f"sample's whole total). Right: the average of all samples. The {named} most abundant phyla have "
                f"their own colour{rest}.{phylum_unnamed(top)} Samples are sorted by {lead} (most first)."
                f"{left_out} {descriptive}",
        details="Relative abundance from the uploaded table, summed to phylum (the p__ rank): each sample's "
                "whole total is 100%. Abundance with no phylum, or a phylum written as unassigned or uncultured, "
                f"is shown as {not_named('phylum')}. This is not the log-transformed table the tools analyse.",
        data=data, kind="composition", x="sample", y="percent", color_by="phylum",
        x_label="Samples", y_label="Relative abundance (%)", summary=(MEAN,),
        colors=tuple(colors.items()), suptitle=dataset,
        hover="%{x}<br>%{fullData.name}: %{y:.1f}%<extra></extra>", section="overview", level="phylum",
    )


def phylum_unnamed(top: pd.DataFrame) -> str:
    """The captions' sentence on "Not named at phylum" ('' if the table names a phylum everywhere)."""
    gap = not_named("phylum")
    if gap not in top:
        return ""
    share = float(top[gap].mean())
    return (f" {gap} (light grey) is the abundance your table does not name to phylum, on average "
            f"{'under 0.1' if share < 0.05 else f'{share:.1f}'}% of a sample; it is part of each sample's total.")


def split_violin_figure(percent: pd.DataFrame, groups: pd.Series, phyla: list[str],
                        colors: dict[str, str], descriptive: str = "", rank: str = "phylum") -> Figure | None:
    """Per phylum (or per taxon of ``rank``), the samples' % in group A (left half) and group B (right
    half); None unless two groups. Only the named taxa get a violin (not "Other" or "Not named")."""
    names = sorted(set(groups.loc[percent.index]))  # group A (left half) first, as in the other figures
    if len(names) != 2:
        return None
    long = percent[[p for p in phyla if p not in (OTHER, OTHER_NAMED) and not _is_not_named(p)]]
    long = long.assign(group=groups.loc[percent.index])
    data = (long.rename_axis("sample").reset_index()
            .melt(id_vars=["sample", "group"], var_name=rank, value_name="percent"))
    data = data.sort_values("group", key=lambda g: g.map({n: i for i, n in enumerate(names)}), kind="stable")
    data = data.reset_index(drop=True)
    data["percent"] = data["percent"].round(4)
    plural = PLURALS.get(rank, rank)
    return Figure(
        name="composition_violin" + ("" if rank == "phylum" else f"_{rank}"),
        title=f"{rank.capitalize()} abundance by group ({names[0]} vs {names[1]})",
        caption=f"For each of the most abundant {plural}, the spread of its share across samples: {names[0]} on "
                f"the left half, {names[1]} on the right; the line is the median. {descriptive}",
        level=rank,
        data=data, kind="violin", x=rank, y="percent", group_by="group",
        x_label="", y_label="Relative abundance (%)", colors=tuple(colors.items()),
        hover="%{x}<br>%{y:.1f}%<extra>%{fullData.name}</extra>", section="overview",
    )


# --- finer levels: class to species ------------------------------------------------------------
# The same figures at a finer rank, from the same uploaded table. The phylum figures above are
# unchanged; like them, abundance the table does not name at the rank is kept, as its own segment
# "Not named at <rank>" (small at phylum; at species it can be a third of a sample or more), so a bar stays an
# honest 100% of its sample. They differ in two ways:
# - The rest of the named taxa are summed as "Other (named, rare)": a second, darker grey.
# - Names are kept as the table writes them (GTDB suffixes included: at these ranks they can be
#   different taxa); taxa are counted by their whole lineage down to the rank, and where two share a
#   name, the family (for class, order and family: the rank above) is added to tell them apart.
FINE_RANKS = ("class", "order", "family", "genus", "species")
LINEAGE_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")  # as profiler.RANKS
PLURALS = {"phylum": "phyla", "class": "classes", "order": "orders", "family": "families", "genus": "genera",
           "species": "species"}
OTHER_NAMED = "Other (named, rare)"
NOT_NAMED_GREY = "#E6E6E6"  # lighter than OTHER_GREY; the legend names the two greys apart too
# Names that say a taxon has no name at that rank (16S databases write them out).
NO_NAME = {"", "unclassified", "unidentified", "uncultured", "unknown", "unassigned", "metagenome",
           "incertae_sedis", "none", "nan", "na"}
TOP_TAXA = TOP_PHYLA  # 12 named taxa, as at phylum
_SEP = ""  # joins a lineage's names into one grouping key
NOISE = 1e-6  # % of a sample: what floating-point sums leave over, not a share


def _same_phylum(name: str) -> str:
    """One name for a phylum's old and new names (``KNOWN_PHYLA``: the current one); others unchanged."""
    slot = _SLOT.get(name.lower())
    return KNOWN_PHYLA[slot][-1] if slot is not None else name


def not_named(rank: str) -> str:
    return f"Not named at {rank}"


def _is_not_named(name) -> bool:
    return str(name).startswith("Not named at ")


def _name_at(ranks: list[str], index: int) -> str:
    """The taxon's name at ``LINEAGE_RANKS[index]``: its ``x__`` rank (or, with no prefixes, the
    index-th part); "" when the lineage does not name it there."""
    prefix = f"{LINEAGE_RANKS[index][0]}__"
    if any(_PREFIX.match(r) for r in ranks):
        found = next((r for r in ranks if r.startswith(prefix)), "")
        name = found[len(prefix):]
    else:
        name = ranks[index] if index < len(ranks) else ""
    name = name.strip()
    lower = name.lower()
    return "" if lower in NO_NAME or lower.endswith("_unclassified") else name


def rank_percent(table: pd.DataFrame, rank: str) -> tuple[pd.DataFrame, dict]:
    """Samples × taxa of ``rank`` in % of each sample's whole total, with "Not named at <rank>" as a
    column (so every row adds up to 100), and what the captions need.

    ``table`` is samples × taxa, columns named by lineage (counts or relative abundances). A table that
    lists every rank (MetaPhlAn) uses its rows of that rank, and the total of its kingdom rows (else its
    phylum rows), so nothing is counted twice. Samples with a total of zero are left out.
    """
    if rank not in FINE_RANKS:
        raise ValueError(f"rank must be one of {FINE_RANKS}")
    by_lineage, totals = _rank_sums(table, rank)
    valid = totals > 0
    if not valid.any():
        raise CompositionError("No sample has any abundance.", "Check that the table holds counts or abundances.")
    percent = by_lineage.loc[valid].div(totals[valid], axis=0) * 100
    gap = 100 - percent.sum(axis=1)
    percent[not_named(rank)] = gap.where(gap > NOISE, 0.0)  # rounding left over is not abundance
    return percent, {"samples_dropped": int((~valid).sum()), "named_taxa": by_lineage.shape[1],
                     "not_named_percent": float(percent[not_named(rank)].mean())}


def _rank_sums(table: pd.DataFrame, rank: str) -> tuple[pd.DataFrame, pd.Series]:
    """Samples × named taxa of ``rank`` (``TAXON_RANKS``), summed by lineage in the table's own values, and
    each sample's whole total (``rank_percent`` says how). At phylum, a phylum's GTDB variants and its old
    and new names are one taxon, as in the phylum figures."""
    index = LINEAGE_RANKS.index(rank)
    values = table.apply(pd.to_numeric, errors="coerce")
    array = values.to_numpy(dtype=float)
    if not (pd.notna(array).all() and abs(array).max(initial=0) < float("inf")) or (array < 0).any():
        raise CompositionError(
            "Composition needs finite, non-negative abundances (counts or relative abundances).",
            "Upload the raw counts or relative abundances, not a log-transformed table.",
        )
    ranks = [_ranks(c) for c in values.columns]
    if _all_ranks(ranks):  # MetaPhlAn style: the rows of this rank, and the total of a coarser rank
        prefix = f"{rank[0]}__"
        keep = [i for i, r in enumerate(ranks) if r and r[-1].startswith(prefix)]
        top_rows = ([i for i, r in enumerate(ranks) if r and r[-1].startswith("k__")]
                    or [i for i, r in enumerate(ranks) if r and r[-1].startswith("p__")])
        totals = values.iloc[:, top_rows].sum(axis=1)
    else:
        keep = list(range(len(ranks)))
        totals = values.sum(axis=1)
    # An old and a new name of one phylum (Firmicutes, Bacillota: MetaPhlAn's database mixes them) are one
    # lineage, as in the phylum figures; otherwise its classes would be split in two.
    def phylum(i: int) -> str:
        name = _name_at(ranks[i], 1)
        return _same_phylum(gtdb_base(name) if rank == "phylum" and MERGE_GTDB_VARIANTS else name)

    lineages = {i: (phylum(i), *(_name_at(ranks[i], j) for j in range(2, index + 1))) for i in keep}
    named = [i for i in keep if lineages[i][-1]]
    if named:
        keys = [_SEP.join(lineages[i]) for i in named]  # one key per lineage (a list of tuples would not do)
        by_lineage = values.iloc[:, named].T.groupby(keys, sort=False).sum().T
        by_lineage.columns = _labels([tuple(k.split(_SEP)) for k in by_lineage.columns], rank)
    else:
        by_lineage = pd.DataFrame(index=values.index)
    return by_lineage, totals


# --- one taxon: the results page's "Look at one taxon" ------------------------------------------
TAXON_RANKS = ("phylum", *FINE_RANKS)


def taxon_choices(table: pd.DataFrame, rank: str) -> pd.Series:
    """Every named taxon of ``rank`` found in the table (the labels of ``rank_percent``), with its mean % of a
    sample's whole total, most abundant first. Abundance not named at the rank ("uncultured", "unclassified",
    …) is not a taxon, so it is not offered."""
    if rank not in TAXON_RANKS:
        raise ValueError(f"rank must be one of {TAXON_RANKS}")
    sums, totals = _taxon_sums(table, rank)
    valid = totals > 0
    percent = sums.loc[valid].div(totals[valid], axis=0) * 100
    mean = percent.mean()
    found = mean[(sums.loc[valid] > 0).any()]
    return found.sort_values(ascending=False, kind="stable").rename("mean_percent")


def taxon_abundance(table: pd.DataFrame, rank: str, taxon: str) -> pd.DataFrame | None:
    """Per sample (with a total above zero): ``percent`` (of the sample's whole total) and ``value`` (the
    table's own value: reads, for a table of counts) of one taxon (a label from ``taxon_choices``); None if the
    table has no such taxon."""
    if rank not in TAXON_RANKS:
        raise ValueError(f"rank must be one of {TAXON_RANKS}")
    sums, totals = _taxon_sums(table, rank)
    if taxon not in sums.columns:
        return None
    valid = totals > 0
    value = sums.loc[valid, taxon]
    return pd.DataFrame({"percent": value / totals[valid] * 100, "value": value})


def _taxon_sums(table: pd.DataFrame, rank: str) -> tuple[pd.DataFrame, pd.Series]:
    """The named taxa of ``rank`` and the whole totals; at phylum, exactly the phylum figures' phyla (their
    names and merges; "Other", the phylum-rank unclassified, is not a taxon)."""
    if rank != "phylum":
        return _rank_sums(table, rank)
    phylum, totals, _ = _phylum_sums(table, MERGE_GTDB_VARIANTS)
    return phylum.drop(columns=[OTHER], errors="ignore"), totals


def _labels(lineages: list[tuple[str, ...]], rank: str) -> list[str]:
    """Each lineage's name at the rank; where names clash, the family (genus, species) or the rank above
    (class, order, family) is added, and the whole lineage if that still clashes."""
    names = [lin[-1] for lin in lineages]
    counts = pd.Series(names).value_counts()
    # lineages hold phylum..rank, so position k is LINEAGE_RANKS[k + 1]
    parent = (LINEAGE_RANKS.index("family") if rank in ("genus", "species") else LINEAGE_RANKS.index(rank) - 1) - 1
    labels = []
    for lin, name in zip(lineages, names):
        if counts[name] == 1:
            labels.append(name)
        else:
            above = lin[parent] if 0 <= parent < len(lin) - 1 else ""
            labels.append(f"{name} ({above or LINEAGE_RANKS[parent + 1] + ' not named'})")
    seen = pd.Series(labels).value_counts()
    if (seen > 1).any():  # still clashing: the whole lineage
        labels = [lab if seen[lab] == 1 else " · ".join(p or "?" for p in lin) for lab, lin in zip(labels, lineages)]
    return labels


def top_taxa(percent: pd.DataFrame, rank: str, top_n: int = TOP_TAXA) -> pd.DataFrame:
    """The ``top_n`` named taxa by mean (most abundant first), then "Other (named, rare)" (the other named
    taxa summed) and "Not named at <rank>" (if any sample has some)."""
    gap = not_named(rank)
    named = [t for t in percent.columns if t != gap]
    order = list(percent[named].mean().sort_values(ascending=False, kind="stable").index)[:top_n]
    out = percent[order].copy()
    rest = [t for t in named if t not in order]
    if rest:
        out[OTHER_NAMED] = percent[rest].sum(axis=1)
    if gap in percent and (percent[gap] > 0).any():
        out[gap] = percent[gap]
    return out


def fine_colors(taxa, palette: str = PALETTE_DEFAULT) -> dict[str, str]:
    """A colour per named taxon in the given order (there are no fixed slots below phylum); the lab palette
    takes tab20's strong hues first, then its light ones, so neighbours in a stack differ. "Other (named,
    rare)" is grey and "Not named at <rank>" a lighter grey."""
    if palette not in PALETTES:
        raise ValueError(f"palette must be one of {sorted(PALETTES)}")
    colors = PALETTES[palette]
    if palette == "lab":
        tab20 = [c for c in map(to_hex, colormaps["tab20"].colors) if _not_grey(c)]
        colors = tuple(dict.fromkeys(tab20[0::2] + tab20[1::2] + list(colors)))
    out: dict[str, str] = {}
    shown = 0
    for t in dict.fromkeys(taxa):
        if t == OTHER_NAMED:
            out[t] = OTHER_GREY
        elif _is_not_named(t):
            out[t] = NOT_NAMED_GREY
        else:
            out[t] = colors[shown % len(colors)]
            shown += 1
    return out


def fine_composition_figure(top: pd.DataFrame, rank: str, colors: dict[str, str], info: dict,
                            dataset: str | None = None, left_out: str = "", descriptive: str = "") -> Figure:
    """The bars of ``composition_figure`` at a finer ``rank`` (``top`` from ``top_taxa``)."""
    named = [t for t in top.columns if t != OTHER_NAMED and not _is_not_named(t)]
    ordered = top.sort_values(named[0], ascending=False, kind="stable") if named else top
    mean = ordered.mean().to_frame(MEAN).T
    rows = pd.concat([ordered, mean]).stack(future_stack=True).round(4).rename("percent")
    data = rows.rename_axis(["sample", rank]).reset_index()
    data.insert(0, "panel", ["Average" if s == MEAN else "Samples" for s in data["sample"]])
    return Figure(
        name=f"composition_bars_{rank}",
        title=f"Composition of each sample ({rank})",
        caption=f"Left: each bar is one sample, split into the share of each {rank} (adding up to 100%: the "
                f"sample's whole total). Right: the average of all samples. {fine_greys(top, rank, info)}"
                + (f" Samples are sorted by {named[0]} (most first)." if named else "")
                + f"{left_out} {descriptive}",
        details=_fine_details(rank),
        data=data, kind="composition", x="sample", y="percent", color_by=rank,
        x_label="Samples", y_label="Relative abundance (%)", summary=(MEAN,),
        colors=tuple(colors.items()), suptitle=dataset,
        hover="%{x}<br>%{fullData.name}: %{y:.1f}%<extra></extra>", section="overview", level=rank,
    )


def fine_greys(top: pd.DataFrame, rank: str, info: dict) -> str:
    """The captions' account of the coloured taxa and the two greys."""
    plural = PLURALS[rank]
    named = [t for t in top.columns if t != OTHER_NAMED and not _is_not_named(t)]
    cover = top[named].sum(axis=1).mean() if named else 0.0
    text = (f"The {len(named)} most abundant {plural} have their own colour; together they make up on average "
            f"{cover:.0f}% of a sample.")
    rest = info["named_taxa"] - len(named)
    greys = []
    if OTHER_NAMED in top:
        greys.append(f"{OTHER_NAMED} (grey) sums the other {rest:,} named {plural if rest != 1 else rank}")
    if not_named(rank) in top:
        share = info["not_named_percent"]
        greys.append(f"{not_named(rank)} (light grey) is the abundance your table does not name down to "
                     f"{rank}, on average {'under 0.1' if share < 0.05 else f'{share:.1f}'}% of a sample")
    if greys:
        text += (" " + "; ".join(greys) + ". At fine levels the greys can be large: many rare taxa (a long tail) "
                 "and taxa named only to a coarser rank. They are part of each sample's total, not missing data.")
    return text


def _fine_details(rank: str) -> str:
    plural = PLURALS[rank]
    parent = "family" if rank in ("genus", "species") else LINEAGE_RANKS[LINEAGE_RANKS.index(rank) - 1]
    return (f"Relative abundance from the uploaded table: each sample's whole total is 100%. Taxa are summed by "
            f"their lineage down to {rank} (the {rank[0]}__ rank); names are as in the table (GTDB suffixes "
            f"included), and where two {plural} share a name the {parent} is added. Abundance named only to a "
            f"coarser rank, or as unclassified or uncultured, is shown as {not_named(rank)}. This is not the "
            "log-transformed table the tools analyse.")


def _ranks(column) -> list[str]:
    return [p for p in _RANKS.split(str(column).strip()) if p]


def _all_ranks(ranks: list[list[str]]) -> bool:
    """True for a table that lists every rank, as MetaPhlAn does: every clade's parent is a column too.

    (A 16S table where some taxa are named only down to their phylum is not one: their
    deeper relatives have no column per parent.)
    """
    seen = {tuple(r) for r in ranks}
    deeper = [r for r in ranks if len(r) > 1]
    return bool(deeper) and all(tuple(r[:-1]) in seen for r in deeper)


def _phylum(ranks: list[str]) -> str:
    """The phylum's name: the ``p__`` rank, else the second rank (after the kingdom); "Unassigned" if none."""
    tagged = [r for r in ranks if r.startswith("p__")]
    if tagged:
        name = _PREFIX.sub("", tagged[0])
    elif len(ranks) > 1 and not any(_PREFIX.match(r) for r in ranks):
        name = ranks[1]
    elif len(ranks) == 1:
        name = _PREFIX.sub("", ranks[0])  # a lone name: a phylum, or e.g. "Bacteria" (dropped)
    else:
        name = ""
    return name.strip() or "Unassigned"
