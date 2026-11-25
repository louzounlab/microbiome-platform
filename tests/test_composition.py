"""Phylum composition in the lab's fig_1B style (reports/composition.py)."""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microplatform.reports.composition import (
    KNOWN_PHYLA,
    NOT_NAMED_GREY,
    OTHER_GREY,
    PALETTES,
    TOP_PHYLA,
    CompositionError,
    composition_figure,
    gtdb_base,
    merged_note,
    phylum_colors,
    phylum_percent,
    sort_samples,
    top_phyla,
)
from microplatform.reports.figures import OTHER


def _relative(columns: dict[str, list[float]]) -> pd.DataFrame:
    n = len(next(iter(columns.values())))
    return pd.DataFrame(columns, index=[f"s{i}" for i in range(n)], dtype=float)


# --- (a) tables with taxonomy strings -------------------------------------------------------------

def test_lineages_are_collapsed_to_phylum_in_percent_of_the_whole_total():
    table = _relative({
        "k__Bacteria;p__Firmicutes;c__Bacilli;o__;f__;g__;s__": [30, 10],
        "k__Bacteria; p__Firmicutes ;c__Clostridia;o__;f__;g__;s__": [30, 10],  # spaces around ';' are cleaned
        "k__Bacteria;p__Bacteroidetes;c__;o__;f__;g__;s__": [20, 70],
        "k__Bacteria;p__;c__;o__;f__;g__;s__": [20, 10],  # no phylum: its own segment, as at the finer ranks
    })
    percent, removed = phylum_percent(table)
    assert list(percent.columns) == ["Firmicutes", "Bacteroidetes", "Not named at phylum"]
    assert percent.loc["s0"].tolist() == pytest.approx([60, 20, 20]) and percent.loc["s1"].tolist() == pytest.approx([20, 70, 10])
    assert np.allclose(percent.sum(axis=1), 100)
    assert removed["unassigned_percent"] == pytest.approx(15) == removed["not_named_percent"]  # (20% + 10%) / 2


@pytest.mark.parametrize("dropped", ["Bacteria", "Archaea", "Unassigned", "Unclassified", "k__Bacteria", "unassigned",
                                     "k__Bacteria;p__uncultured", "k__Bacteria;p__unknown"])
def test_columns_without_a_phylum_are_not_named_at_phylum(dropped):
    percent, _ = phylum_percent(_relative({"p__Firmicutes": [50, 20], dropped: [50, 80]}))
    assert list(percent.columns) == ["Firmicutes", "Not named at phylum"]
    assert percent["Firmicutes"].tolist() == [50, 20] and percent["Not named at phylum"].tolist() == [50, 80]


def test_a_table_that_names_every_phylum_has_no_grey_segment_and_keeps_its_shares():
    percent, removed = phylum_percent(_relative({"p__Firmicutes": [1, 3], "p__Bacteroidetes": [3, 1]}))
    assert list(percent.columns) == ["Firmicutes", "Bacteroidetes"] and removed["not_named_percent"] == 0
    assert percent["Firmicutes"].tolist() == [25, 75]


def test_not_named_at_phylum_comes_last_in_light_grey_and_is_never_other():
    cols = {f"p__P{i:02d}": [i + 1, i + 1] for i in range(14)} | {"k__Bacteria;p__": [200, 200]}
    percent, _ = phylum_percent(_relative(cols))
    top = top_phyla(percent)
    assert top.columns[-1] == "Not named at phylum" and OTHER not in top  # the big grey is not ranked as a phylum
    assert len(top.columns) == 15  # every phylum, then the light grey
    assert np.allclose(top.sum(axis=1), 100)
    colors = phylum_colors(top.columns)
    assert colors["Not named at phylum"] == NOT_NAMED_GREY and len({colors[f"P{i:02d}"] for i in range(14)}) == 14


def test_beyond_max_taxa_the_rest_is_other():
    table = _relative({f"p__P{i:02d}": [i + 1, i + 1] for i in range(20)})
    top = top_phyla(phylum_percent(table)[0], top_n=12)
    assert len(top.columns) == 13 and top.columns[-1] == OTHER
    assert np.allclose(top.sum(axis=1), 100)


def test_the_bars_caption_names_the_light_grey_and_leaves_nothing_out():
    percent, _ = phylum_percent(_relative({"p__Firmicutes": [90, 80], "k__Bacteria;p__": [10, 20]}))
    top = top_phyla(percent)
    caption = composition_figure(top, phylum_colors(top.columns)).caption
    assert "the sample's whole total" in caption and "left out" not in caption
    assert "Not named at phylum (light grey)" in caption and "on average 15.0% of a sample" in caption


def test_phylum_level_tables_with_or_without_prefixes():
    a, _ = phylum_percent(_relative({"Bacteria;Firmicutes": [1, 3], "Bacteria;Bacteroidetes": [3, 1]}))
    b, _ = phylum_percent(_relative({"p__Firmicutes": [1, 3], "p__Bacteroidetes": [3, 1]}))
    pd.testing.assert_frame_equal(a, b)
    assert a["Firmicutes"].tolist() == [25, 75]


def test_a_16s_table_with_taxa_named_only_to_phylum_is_not_taken_for_metaphlan():
    table = _relative({"k__Bacteria;p__Firmicutes": [10, 10],  # an ASV classified only to phylum
                       "k__Bacteria;p__Firmicutes;c__Bacilli": [10, 0],
                       "k__Bacteria;p__Bacteroidetes;c__Bacteroidia": [20, 10]})
    percent, _ = phylum_percent(table)
    assert percent.loc["s0", "Firmicutes"] == pytest.approx(50)  # both Firmicutes columns are counted


# --- (b) MetaPhlAn tables with every rank ---------------------------------------------------------

@pytest.mark.parametrize("sep", ["|", ","])
def test_metaphlan_uses_only_the_phylum_rows(sep):
    j = sep.join
    table = _relative({
        "k__Bacteria": [100, 100],
        j(["k__Bacteria", "p__Firmicutes"]): [60, 30],
        j(["k__Bacteria", "p__Firmicutes", "c__Bacilli"]): [60, 30],  # not counted again
        j(["k__Bacteria", "p__Bacteroidetes"]): [40, 70],
        j(["k__Bacteria", "p__Bacteroidetes", "c__Bacteroidia"]): [40, 70],
    })
    percent, _ = phylum_percent(table)
    assert percent.loc["s0"].tolist() == [60, 40] and percent.loc["s1"].tolist() == [30, 70]


# --- one phylum, one bar (display only) ------------------------------------------------------------

@pytest.mark.parametrize("variant, base", [
    ("Bacillota_A_368345", "Bacillota"), ("Bacillota_C", "Bacillota"), ("Desulfobacterota_G_459543", "Desulfobacterota"),
    ("Methanobacteriota_A_1229", "Methanobacteriota"), ("Bacteroidota", "Bacteroidota"),
    ("Candidatus_Saccharibacteria", "Candidatus_Saccharibacteria"), ("SAR324", "SAR324"),
])
def test_gtdb_base(variant, base):
    assert gtdb_base(variant) == base


def test_gtdb_variants_are_merged_for_display_and_the_note_says_so():
    table = _relative({"k__Bacteria;p__Bacillota_A_368345": [30, 10], "k__Bacteria;p__Bacillota_C": [10, 10],
                       "k__Bacteria;p__Bacteroidota": [60, 80]})
    percent, removed = phylum_percent(table)
    assert list(percent.columns) == ["Bacillota", "Bacteroidota"]
    assert percent.loc["s0", "Bacillota"] == pytest.approx(40)
    assert removed["merged"]["variants"] == {"Bacillota_A_368345": "Bacillota", "Bacillota_C": "Bacillota"}
    assert "GTDB variants (e.g. Bacillota_A_368345 → Bacillota)" in merged_note(removed["merged"])
    kept, removed = phylum_percent(table, merge_gtdb_variants=False)  # easy to switch off
    assert "Bacillota_A_368345" in kept.columns and merged_note(removed["merged"]) == ""


def test_old_and_new_names_of_one_phylum_become_the_one_the_table_uses_more():
    table = _relative({"k__Bacteria;p__Firmicutes": [50, 50], "k__Bacteria;p__Bacillota": [1, 1],
                       "k__Bacteria;p__Bacteroidota": [49, 49]})
    percent, removed = phylum_percent(table)
    assert list(percent.columns) == ["Firmicutes", "Bacteroidota"]
    assert removed["merged"]["synonyms"] == {"Bacillota": "Firmicutes"}
    assert "old and new names (e.g. Bacillota → Firmicutes)" in merged_note(removed["merged"])


def test_unclassified_at_phylum_rank_is_counted_in_other_not_left_out():
    percent, _ = phylum_percent(_relative({"k__Bacteria;p__Bacteria_unclassified": [1, 1],
                                           "k__Bacteria;p__Bacteroidota": [9, 9]}))
    assert sorted(percent.columns) == sorted(["Bacteroidota", OTHER])
    assert percent.loc["s0", OTHER] == pytest.approx(10) and percent.loc["s0"].sum() == pytest.approx(100)
    top = top_phyla(percent)
    assert list(top.columns) == ["Bacteroidota", OTHER]  # "Other" stays last and grey


# --- guards ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [-1.0, np.inf, np.nan])
def test_values_must_be_finite_and_not_negative(bad):
    with pytest.raises(CompositionError):
        phylum_percent(_relative({"p__Firmicutes": [5, bad], "p__Bacteroidetes": [5, 5]}))


def test_samples_with_no_abundance_are_dropped_and_unnamed_ones_kept():
    percent, removed = phylum_percent(_relative({"p__Firmicutes": [5, 0, 1, 0], "Unassigned": [5, 9, 0, 0]}))
    assert list(percent.index) == ["s0", "s1", "s2"] and removed["samples_dropped"] == 1
    assert percent.loc["s1", "Not named at phylum"] == 100  # reads, but none with a phylum: one grey bar
    with pytest.raises(CompositionError):
        phylum_percent(_relative({"Unassigned": [1, 2]}))


# --- top phyla, order and colours -----------------------------------------------------------------

def test_every_phylum_by_mean_with_its_own_colour_and_no_other():
    table = _relative({f"p__P{i:02d}": [i + 1, i + 1] for i in range(20)})
    top = top_phyla(phylum_percent(table)[0])
    assert len(top.columns) == 20 and OTHER not in top
    assert list(top.columns[:3]) == ["P19", "P18", "P17"]
    assert np.allclose(top.sum(axis=1), 100)
    colors = phylum_colors(top.columns)
    assert len({colors[p] for p in top.columns}) == 20


def test_a_phylum_absent_from_every_sample_is_not_drawn():
    table = _relative({"p__Firmicutes": [90, 80], "p__Bacteroidetes": [10, 20], "p__Ghost": [0, 0]})
    assert "Ghost" not in top_phyla(phylum_percent(table)[0])


def test_past_the_palette_the_colours_are_generated_distinct_and_never_grey():
    from microplatform.reports.composition import LEGEND_MAX, extra_colors, fine_colors
    names = [f"T{i:04d}" for i in range(1300)]
    colors = fine_colors(names)
    assert len(set(colors.values())) == 1300  # a thousand species: still one colour each
    for color in extra_colors(500):
        r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
        assert max(r, g, b) - min(r, g, b) > 16, color
    many = phylum_colors([f"P{i}" for i in range(60)])
    assert len({c for p, c in many.items() if p.startswith("P")}) == 60
    assert LEGEND_MAX == 20


def test_the_legend_is_dropped_above_twenty_taxa_and_the_caption_says_so(tmp_path):
    from microplatform.reports.figures import export, to_plotly
    few = _relative({f"p__P{i:02d}": [i + 1, i + 1] for i in range(20)})
    many = _relative({f"p__P{i:02d}": [i + 1, i + 1] for i in range(21)})
    for table, shown in ((few, True), (many, False)):
        top = top_phyla(phylum_percent(table)[0])
        fig = composition_figure(top, phylum_colors(top.columns))
        assert fig.legend is shown
        assert ("there is no legend" in fig.caption) is not shown
        assert any(t.showlegend for t in to_plotly(fig).data) is shown
        export([fig], tmp_path / str(shown))
        svg = (tmp_path / str(shown) / "composition_bars.svg").read_text(encoding="utf-8")
        assert (">P00<" in svg or "P00</text>" in svg or "P00" in re.sub(r"<[^>]+>", " ", svg)) is shown


def test_samples_are_sorted_by_firmicutes_else_the_top_phylum():
    top = pd.DataFrame({"Bacteroidetes": [60, 20, 40], "Firmicutes": [40, 80, 60]}, index=["a", "b", "c"])
    assert list(sort_samples(top).index) == ["b", "c", "a"]
    top = pd.DataFrame({"Bacteroidota": [60, 20, 40], "Bacillota": [40, 80, 60]}, index=["a", "b", "c"])
    assert list(sort_samples(top).index) == ["b", "c", "a"]  # the same phylum, under its new name
    top = pd.DataFrame({"Bacteroidetes": [60, 20, 40], "Other": [40, 80, 60]}, index=["a", "b", "c"])
    assert list(sort_samples(top).index) == ["a", "c", "b"]


@pytest.mark.parametrize("palette", sorted(PALETTES))
def test_a_phylum_keeps_its_colour_across_datasets(palette):
    one = phylum_colors(["Firmicutes", "Bacteroidetes", "Mystery", OTHER], palette)
    two = phylum_colors(["Bacteroidetes", "Proteobacteria", "Firmicutes"], palette)
    assert one["Firmicutes"] == two["Firmicutes"] and one["Bacteroidetes"] == two["Bacteroidetes"]
    assert one[OTHER] == OTHER_GREY and phylum_colors(["Bacillota"], palette)["Bacillota"] == one["Firmicutes"]


@pytest.mark.parametrize("palette", sorted(PALETTES))
def test_colours_are_distinct_within_a_figure_and_never_grey(palette):
    names = [n[0] for n in KNOWN_PHYLA][:6] + [f"Unknown{i}" for i in range(6)]  # 12 shown, as at most
    colors = phylum_colors(names + [OTHER], palette)
    shown = [colors[n] for n in names]
    assert len(set(shown)) == len(shown) and OTHER_GREY not in shown
    for color in shown:
        r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
        assert max(r, g, b) - min(r, g, b) > 16, color  # no grey that reads as "Other"


def test_the_default_is_the_lab_look():
    assert phylum_colors(["Firmicutes"])["Firmicutes"] == "#1f77b4"  # tab20's first colour, as in fig_1B
    with pytest.raises(ValueError):
        phylum_colors(["Firmicutes"], "rainbow")


# --- one map everywhere, on a real dataset with more than 12 phyla --------------------------------

GIMIC = Path(__file__).parent / "data" / "gimic_example"  # the authors' published GIMIC example: 16 phyla


@pytest.fixture(scope="module")
def many_phyla(tmp_path_factory):
    from microplatform.core.preprocessing import CANONICAL
    from microplatform.core.profiler import count_matrix, profile_labels, read_table
    from microplatform.reports import figures, overview

    counts = count_matrix(read_table(GIMIC / "for_preprocess.csv", low_memory=False))
    tags = read_table(GIMIC / "tag.csv")
    tags = tags.set_index(tags.columns[0])
    tags.index = tags.index.astype(str)
    labels = profile_labels(tags, list(counts.index))
    figs = {f.name: f for f in overview.overview_figures(counts, CANONICAL, labels.values, labels.kind,
                                                        split_violin=True)}
    out = tmp_path_factory.mktemp("many_phyla")
    figures.export(list(figs.values()), out)
    return counts, figs, out


def _legend_names(svg: str, names: set[str]) -> list[str]:
    """The phylum names in an SVG's text, in drawing order (the legend lists every one)."""
    texts = re.findall(r"<text[^>]*>(.*?)</text>", svg, re.S)
    return [t.strip() for t in texts if t.strip() in names]


def test_every_per_phylum_plot_uses_the_one_global_map_and_the_greys_are_last(many_phyla):
    from microplatform.reports.figures import to_plotly

    counts, figs, out = many_phyla
    percent, _ = phylum_percent(counts)
    assert percent.shape[1] > TOP_PHYLA  # the point of this dataset: more phyla than the old 12
    top = top_phyla(percent)
    gap = "Not named at phylum"  # this table has abundance with no phylum: the light grey, last
    expected = phylum_colors(top.columns)  # the global map, lab palette
    named = [p for p in expected if p not in (OTHER, gap)]
    assert OTHER not in top and len(named) == percent.shape[1] - 1  # every phylum has its own colour
    assert len({expected[p] for p in named}) == len(named) and expected[gap] == NOT_NAMED_GREY
    for name in ("composition_bars", "composition_pie"):
        fig = figs[name]
        assert {k: v for k, v in fig.colors if k != OTHER} == {k: v for k, v in expected.items() if k != OTHER}, name
        # The data: the light grey is the last segment of every bar and the last slice of every pie.
        by = "sample" if name == "composition_bars" else "group"
        assert set(fig.data.groupby(by, sort=False)["phylum"].last()) == {gap}, name
        assert OTHER not in set(fig.data["phylum"]), name
        # The drawing (SVG, same as the PNG): 16 phyla is within the legend's 20, the grey is its last entry.
        assert fig.legend, name
        svg = (out / f"{name}.svg").read_text(encoding="utf-8")
        assert _legend_names(svg, set(expected))[-1] == gap, name
        assert NOT_NAMED_GREY.lower() in svg.lower(), name
    bars = to_plotly(figs["composition_bars"])
    assert {t.name for t in bars.data[-2:]} == {gap}  # drawn last in both panels
    assert {t.name: t.marker.color for t in bars.data} == {k: v for k, v in expected.items() if k != OTHER}
    for pie in to_plotly(figs["composition_pie"]).data:
        assert list(pie.labels)[-1] == gap and list(pie.marker.colors)[-1] == NOT_NAMED_GREY
        assert dict(zip(pie.labels, pie.marker.colors)) == {p: expected[p] for p in pie.labels}
    # The split violins compare groups (coloured by group); they show the 12 most abundant named phyla, in order.
    violin = figs["composition_violin"]
    assert list(dict.fromkeys(violin.data["phylum"])) == named[:TOP_PHYLA]


def test_a_figure_without_explicit_colours_still_gets_the_global_map():
    from microplatform.reports.figures import Figure, _color_map

    data = pd.DataFrame({"phylum": ["Bacteroidetes", "Firmicutes", OTHER], "percent": [30.0, 60.0, 10.0]})
    fig = Figure(name="x", title="x", caption="x", data=data, kind="pie", x="phylum", y="percent",
                 x_label="", y_label="")
    assert _color_map(fig, "phylum") == phylum_colors(["Bacteroidetes", "Firmicutes", OTHER])
