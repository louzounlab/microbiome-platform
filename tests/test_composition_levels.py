"""Composition at finer ranks (class to species): the same uploaded table, each bar still 100% of its sample."""

import io
import re

import numpy as np
import pandas as pd
import pytest

from microplatform import workflows
from microplatform.core.profiler import count_matrix, profile_labels, read_table
from microplatform.examples import EXAMPLE_COUNTS, EXAMPLE_TAGS
from microplatform.reports import composition, overview
from microplatform.reports.composition import (
    NOT_NAMED_GREY,
    OTHER_GREY,
    OTHER_NAMED,
    fine_colors,
    not_named,
    rank_percent,
    top_taxa,
)
from microplatform.reports.overview import FIGURE_NAMES
from microplatform.web import create_app

FAMILY_A = "k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae"
FAMILY_B = "k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae"


def _table(columns: dict[str, list[float]]) -> pd.DataFrame:
    return pd.DataFrame(columns, index=[f"s{i}" for i in range(len(next(iter(columns.values()))))])


def test_abundance_not_named_at_the_rank_is_kept_as_its_own_segment():
    table = _table({
        f"{FAMILY_A};g__Clostridium;s__": [6, 2],
        f"{FAMILY_B};g__": [2, 0],  # named only to family
        f"{FAMILY_B};g__uncultured;s__": [1, 0],  # written out as unnamed
        "k__Bacteria": [1, 8],  # named only to kingdom
    })
    percent, info = rank_percent(table, "genus")
    assert list(percent.columns) == ["Clostridium", not_named("genus")]
    assert percent.sum(axis=1).round(9).eq(100).all()  # the sample's whole total, nothing dropped
    assert percent.loc["s0", not_named("genus")] == pytest.approx(40)
    assert percent.loc["s1", not_named("genus")] == pytest.approx(80)
    assert info["named_taxa"] == 1 and info["not_named_percent"] == pytest.approx(60)


def test_names_that_clash_get_their_family_and_are_counted_by_lineage():
    table = _table({f"{FAMILY_A};g__Clostridium": [1, 1], f"{FAMILY_B};g__Clostridium": [3, 1],
                    f"{FAMILY_B};g__Blautia": [1, 3]})
    percent, _ = rank_percent(table, "genus")
    assert set(percent.columns) == {"Clostridium (Clostridiaceae)", "Clostridium (Lachnospiraceae)", "Blautia",
                                    not_named("genus")}
    assert percent.loc["s0", "Clostridium (Lachnospiraceae)"] == pytest.approx(60)


def test_gtdb_names_are_kept_as_they_are_below_phylum():
    table = _table({f"{FAMILY_A};g__Klebsiella_724518": [1, 1], f"{FAMILY_A};g__Klebsiella": [1, 3]})
    percent, _ = rank_percent(table, "genus")
    assert {"Klebsiella_724518", "Klebsiella"} <= set(percent.columns)


def test_an_old_and_a_new_phylum_name_do_not_split_a_class():
    table = _table({"k__Bacteria;p__Firmicutes;c__Clostridia": [1, 1], "k__Bacteria;p__Bacillota;c__Clostridia": [1, 1],
                    "k__Bacteria;p__Bacteroidota;c__Bacteroidia": [2, 2]})
    percent, _ = rank_percent(table, "class")
    assert set(percent.columns) == {"Clostridia", "Bacteroidia", not_named("class")}
    assert percent.loc["s0", "Clostridia"] == pytest.approx(50)


def test_a_metaphlan_table_uses_the_rows_of_the_rank_and_counts_nothing_twice():
    table = _table({
        "k__Bacteria": [100, 100],
        "k__Bacteria|p__Bacillota": [100, 100],
        "k__Bacteria|p__Bacillota|c__Clostridia": [100, 100],
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales": [100, 100],
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae": [100, 100],
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Blautia": [70, 40],
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia": [20, 60],
    })
    percent, _ = rank_percent(table, "genus")
    assert percent.loc["s0", "Blautia"] == pytest.approx(70) and percent.loc["s1", "Roseburia"] == pytest.approx(60)
    assert percent.loc["s0", not_named("genus")] == pytest.approx(10)  # the kingdom's 100 minus the genera


def test_top_taxa_keeps_twelve_then_the_two_greys_and_adds_up_to_100():
    rng = np.random.default_rng(0)
    table = _table({f"{FAMILY_A};g__G{i}": list(rng.integers(1, 50, 4)) for i in range(20)}
                   | {f"{FAMILY_A};g__": [30, 30, 30, 30]})
    top = top_taxa(rank_percent(table, "genus")[0], "genus")
    assert list(top.columns[-2:]) == [OTHER_NAMED, not_named("genus")] and len(top.columns) == 14
    assert top.sum(axis=1).round(9).eq(100).all()
    colors = fine_colors(top.columns)
    assert colors[OTHER_NAMED] == OTHER_GREY and colors[not_named("genus")] == NOT_NAMED_GREY != OTHER_GREY
    assert len({colors[t] for t in top.columns[:12]}) == 12  # every named taxon its own colour


def test_composition_levels_draws_each_finer_rank_and_leaves_phylum_unchanged():
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    figs = overview.composition_levels(counts, dataset="example")
    assert [f.name for f in figs] == [n for rank in composition.FINE_RANKS
                                      for n in (f"composition_bars_{rank}", f"composition_pie_{rank}")]
    assert all(f.level and f.section == "overview" for f in figs)
    genus = next(f for f in figs if f.name == "composition_bars_genus")
    assert "12 most abundant genera" in genus.caption and "not missing data" in genus.caption
    assert f"{OTHER_NAMED} (grey)" in genus.caption and "Not named at genus (light grey)" in genus.caption
    shown = genus.data[genus.data["panel"] == "Samples"].groupby("sample")["percent"].sum()
    assert shown.round(2).eq(100).all()
    pie = next(f for f in figs if f.name == "composition_pie_genus")
    assert dict(pie.colors) == dict(genus.colors) and "Not named at genus" in pie.caption
    # The phylum figures are as they were.
    phylum = overview.overview_figures(counts, workflows.CANONICAL)
    assert tuple(f.name for f in phylum) == FIGURE_NAMES
    assert phylum[0].title == "Composition of each sample (phylum)" and phylum[0].level == "phylum"


def test_the_violins_follow_the_rank_when_asked():
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    labels = profile_labels(read_table(EXAMPLE_TAGS, index_col=0), list(counts.index))
    figs = overview.composition_levels(counts, labels.values, labels.kind, split_violin=True)
    violin = next(f for f in figs if f.name == "composition_violin_genus")
    assert violin.level == "genus" and OTHER_NAMED not in set(violin.data["genus"])
    assert not any(n.startswith("Not named") for n in violin.data["genus"])


def test_no_finer_composition_for_a_table_with_negative_values():
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    assert overview.composition_levels(counts - counts.to_numpy().max()) == []


# --- the results page -------------------------------------------------------------------------

@pytest.fixture(scope="module")
def page(tmp_path_factory):
    instance = tmp_path_factory.mktemp("instance")
    client = create_app({"TESTING": True, "JOBS": "inline", "RUNS_DIR": instance / "runs",
                         "CACHE_DIR": instance / "cache"}).test_client()
    data = {"files": [(io.BytesIO(EXAMPLE_COUNTS.read_bytes()), EXAMPLE_COUNTS.name),
                      (io.BytesIO(EXAMPLE_TAGS.read_bytes()), EXAMPLE_TAGS.name)]}
    url = client.post("/upload", data=data, content_type="multipart/form-data").headers["Location"]
    client.post(f"{url}/run", data={"goal": "preprocess"})
    return lambda query="": client.get(url + query).get_data(as_text=True)


def _overview_names(html: str) -> list[str]:
    section = html[html.index('id="overview-heading"'):]
    return re.findall(r'id="fig-([a-z_]+)"', section[:section.index("</section>")])


def test_the_page_offers_the_levels_with_phylum_shown(page):
    html = page()
    switch = html[html.index('class="level-switch"'):html.index("</nav>", html.index('class="level-switch"'))]
    assert '<strong aria-current="true">Phylum</strong>' in switch
    assert re.findall(r"composition=([a-z]+)#composition", switch) == list(composition.FINE_RANKS)
    assert tuple(_overview_names(html)) == FIGURE_NAMES  # as before


def test_choosing_genus_swaps_only_the_composition_figures(page):
    names = _overview_names(page("?composition=genus"))
    assert names == ["composition_bars_genus", "composition_pie_genus", *FIGURE_NAMES[2:]]
    assert '<strong aria-current="true">Genus</strong>' in page("?composition=genus")


@pytest.mark.parametrize("asked", ["kingdom", "strain", "", "PHYLUM"])
def test_a_level_that_was_not_drawn_shows_phylum(page, asked):
    assert tuple(_overview_names(page(f"?composition={asked}"))) == FIGURE_NAMES
