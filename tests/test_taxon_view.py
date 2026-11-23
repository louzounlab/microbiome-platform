"""Look at one taxon: one bar per sample (its % of the sample's whole total), a box per group, no test."""

import io
import re
import zipfile

import numpy as np
import pandas as pd
import pytest

from microplatform import runs
from microplatform.core.profiler import count_matrix, profile_labels, read_table
from microplatform.examples import EXAMPLE_COUNTS, EXAMPLE_TAGS
from microplatform.reports import composition, overview
from microplatform.reports.composition import not_named, taxon_abundance, taxon_choices
from microplatform.web import create_app

FAMILY_A = "k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae"
FAMILY_B = "k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae"


def _table(columns: dict[str, list[float]]) -> pd.DataFrame:
    return pd.DataFrame(columns, index=[f"s{i}" for i in range(len(next(iter(columns.values()))))])


def test_the_list_holds_named_taxa_only_most_abundant_first():
    table = _table({
        f"{FAMILY_A};g__Clostridium": [1, 1],
        f"{FAMILY_B};g__Blautia": [5, 3],
        f"{FAMILY_B};g__uncultured": [40, 40],  # not a taxon: "Not named at genus"
        f"{FAMILY_B};g__": [20, 20],
        f"{FAMILY_B};g__Roseburia": [0, 0],  # in no sample
        "k__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__CFGB1765": [2, 2],
    })
    choices = taxon_choices(table, "genus")
    assert list(choices.index) == ["Blautia", "CFGB1765", "Clostridium"]
    assert not any(n == not_named("genus") or "uncultured" in n for n in choices.index)


def test_clashing_names_get_their_family_as_in_the_composition():
    table = _table({f"{FAMILY_A};g__Clostridium": [1, 1], f"{FAMILY_B};g__Clostridium": [3, 1]})
    assert set(taxon_choices(table, "genus").index) == {"Clostridium (Clostridiaceae)", "Clostridium (Lachnospiraceae)"}
    found = taxon_abundance(table, "genus", "Clostridium (Lachnospiraceae)")
    assert found.loc["s0", "percent"] == pytest.approx(75) and found.loc["s0", "value"] == 3


def test_the_share_is_of_the_sample_whole_total_and_a_lineage_is_summed():
    table = _table({f"{FAMILY_B};g__Blautia;s__Blautia obeum": [2, 0], f"{FAMILY_B};g__Blautia;s__": [3, 0],
                    f"{FAMILY_A};g__Clostridium": [5, 10]})
    found = taxon_abundance(table, "genus", "Blautia")
    assert found["percent"].tolist() == pytest.approx([50, 0]) and found["value"].tolist() == [5, 0]
    assert taxon_abundance(table, "genus", "Nothing") is None


def test_at_phylum_variants_and_old_and_new_names_are_one_phylum():
    table = _table({"k__Bacteria;p__Firmicutes": [1, 1], "k__Bacteria;p__Bacillota_A_368345": [1, 1],
                    "k__Bacteria;p__Bacteroidota": [2, 2]})
    assert set(taxon_choices(table, "phylum").index) == {"Firmicutes", "Bacteroidota"}  # a tie: the bars' choice


@pytest.mark.parametrize("more, less", [("Patescibacteria", "TM7"), ("TM7", "Patescibacteria"),
                                        ("Firmicutes", "Bacillota"), ("Bacillota", "Firmicutes")])
def test_at_phylum_a_phylum_has_the_bars_name_and_the_bars_share(more, less):
    table = _table({f"k__Bacteria;p__{more}": [6, 4], f"k__Bacteria;p__{less}": [1, 1],
                    "k__Bacteria;p__Bacteroidota": [2, 4], "k__Bacteria;p__;c__": [1, 1],  # no phylum
                    "k__Bacteria;p__uncultured": [0, 2]})
    bars, _ = composition.phylum_percent(table)
    choices = taxon_choices(table, "phylum")
    assert list(choices.index) == [more, "Bacteroidota"]  # the name the table uses more, as on the bars
    for taxon in choices.index:
        one = composition.taxon_abundance(table, "phylum", taxon)
        assert np.allclose(one["percent"], bars.loc[one.index, taxon])  # % of the whole total, as on the bars
    assert one["percent"].tolist() == pytest.approx([20, 33.333333])


@pytest.fixture(scope="module")
def example():
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    labels = profile_labels(read_table(EXAMPLE_TAGS, index_col=0), list(counts.index))
    return counts, labels


def test_bars_per_sample_with_reads_and_a_box_per_group(example):
    counts, labels = example
    taxon = taxon_choices(counts, "genus").index[0]
    bars, boxes = overview.taxon_figures(counts, "genus", taxon, labels.values, labels.kind)
    assert bars.kind == "bar" and list(bars.data.columns) == ["sample", "group", "percent", "reads"]
    assert len(bars.data) == len(counts) and bars.data["percent"].is_monotonic_decreasing
    assert bars.caption.startswith("Found in ") and "Descriptive only" in bars.caption
    assert boxes.kind == "box" and set(boxes.data["group"]) == {"Tag 0", "Tag 1"}
    assert "not tested" in boxes.caption and boxes.stats is None and not boxes.brackets
    present = int((bars.data["reads"] > 0).sum())
    assert f"Found in {present} of {len(counts)} samples" in bars.caption


def test_without_labels_only_the_bars_and_relative_tables_have_no_reads(example):
    counts, _ = example
    relative = counts.div(counts.sum(axis=1), axis=0)
    taxon = taxon_choices(relative, "family").index[0]
    (bars,) = overview.taxon_figures(relative, "family", taxon)
    assert "abundance_in_table" in bars.data and "reads" not in bars.data
    assert overview.taxon_figures(counts, "genus", "No such genus") == []


# --- the results page -------------------------------------------------------------------------

@pytest.fixture(scope="module")
def site(tmp_path_factory):
    instance = tmp_path_factory.mktemp("instance")
    client = create_app({"TESTING": True, "JOBS": "inline", "RUNS_DIR": instance / "runs",
                         "CACHE_DIR": instance / "cache"}).test_client()
    data = {"files": [(io.BytesIO(EXAMPLE_COUNTS.read_bytes()), EXAMPLE_COUNTS.name),
                      (io.BytesIO(EXAMPLE_TAGS.read_bytes()), EXAMPLE_TAGS.name)]}
    url = client.post("/upload", data=data, content_type="multipart/form-data").headers["Location"]
    client.post(f"{url}/run", data={"goal": "preprocess"})
    run = runs.load(instance / "runs", url.rstrip("/").split("/")[-1])
    return client, url, run


def _figure_ids(html: str) -> list[str]:
    """The figures of the data overview, in page order."""
    return re.findall(r'id="fig-([A-Za-z0-9_]+)"', html[html.index('id="overview-heading"'):])  # the last section


def test_the_control_offers_the_levels_and_the_named_taxa(site):
    client, url, run = site
    html = client.get(url).get_data(as_text=True)
    assert 'id="taxon"' in html and "Look at one taxon" in html
    form = html[html.index('class="taxon-pick"'):html.index("</form>", html.index('class="taxon-pick"'))]
    assert re.findall(r'<option value="([a-z]+)"', form)[:6] == list(composition.TAXON_RANKS)
    assert 'value="phylum" selected' in form  # the composition rank shown
    assert "Firmicutes" in form and not_named("phylum") not in form  # the example's own name, as on its bars
    # Nothing drawn until a taxon is asked for, and the overview's figures are all still there.
    assert not run.taxa_dir.exists() or not any(p.is_dir() for p in run.taxa_dir.iterdir())
    assert [n for n in _figure_ids(html) if not n.startswith("taxon_")] == list(overview.FIGURE_NAMES)


def test_a_new_level_brings_its_list(site):
    client, url, run = site
    names = client.get(f"{url}/taxa.json?rank=genus").get_json()
    assert names[:2] == list(taxon_choices(count_matrix(run.counts()), "genus").index[:2])
    assert client.get(f"{url}/taxa.json?rank=strain").get_json() == []


def test_picking_a_taxon_draws_it_and_serves_every_format(site):
    client, url, run = site
    genus = client.get(f"{url}/taxa.json?rank=genus").get_json()[0]
    html = client.get(url, query_string={"taxon_level": "genus", "taxon": genus}).get_data(as_text=True)
    drawn = [n for n in _figure_ids(html) if n.startswith("taxon_genus_")]
    assert len(drawn) == 2 and drawn[1] == f"{drawn[0]}_groups"
    assert f"{genus} (genus) in each sample" in html and "Found in " in html
    links = re.findall(r'href="(/runs/[^"]+/taxa/[^"]+\?download=1)"', html)
    assert {link.split("/")[-1].split("?")[0].rsplit(".", 1)[1] for link in links} == {"png", "svg", "html", "csv"}
    for link in links:
        assert client.get(link).status_code == 200
    csv = client.get(next(link for link in links if ".csv" in link and "_groups" not in link)).get_data(as_text=True)
    assert csv.splitlines()[0] == "sample,group,percent,reads"
    # The run's results are unchanged: the overview, the manifest and the ZIP know nothing of it.
    assert not any(f.name.startswith("taxon_") for f in runs.load(run.path.parent, run.id) and
                   __import__("microplatform").workflows.describe(run).figures)
    with zipfile.ZipFile(io.BytesIO(client.get(f"{url}/bundle").data)) as zf:
        assert not any("taxon_" in n for n in zf.namelist())


def test_an_unknown_taxon_says_so_and_draws_nothing(site):
    client, url, run = site
    html = client.get(url, query_string={"taxon_level": "genus", "taxon": "uncultured"}).get_data(as_text=True)
    assert "is not a genus in your table" in html and 'role="alert"' in html
    assert not [n for n in _figure_ids(html) if n.startswith("taxon_genus_unc")]


def test_only_the_latest_taxa_are_kept(site, monkeypatch):
    client, url, run = site
    from microplatform import workflows

    monkeypatch.setattr(workflows, "TAXA_KEPT", 2)
    for name in client.get(f"{url}/taxa.json?rank=phylum").get_json()[:3]:
        client.get(url, query_string={"taxon_level": "phylum", "taxon": name})
    assert len([p for p in run.taxa_dir.iterdir() if p.is_dir()]) == 2


def test_files_outside_a_drawn_taxon_are_not_served(site):
    client, url, _ = site
    assert client.get(f"{url}/taxa/..%2Fmeta.json/x").status_code == 404
    assert client.get(f"{url}/taxa/nothing/figures.json").status_code == 404
