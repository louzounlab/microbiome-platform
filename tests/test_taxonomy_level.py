"""The run's taxonomy level: one level per run, for every tool, chosen from the levels the table reaches."""

import html
import io
import re

import pytest

import microplatform as mp
from microplatform import workflows
from microplatform.core.inputs import InputError
from microplatform.core.preprocessing import CANONICAL, REPRODUCED_AT, SPECIES
from microplatform.engines.mimic import SPECIES_ONLY_PLOTS
from microplatform.examples import EXAMPLE_COUNTS, EXAMPLE_TAGS
from microplatform.web import create_app

GENUS, PHYLUM = 6, 2


@pytest.fixture
def run(tmp_path):
    return workflows.start_example_run(tmp_path / "runs")


def test_the_levels_offered_are_those_the_table_reaches_with_species_the_default(run):
    view = workflows.describe(run)
    assert view.taxonomy_level == SPECIES == CANONICAL.taxonomy_level  # the default is unchanged
    # The example is named down to genus for most taxa: phylum to genus, and species (the default).
    assert [(o.level, o.name) for o in view.taxonomy_levels] == [
        (2, "phylum"), (3, "class"), (4, "order"), (5, "family"), (6, "genus"), (7, "species")]
    hints = {o.name: o.hint for o in view.taxonomy_levels}
    assert hints["phylum"] == "few features: less informative for iMic, and GIMIC cannot run at this level"
    assert hints["class"] == "few features: less informative for iMic and GIMIC"
    assert hints["order"] is None and hints["genus"] is None
    assert hints["species"].startswith("default; your table names a species for")


def test_a_level_can_be_given_by_name_or_number(run):
    workflows.set_taxonomy_level(run, "Genus")
    assert workflows.describe(run).taxonomy_level == GENUS
    workflows.set_taxonomy_level(run, "5")
    assert workflows.describe(run).taxonomy_level == 5


@pytest.mark.parametrize("level", ["kingdom", 1, 8, "strain", ""])
def test_levels_the_table_cannot_use_are_refused_with_the_choices(run, level):
    with pytest.raises(InputError) as caught:
        workflows.set_taxonomy_level(run, level)
    assert "phylum, class, order, family, genus, species" in caught.value.fix


def test_a_table_used_as_it_is_has_no_level(run):
    workflows.set_preprocessing(run, "as_is")
    with pytest.raises(InputError, match="not grouped at any taxonomy level"):
        workflows.set_taxonomy_level(run, "genus")


def test_the_run_level_reaches_the_tool_and_the_record(run, tmp_path):
    workflows.set_taxonomy_level(run, "genus")
    workflows.run_goal(run, "preprocess", tmp_path / "cache")
    view = workflows.describe(run)
    assert view.record["settings"]["taxonomy_level"] == GENUS
    assert "at genus level" in view.record["summary"][0]
    assert view.representation_note is None  # the run's own level, not a "different table"


def test_by_default_every_tool_still_uses_species(run, tmp_path):
    workflows.run_goal(run, "preprocess", tmp_path / "cache")
    assert workflows.describe(run).record["settings"]["taxonomy_level"] == SPECIES


def test_locate_keeps_its_z_score_at_the_run_level(run, tmp_path):
    from microplatform.engines import ENGINES

    settings = ENGINES["metabolites"].settings_for({}).with_changes(taxonomy_level=GENUS)
    dataset, _ = workflows.build_dataset(run, settings, tmp_path / "cache", "LOCATE")
    assert dataset.settings.taxonomy_level == GENUS and dataset.settings.z_scoring == "col"


def test_mimic_below_species_runs_its_test_and_skips_the_authors_plots(run, tmp_path):
    workflows.set_taxonomy_level(run, "genus")
    workflows.run_goal(run, "differ", tmp_path / "cache")
    view = workflows.describe(run)
    assert view.status == "done"
    record = view.engine["record"]
    assert record["authors_plots"] == record["rp_sp_check"] == "not drawn (species level only)"
    assert view.engine["figure_note"] == SPECIES_ONLY_PLOTS
    shown = [f.name for f in view.findings]
    assert "mimic_cladogram" not in shown and "mimic_rp_sp" not in shown
    assert {"significant_taxa", "leaf_test_only", "apriori_test"} <= set(view.engine["tables"])


@pytest.mark.parametrize(("goal", "level"), [("group", "class"), ("predict", "phylum")])
def test_image_tools_say_a_coarse_level_holds_little(run, tmp_path, goal, level):
    workflows.set_taxonomy_level(run, level)
    workflows.run_goal(run, goal, tmp_path / "cache")
    notes = " ".join(workflows.describe(run).engine["notes"])
    assert f"At {level} level there are few features" in notes and "less informative" in notes


def test_gimic_at_phylum_says_why_it_cannot_measure_distances(run, tmp_path):
    """The authors' distances are all undefined at phylum (images of 3 rows): a clear message, no crash."""
    from microplatform.core.errors import UserFacingError

    workflows.set_taxonomy_level(run, "phylum")
    with pytest.raises(UserFacingError, match="could not measure the distances at phylum level") as caught:
        workflows.run_goal(run, "group", tmp_path / "cache")
    assert "Choose class or a finer taxonomy level" in caught.value.fix
    assert workflows.describe(run).error["message"].startswith("GIMIC could not measure")


def test_no_coarse_note_at_genus(run, tmp_path):
    workflows.set_taxonomy_level(run, "genus")
    workflows.run_goal(run, "group", tmp_path / "cache")
    assert "few features" not in " ".join(workflows.describe(run).engine["notes"])


def test_analyze_takes_a_taxonomy_level(tmp_path):
    view = mp.analyze(EXAMPLE_COUNTS, goal="preprocess", out_dir=tmp_path / "runs", taxonomy_level="family")
    assert view.taxonomy_level == 5 and view.record["settings"]["taxonomy_level"] == 5


# --- the website ------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def client(tmp_path_factory):
    instance = tmp_path_factory.mktemp("instance")
    app = create_app({"TESTING": True, "JOBS": "inline", "RUNS_DIR": instance / "runs",
                      "CACHE_DIR": instance / "cache"})
    return app.test_client()


def _upload(client) -> str:
    files = [(EXAMPLE_COUNTS.name, EXAMPLE_COUNTS.read_bytes()), (EXAMPLE_TAGS.name, EXAMPLE_TAGS.read_bytes())]
    data = {"files": [(io.BytesIO(content), name) for name, content in files]}
    return client.post("/upload", data=data, content_type="multipart/form-data").headers["Location"]


def test_the_plan_offers_the_levels_with_species_selected(client):
    page = html.unescape(client.get(_upload(client)).get_data(as_text=True))
    chooser = re.search(r'<select name="taxonomy_level".*?</select>', page, re.S).group(0)
    values = re.findall(r'<option value="(\d)"', chooser)
    assert values == ["2", "3", "4", "5", "6", "7"]
    assert re.search(r'value="7"\s+selected', chooser) and "Species (default" in chooser
    assert "Phylum (few features: less informative for iMic, and GIMIC cannot run at this level)" in chooser


def test_the_level_is_sent_with_run_and_shown_in_the_record(client):
    url = _upload(client)
    page = html.unescape(client.post(f"{url}/run", data={"goal": "preprocess", "taxonomy_level": "6"},
                                     follow_redirects=True).get_data(as_text=True))
    assert "at genus level" in page
    record = page[page.index("Reproducibility record"):]
    assert f"Reproduction against the authors' published results</dt><dd>{REPRODUCED_AT}" in record
    assert page.count("Reproduction against the authors") == 1  # only in the folded record
    assert "Different table" not in page


def test_a_level_the_table_does_not_reach_is_refused_on_the_page(client):
    url = _upload(client)
    page = html.unescape(client.post(f"{url}/run", data={"goal": "preprocess", "taxonomy_level": "1"},
                                     follow_redirects=True).get_data(as_text=True))
    assert "is not a taxonomy level this table can use" in page


def test_mimic_below_species_shows_the_line_where_its_plots_would_be(client):
    url = _upload(client)
    page = html.unescape(client.post(f"{url}/run", data={"goal": "differ", "taxonomy_level": "6"},
                                     follow_redirects=True).get_data(as_text=True))
    findings = page[page.index("<h3>Findings</h3>") + len("<h3>Findings</h3>"):]
    assert findings.lstrip().startswith('<p class="muted figure-note">')  # first, where the plots would be
    assert SPECIES_ONLY_PLOTS in findings
    assert "mimic_cladogram" not in page


def test_mimic_at_species_has_no_such_line(client):
    url = _upload(client)
    page = html.unescape(client.post(f"{url}/run", data={"goal": "differ"}, follow_redirects=True)
                         .get_data(as_text=True))
    assert SPECIES_ONLY_PLOTS not in page and "mimic_cladogram" in page


def test_methods_page_names_the_default_level(client):
    page = client.get("/methods").get_data(as_text=True)
    assert "species (default; each run can choose another level)" in page
