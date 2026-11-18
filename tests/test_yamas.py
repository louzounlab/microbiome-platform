"""YaMAS output as an upload (core/yamas.py): otu.csv + taxonomy.csv become one counts table, counts kept."""

import ast
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microplatform import workflows
from microplatform.core import yamas
from microplatform.core.errors import UserFacingError
from microplatform.core.inputs import MAX_FILES, InputError, sort_uploads
from microplatform.core.preprocessing import preprocess
from microplatform.core.profiler import count_matrix, profile_counts, read_table
from microplatform.examples import EXAMPLE_COUNTS, EXAMPLE_TAGS
from microplatform.web import create_app

BRIDGE = Path(__file__).parents[1] / "docs" / "preprocess_16S.py"  # the lab's bridge script
TAGS = (EXAMPLE_TAGS.name, EXAMPLE_TAGS.read_bytes())


def _messy(lineage: str, i: int) -> str:
    """A lineage as a classifier writes it: d__ for the kingdom, spaces, and a confidence on each named rank."""
    parts = lineage.replace("k__", "d__", 1).split(";")
    return "; ".join(p + (f"({0.5 + (i % 5) / 10:.2f})" if not p.endswith("__") else "") for p in parts)


def _example_as_yamas() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The bundled example (40 samples × 263 taxa) as YaMAS would write it: OTUs × samples, and a taxonomy."""
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    ids = [f"otu{i:04d}" for i in range(counts.shape[1])]
    otu = pd.DataFrame(counts.to_numpy().T.astype(int), index=pd.Index(ids, name="OTU"), columns=counts.index)
    taxonomy = pd.DataFrame({"Taxon": [_messy(c, i) for i, c in enumerate(counts.columns)], "Confidence": 0.9},
                            index=pd.Index(ids, name="Feature ID"))
    return otu, taxonomy


def _files(otu: pd.DataFrame, taxonomy: pd.DataFrame, names=("otu.csv", "taxonomy.csv")):
    return [(names[0], otu.to_csv().encode()), (names[1], taxonomy.to_csv().encode())]


OTU, TAXONOMY = _example_as_yamas()


# --- the lab's standardisation ---------------------------------------------------------------------

def test_normalize_taxonomy_name_is_the_labs_own_unchanged():
    if not BRIDGE.exists():
        pytest.skip("the lab's bridge script docs/preprocess_16S.py is not in this checkout")
    source = BRIDGE.read_text(encoding="utf-8")
    theirs = {n.name: ast.dump(n) for n in ast.parse(source).body if isinstance(n, ast.FunctionDef)}
    ours = {n.name: ast.dump(n) for n in ast.parse(Path(yamas.__file__).read_text(encoding="utf-8")).body
            if isinstance(n, ast.FunctionDef)}
    assert ours["normalize_taxonomy_name"] == theirs["normalize_taxonomy_name"]
    assert (yamas.RANKS, yamas.MISSING_TAXON_VALUES) == (("k", "p", "c", "o", "f", "g", "s"),
                                                          {"", "nan", "none", "unassigned", "unknown"})


def test_the_labs_names_use_the_contracts_rank_separator():
    """The copied function writes ';' itself; the shared format must use the same separator."""
    from microplatform.core.contract import RANK_SEPARATOR

    assert yamas.normalize_taxonomy_name("d__Bacteria;p__Firmicutes").split(RANK_SEPARATOR) == [
        "k__Bacteria", "p__Firmicutes"]


@pytest.mark.parametrize("raw, standard", [
    ("D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli", "k__Bacteria;p__Firmicutes;c__Bacilli"),
    ("d__Bacteria(0.99); p__Firmicutes(0.71); c__Bacilli", "k__Bacteria;p__Firmicutes;c__Bacilli"),
    ("k__Bacteria; c__Clostridia", "k__Bacteria;p__;c__Clostridia"),  # a missing rank keeps its place
    ("Unassigned", "Unassigned"), ("", "Unassigned"),
])
def test_lineages_are_written_in_the_standard_form(raw, standard):
    assert yamas.normalize_taxonomy_name(raw) == standard


# --- combining the two files -----------------------------------------------------------------------

def test_the_pair_becomes_a_counts_table_with_the_read_counts_kept():
    table, details = yamas.counts_from_16s(OTU.reset_index(), TAXONOMY.reset_index())
    assert list(table.columns[:3]) == ["ID", "otu0000", "otu0001"] and table.iloc[-1, 0] == "taxonomy"
    counts = count_matrix(table)
    # Samples as rows, the same whole numbers (no percentages): every sample keeps its sequencing depth.
    np.testing.assert_array_equal(counts.to_numpy(), OTU.to_numpy().T)
    assert counts.sum(axis=1).tolist() == OTU.sum(axis=0).tolist()
    assert profile_counts(table).value_kind == "counts"
    # Lineages standardised: exactly the example's own (already standard) names.
    assert list(counts.columns) == list(count_matrix(read_table(EXAMPLE_COUNTS)).columns)
    assert details == {"otus": 263, "samples": 40, "unassigned": 0, "lineages": counts.columns.nunique()}


def test_otus_with_the_same_lineage_stay_separate():
    otu = pd.DataFrame({"s1": [3, 4, 0], "s2": [1, 0, 7]}, index=pd.Index(["a", "b", "c"], name="OTU"))
    taxonomy = pd.DataFrame({"Taxon": ["d__Bacteria; p__Firmicutes; g__Blautia(0.9)",
                                       "D_0__Bacteria;D_1__Firmicutes;D_5__Blautia", "Unassigned"]},
                            index=otu.index)
    table, details = yamas.counts_from_16s(otu.reset_index(), taxonomy.reset_index())
    counts = count_matrix(table)
    assert counts.shape == (2, 3)  # two OTUs of Blautia are two columns, not summed
    assert list(counts.columns) == ["k__Bacteria;p__Firmicutes;c__;o__;f__;g__Blautia"] * 2 + ["Unassigned"]
    assert details == {"otus": 3, "samples": 2, "unassigned": 1, "lineages": 2}


def test_mipmlp_gets_the_same_table_as_from_the_example_itself(tmp_path):
    """Faithfulness: the example as YaMAS files (messy names) preprocesses exactly like the example."""
    table, _ = yamas.counts_from_16s(OTU.reset_index(), TAXONOMY.reset_index())
    ours = preprocess(table).table
    reference = preprocess(read_table(EXAMPLE_COUNTS)).table
    pd.testing.assert_frame_equal(ours, reference, check_names=False)


# --- recognised among the uploads ------------------------------------------------------------------

@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_the_pair_is_recognised_in_any_order_and_under_any_names(order):
    files = _files(OTU, TAXONOMY, names=("run7_table.csv", "lineages.csv"))
    inputs = sort_uploads([files[i] for i in order])
    assert inputs.counts_layout == "yamas-16s"
    assert inputs.counts_name == "run7_table.csv + lineages.csv"
    assert inputs.counts_details["files"] == ["run7_table.csv", "lineages.csv"]
    assert inputs.labels is None and inputs.metabolites is None


def test_the_pair_with_labels_and_metabolites_is_four_files():
    ids = OTU.columns
    metabolites = pd.DataFrame(np.random.default_rng(0).uniform(1, 100, (len(ids), 4)),
                               columns=[f"met{i}" for i in range(4)], index=pd.Index(ids, name="ID"))
    files = [TAGS, *_files(OTU, TAXONOMY), ("metabolome.csv", metabolites.to_csv().encode())]
    assert len(files) == MAX_FILES == 4
    inputs = sort_uploads(files)
    assert inputs.counts_layout == "yamas-16s"
    assert inputs.labels_name == EXAMPLE_TAGS.name and inputs.metabolites_name == "metabolome.csv"


def test_five_files_are_too_many():
    with pytest.raises(InputError, match="at most 4"):
        sort_uploads([TAGS, *_files(OTU, TAXONOMY), TAGS, TAGS])


@pytest.mark.parametrize("problem, message", [
    ("taxonomy alone", "no uploaded table lists its OTUs"),
    ("missing OTU", "not in 'taxonomy.csv'"),
    ("percentages", "does not hold read counts"),
    ("text in counts", "not numbers"),
    ("two taxonomies", "both look like YaMAS taxonomy files"),
])
def test_problems_are_explained(problem, message):
    otu, taxonomy = OTU.copy(), TAXONOMY.copy()
    if problem == "taxonomy alone":
        files = [_files(otu, taxonomy)[1], TAGS]
    elif problem == "missing OTU":
        files = _files(otu, taxonomy.iloc[1:])
    elif problem == "percentages":
        files = _files(otu.div(otu.sum(axis=0), axis=1) * 100, taxonomy)
    elif problem == "text in counts":
        otu = otu.astype(object)
        otu.iloc[2, 3] = "low"  # (an empty cell, or "n/a", is read as 0, as in any counts table)
        files = _files(otu, taxonomy)
    else:
        files = [*_files(otu, taxonomy), ("taxonomy2.csv", taxonomy.to_csv().encode())]
    with pytest.raises(UserFacingError, match=message) as caught:
        sort_uploads(files)
    assert caught.value.fix  # every message says how to fix it


def test_an_ordinary_counts_table_is_not_affected():
    inputs = sort_uploads([(EXAMPLE_COUNTS.name, EXAMPLE_COUNTS.read_bytes()), TAGS])
    assert inputs.counts_layout == "row" and inputs.counts_details is None


# --- through the website ---------------------------------------------------------------------------

def test_uploaded_through_the_site_the_data_check_says_how_it_was_read(tmp_path):
    client = create_app({"TESTING": True, "JOBS": "inline", "RUNS_DIR": tmp_path / "runs",
                         "CACHE_DIR": tmp_path / "cache"}).test_client()
    data = {"files": [(io.BytesIO(content), name) for name, content in [*_files(OTU, TAXONOMY), TAGS]]}
    response = client.post("/upload", data=data, content_type="multipart/form-data", follow_redirects=True)
    page = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Read as YaMAS 16S output: otu.csv (263 OTUs × 40 samples)" in page
    assert "the read counts kept as they are" in page and "Median reads" in page  # counts, with depth


def test_a_whole_analysis_runs_on_yamas_files(tmp_path):
    run = workflows.start_run(tmp_path / "runs", sort_uploads(_files(OTU, TAXONOMY)))
    workflows.run_goal(run, "preprocess", tmp_path / "cache")
    assert run.meta()["status"] == "done" and run.meta()["counts_layout"] == "yamas-16s"
    names = {p.stem for p in run.figures_dir.glob("*.png")}
    assert {"reads_per_sample", "alpha_diversity", "composition_bars"} <= names  # counts: depth is shown
