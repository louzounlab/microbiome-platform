"""MetaPhlAn 4 output (shotgun) as an upload (core/metaphlan.py): species kept, % or read counts."""

import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microplatform import workflows
from microplatform.core import metaphlan
from microplatform.core.inputs import InputError, sort_uploads
from microplatform.core.profiler import count_matrix, profile_counts
from microplatform.examples import EXAMPLE_COUNTS
from microplatform.web import create_app

REAL = Path(__file__).parents[1] / "docs" / "yamas_examples" / "ERP020710_Pregnant_SG.csv"
DB = "#mpa_vJun23_CHOCOPhlAnSGB_202307"
SPECIES = {  # lineage (MetaPhlAn's '|' form) -> per-sample abundance in %, three samples
    "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Blautia|s__Blautia_obeum":
        [60.0, 10.5, 0.0],
    "k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_uniformis":
        [39.95, 89.5, 99.0],
    "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii":
        [0.05, 0.0, 1.0],
}
SAMPLES = ["S1", "S2", "S3"]


def _all_ranks(species: dict[str, list[float]]) -> dict[str, list[float]]:
    """Every rank, as MetaPhlAn writes it: each parent clade is the sum of its children, plus one SGB per species."""
    out: dict[str, np.ndarray] = {}
    for lineage, values in species.items():
        ranks = lineage.split("|")
        for depth in range(1, len(ranks) + 1):
            clade = "|".join(ranks[:depth])
            out[clade] = out.get(clade, 0) + np.array(values, dtype=float)
        out[lineage + "|t__SGB" + str(len(lineage))] = np.array(values, dtype=float)
    return {k: list(v) for k, v in out.items()}


def merged_table(species=SPECIES, samples=SAMPLES) -> bytes:
    """merge_metaphlan_tables.py output: clades as rows, one tab-separated column per sample."""
    rows = [DB, "clade_name\tclade_taxid\t" + "\t".join(samples)]
    rows += [f"{c}\t123|456\t" + "\t".join(f"{v:g}" for v in vals) for c, vals in _all_ranks(species).items()]
    return ("\n".join(rows) + "\n").encode()


def samples_as_rows(species=SPECIES, samples=SAMPLES) -> bytes:
    """The layout of the lab's ERP020710 file: samples as rows, comma-separated, ranks separated by ','."""
    clades = _all_ranks(species)
    table = pd.DataFrame({c.replace("|", ","): v for c, v in clades.items()}, index=pd.Index(samples, name="clade_name"))
    return (DB + "\n" + table.to_csv()).encode()


def read_stats_profile(reads: dict[str, int], sample: str = "patient7") -> bytes:
    """One sample's profile from ``metaphlan -t rel_ab_w_read_stats``."""
    total = sum(reads.values())
    rows = [DB, "#/usr/bin/metaphlan sample.fastq -t rel_ab_w_read_stats", f"#SampleID\t{sample}",
            "#estimated_reads_mapped_to_known_clades:" + str(total),
            "#clade_name\tclade_taxid\trelative_abundance\tcoverage\testimated_number_of_reads_from_the_clade"]
    for clade, n in _all_ranks({k: [v] for k, v in reads.items()}).items():
        rows.append(f"{clade}\t1|2\t{100 * n[0] / total:.5f}\t0.5\t{int(n[0])}")
    return ("\n".join(rows) + "\n").encode()


def _counts(name: str, content: bytes):
    inputs = sort_uploads([(name, content)])
    return inputs, count_matrix(inputs.counts)


@pytest.mark.parametrize("make", [merged_table, samples_as_rows], ids=["clades-as-rows", "samples-as-rows"])
def test_only_the_species_are_kept_in_either_orientation(make):
    inputs, counts = _counts("profile.csv", make())
    assert inputs.counts_layout == "metaphlan"
    assert list(counts.index) == SAMPLES
    assert counts.shape == (3, 3)  # three species, not the 8 ranks + SGB rows
    assert all(c.startswith("k__") and ";s__" in c for c in counts.columns)  # standard form, '|' and ',' alike
    np.testing.assert_allclose(counts.sum(axis=1), 100)  # nothing counted twice
    assert inputs.counts_details["values"] == "relative"
    assert inputs.counts_details["species"] == 3


def test_relative_abundance_is_passed_on_as_it_is_and_the_limitation_is_measured():
    inputs, counts = _counts("profile.tsv", merged_table())
    blautia = next(c for c in counts.columns if "Blautia" in c)
    assert counts.loc["S2", blautia] == 10.5
    assert profile_counts(inputs.counts).value_kind == "relative"
    present = np.array([v for vals in SPECIES.values() for v in vals if v > 0])
    assert inputs.counts_details["below_low_percent"] == round((present < 0.1).mean() * 100, 1)


def test_read_counts_are_kept_as_counts():
    reads = {"S1": [6000, 105, 0], "S2": [3995, 895, 990], "S3": [5, 0, 10]}
    species = {lin: [reads[s][i] for s in SAMPLES] for i, lin in enumerate(SPECIES)}
    inputs, counts = _counts("reads.tsv", merged_table(species))
    assert inputs.counts_details["values"] == "counts"
    assert (counts.dtypes == "int64").all()
    assert counts.loc["S1"].tolist() == [6000, 105, 0]
    profile = profile_counts(inputs.counts)
    assert profile.value_kind == "counts" and profile.median_reads > 0  # the depth checks apply


def test_a_read_stats_profile_uses_its_estimated_reads():
    reads = dict(zip(SPECIES, [1200, 800, 0]))
    inputs, counts = _counts("sample.txt", read_stats_profile(reads))
    assert inputs.counts_details["values"] == "counts"
    assert list(counts.index) == ["patient7"]  # from MetaPhlAn's #SampleID line
    assert sorted(counts.loc["patient7"].tolist()) == [0, 800, 1200]


def test_unclassified_is_left_out_and_its_share_reported():
    content = merged_table().decode() + "UNCLASSIFIED\t-1\t20\t0\t50\n"
    inputs, counts = _counts("profile.tsv", content.encode())
    assert not any("UNCLASSIFIED" in c for c in counts.columns)
    assert inputs.counts_details["unclassified_percent"] > 0


@pytest.mark.parametrize("content, message", [
    (f"{DB}\nclade_name\tS1\nk__Bacteria\t100\nk__Bacteria|p__Firmicutes\t100\n", "no species-level"),
    (f"{DB}\nclade_name\tS1\tS2\nk__Bacteria|p__F|c__C|o__O|f__F|g__G|s__G_x\tabc\t1\n", "not numbers"),
    (f"{DB}\nclade_name\tS1\tS1\nk__Bacteria|p__F|c__C|o__O|f__F|g__G|s__G_x\t1\t1\n", "more than once"),
])
def test_problems_are_explained(content, message):
    with pytest.raises(InputError, match=message):
        sort_uploads([("profile.tsv", content.encode())])


def test_an_ordinary_counts_table_is_not_taken_for_metaphlan():
    assert not metaphlan.is_metaphlan(EXAMPLE_COUNTS.read_bytes())
    assert sort_uploads([(EXAMPLE_COUNTS.name, EXAMPLE_COUNTS.read_bytes())]).counts_layout == "row"


def test_the_data_check_says_how_it_was_read_and_states_the_limitation(tmp_path):
    client = create_app({"TESTING": True, "JOBS": "inline", "RUNS_DIR": tmp_path / "runs",
                         "CACHE_DIR": tmp_path / "cache"}).test_client()
    data = {"files": [(io.BytesIO(samples_as_rows()), "profile.csv")]}
    response = client.post("/upload", data=data, content_type="multipart/form-data", follow_redirects=True)
    page = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Read as MetaPhlAn output (mpa_vJun23_CHOCOPhlAnSGB_202307)" in page and "3 species are kept" in page
    assert "relative abundances (%), not read counts" in page and "rel_ab_w_read_stats" in page


def test_a_whole_analysis_runs_on_metaphlan_output(tmp_path):
    run = workflows.start_run(tmp_path / "runs", sort_uploads([("profile.tsv", merged_table())]))
    workflows.run_goal(run, "preprocess", tmp_path / "cache")
    assert run.meta()["status"] == "done" and run.meta()["counts_layout"] == "metaphlan"
    assert (run.path / "processed.csv").exists()


@pytest.mark.skipif(not REAL.exists(), reason="the lab's real MetaPhlAn file is not in this checkout")
def test_the_labs_real_file():
    inputs, counts = _counts(REAL.name, REAL.read_bytes())
    assert counts.shape == (54, 1268)
    assert inputs.counts_details["values"] == "relative"
    np.testing.assert_allclose(counts.sum(axis=1), 100, atol=0.01)
