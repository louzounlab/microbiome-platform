import io
import json
import re
import zipfile

import pandas as pd
import pytest

from microplatform import workflows
from microplatform.core.preprocessing import CANONICAL
from microplatform.core.profiler import count_matrix, read_table
from microplatform.examples import EXAMPLE_COUNTS
from microplatform.reports.figures import (
    FORMATS,
    PLOTLY_JS,
    export,
    load_manifest,
    preprocessing_figures,
    reads_vs_taxa,
    taxa_prevalence,
    to_plotly,
)


@pytest.fixture(scope="module")
def figures():
    counts = count_matrix(read_table(EXAMPLE_COUNTS))
    return preprocessing_figures(counts, CANONICAL) + [taxa_prevalence(counts, CANONICAL)]


@pytest.fixture(scope="module")
def exported(figures, tmp_path_factory):
    out = tmp_path_factory.mktemp("figures")
    return out, export(figures, out)


def test_reads_per_sample(figures):
    reads = figures[0]
    assert reads.name == "reads_per_sample" and reads.log_y
    assert len(reads.data) == 40
    assert reads.data["reads"].is_monotonic_decreasing


def test_reads_vs_taxa_comes_right_after_the_depth_figure(figures):
    assert [f.name for f in figures[:2]] == ["reads_per_sample", "reads_vs_taxa"]
    scatter = figures[1]
    assert scatter.kind == "depth_taxa" and len(scatter.data) == 40
    assert set(scatter.data.columns) == {"sample", "reads", "taxa", "flag", "named_on_plot"}
    assert "check before removing" in scatter.caption and "25% of the median, or under 1,000 reads" in scatter.caption


def _counts_with_weak_samples() -> pd.DataFrame:
    """30 ordinary samples (~10,000 reads, 40 taxa), plus a failed one, one shallow and one with few taxa."""
    rows = {f"s{i:02d}": [250 + i] * 40 + [0] * 10 for i in range(30)}
    rows["failed"] = [20] * 5 + [0] * 45  # 100 reads, 5 taxa: low on both
    rows["shallow"] = [30] * 40 + [0] * 10  # 1,200 reads (below a quarter of ~10,000), 40 taxa
    rows["few_taxa"] = [1250] * 8 + [0] * 42  # 10,000 reads, 8 taxa (below a quarter of 40)
    rows["empty"] = [0] * 50  # no reads: left out (the data check warns separately)
    return pd.DataFrame.from_dict(rows, orient="index", columns=[f"t{j}" for j in range(50)])


def test_reads_vs_taxa_flags_and_names_the_weak_samples():
    fig = reads_vs_taxa(_counts_with_weak_samples())
    flags = fig.data.set_index("sample")["flag"]
    assert flags["failed"] == "low on both" and flags["shallow"] == "low reads" and flags["few_taxa"] == "low taxa"
    assert (flags.drop(["failed", "shallow", "few_taxa"]) == "within range").all()
    assert "empty" not in flags.index and "1 sample with no reads is not shown" in fig.caption
    named = set(fig.data.loc[fig.data["named_on_plot"], "sample"])
    assert named == {"failed", "shallow", "few_taxa"}  # only flagged samples are named
    assert fig.marker[0] > 1000 and fig.hline[0] == pytest.approx(10)  # a quarter of the medians


def test_reads_vs_taxa_names_at_most_ten_samples_the_most_extreme_first():
    counts = _counts_with_weak_samples()
    for i in range(14):
        counts.loc[f"weak{i:02d}"] = [5 + i] * 4 + [0] * 46
    fig = reads_vs_taxa(counts)
    named = fig.data.loc[fig.data["named_on_plot"], "sample"].tolist()
    assert named == [f"weak{i:02d}" for i in range(10)]  # the 10 furthest below the cut-offs


def test_reads_vs_taxa_cut_offs_can_be_tuned():
    counts = _counts_with_weak_samples()
    strict = reads_vs_taxa(counts, fraction=0.05, min_reads=50)
    assert strict.data.set_index("sample")["flag"]["shallow"] == "within range"
    assert reads_vs_taxa(counts, min_reads=20_000).marker[0] == 20_000  # the floor wins when higher


def test_reads_vs_taxa_interactive_has_bands_a_log_axis_and_the_legend_above():
    plot = to_plotly(reads_vs_taxa(_counts_with_weak_samples()))
    assert plot.layout.xaxis.type == "log" and plot.layout.legend.y > 1
    assert [t.name for t in plot.data] == ["Within range", "Low on reads or on taxa (2)",
                                          "Low on both: likely failed (1)"]
    rects = [s for s in plot.layout.shapes if s.type == "rect"]
    assert len(rects) == 3  # the two low bands and the corner where both are low (drawn apart, not overlaid)


def test_relative_abundances_get_no_reads_vs_taxa():
    counts = _counts_with_weak_samples().drop(index="empty")
    relative = counts.div(counts.sum(axis=1), axis=0)
    assert preprocessing_figures(relative, CANONICAL) == []


def test_prevalence_marks_the_mipmlp_cut_off(figures):
    prevalence = figures[-1]
    assert len(prevalence.data) == 263
    assert prevalence.data["prevalence_percent"].between(0, 100).all()
    assert prevalence.marker[0] == CANONICAL.rare_bacteria_threshold * 100


def test_every_format_is_written(exported):
    out, figs = exported
    for fig in figs:
        for fmt in (*FORMATS, "json"):
            assert (out / fig.files[fmt]).stat().st_size > 0


def test_data_csv_is_exactly_what_is_plotted(figures, exported):
    out, _ = exported
    for fig in figures:
        pd.testing.assert_frame_equal(pd.read_csv(out / f"{fig.name}.csv"), fig.data, check_dtype=False)


def test_svg_text_stays_editable(exported):
    out, _ = exported
    assert "<text" in (out / "reads_per_sample.svg").read_text(encoding="utf-8")


def test_interactive_html_works_offline(exported):
    out, _ = exported
    html = (out / "taxa_prevalence.html").read_text(encoding="utf-8")
    assert not re.search(r"<script[^>]*\bsrc=", html)  # nothing loaded from outside
    assert "Plotly.newPlot" in html and len(html) > 1_000_000  # plotly.js is embedded


def test_manifest_round_trip(exported):
    out, figs = exported
    assert load_manifest(out) == figs
    assert load_manifest(out / "missing") == []


def test_plotly_figure_has_log_axis_and_marker(figures):
    assert to_plotly(figures[0]).layout.yaxis.type == "log"
    assert to_plotly(figures[1]).layout.shapes  # the cut-off line


def test_plotly_js_is_bundled():
    assert PLOTLY_JS.is_file()


# --- through the workflow ------------------------------------------------------

@pytest.fixture(scope="module")
def finished_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("runs")
    run = workflows.start_example_run(root)
    workflows.run_goal(run, "preprocess", cache_dir=root / "cache")
    return run


def test_preprocessing_produces_figures(finished_run):
    view = workflows.describe(finished_run)
    assert [f.name for f in view.quality_figures] == ["reads_per_sample", "reads_vs_taxa"]
    assert [f.name for f in view.overview_figures] == ["composition_bars", "composition_pie", "alpha_diversity",
                                                       "beta_pcoa", "taxa_per_sample", "taxa_prevalence"]


def test_figure_files_cannot_escape_the_run(finished_run):
    assert finished_run.figure_file("reads_per_sample.png") is not None
    # "figures.json.." matters on Windows, which silently drops trailing dots.
    for name in ("../meta.json", "..\\meta.json", ".hidden", "nope.png", "figures.json..",
                 "READS_PER_SAMPLE.PNG"):
        assert finished_run.figure_file(name) is None


def test_bundle_has_everything(finished_run):
    with zipfile.ZipFile(io.BytesIO(workflows.build_bundle(finished_run))) as zf:
        names = set(zf.namelist())
        assert {"README.txt", "record.json", "results/processed_table.csv",
                "inputs/counts.csv", "inputs/labels.csv"} <= names
        assert {f"figures/2_data_quality/reads_per_sample.{fmt}" for fmt in FORMATS} <= names
        for fig in ("taxa_prevalence", "composition_bars", "beta_pcoa"):
            assert {f"figures/3_descriptive_overview/{fig}.{fmt}" for fmt in FORMATS} <= names
        assert not any(n.endswith(".json") and n.startswith("figures/") for n in names)
        assert json.loads(zf.read("record.json"))["tool"] == "MIPMLP"
        assert finished_run.id in zf.read("README.txt").decode()
