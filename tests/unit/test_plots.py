"""Unit tests for funlr.plots (matplotlib port of the stage-7 R plotting block).

Synthetic fixtures mimic stage2/all_domain_hits.tsv and
stage3/tables/tiered_candidates.tsv. The summary TSV assertions are
hand-computed from the fixture below.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest

from funlr import plots
from funlr.plots import classify_domain_type, generate_domain_plots

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


@pytest.fixture
def log() -> logging.Logger:
    logger = logging.getLogger("test_plots")
    logger.addHandler(logging.NullHandler())
    return logger


def _write_fixture(
    tmp_path: Path,
    hits_rows: list[tuple],
    tiered_rows: list[tuple],
    tiered_columns: list[str] | None = None,
) -> tuple[Path, Path]:
    hits = pd.DataFrame(hits_rows, columns=["protein_id", "domain", "start", "end"])
    hits_path = tmp_path / "all_domain_hits.tsv"
    hits.to_csv(hits_path, sep="\t", index=False)
    tiered = pd.DataFrame(
        tiered_rows,
        columns=tiered_columns or ["protein_id", "tier", "rescue_priority"],
    )
    tiered["protein_length"] = 1200
    tiered_path = tmp_path / "tiered_candidates.tsv"
    tiered.to_csv(tiered_path, sep="\t", index=False)
    return hits_path, tiered_path


# Hand-computed fixture:
#   p1/p3 (TIER_1A): NACHT + WD40            -> architecture "NACHT;WD40" x2
#   p2    (TIER_1A): NB-ARC + WD40           -> "NB-ARC;WD40" x1
#   p4    (TIER_2A): NACHT + LRR + HRAM      -> "HRAM;LRR;NACHT" x1
#   p5    (TIER_3B): HeLo + NACHT + TPR      -> "HeLo;NACHT;TPR" x1
#   p6    (TIER_4A): excluded from plots entirely (tier not in ^TIER_[123])
#   p7    (TIER_2D): candidate with no domain hits -> absent everywhere
#   p8    (FUSION_RESCUE, is_fusion_model): excluded despite domain hits
HITS_ROWS = [
    ("p1", "NACHT", 100, 300),
    ("p1", "WD40", 500, 700),
    ("p2", "NB-ARC", 100, 300),
    ("p2", "WD40", 500, 700),
    ("p3", "NACHT", 100, 300),
    ("p3", "WD40", 500, 700),
    ("p4", "NACHT", 100, 300),
    ("p4", "LRR", 500, 900),
    ("p4", "HRAM", 950, 1000),
    ("p5", "HeLo", 10, 80),
    ("p5", "NACHT", 100, 300),
    ("p5", "TPR", 400, 600),
    ("p6", "NACHT", 100, 300),
    ("p8", "NACHT", 100, 300),
]
TIERED_ROWS = [
    ("p1", "TIER_1A_HIGH_CONFIDENCE", 10, "False"),
    ("p2", "TIER_1A_HIGH_CONFIDENCE", 20, "False"),
    ("p3", "TIER_1A_HIGH_CONFIDENCE", 30, "False"),
    ("p4", "TIER_2A_HIGH_PRIORITY_RESCUE", 40, "False"),
    ("p5", "TIER_3B_ARCHITECTURAL_VARIANT", 50, "False"),
    ("p6", "TIER_4A_REPEAT_ONLY_NO_NBD", 90, "False"),
    ("p7", "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE", 60, "False"),
    ("p8", "FUSION_RESCUE", 100, "True"),
]
TIERED_COLUMNS = ["protein_id", "tier", "rescue_priority", "is_fusion_model"]

EXPECTED_SUMMARY = (
    "tier\ttotal_proteins\tarchitecture\tcount\n"
    "TIER_1A_HIGH_CONFIDENCE\t3\tNACHT;WD40\t2\n"
    "TIER_1A_HIGH_CONFIDENCE\t3\tNB-ARC;WD40\t1\n"
    "TIER_2A_HIGH_PRIORITY_RESCUE\t1\tHRAM;LRR;NACHT\t1\n"
    "TIER_3B_ARCHITECTURAL_VARIANT\t1\tHeLo;NACHT;TPR\t1\n"
)

SUMMARY_PNGS = [
    "sensor_distribution.png",
    "nbd_distribution.png",
    "asm_presence.png",
    "length_vs_priority.png",
    "order_sanity.png",
    "nlr_domains_by_tier.png",
]


def _main_fixture(tmp_path: Path) -> tuple[Path, Path]:
    return _write_fixture(tmp_path, HITS_ROWS, TIERED_ROWS, TIERED_COLUMNS)


# ---------------------------------------------------------------------------
# Domain token cascade
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "token,expected",
    [
        ("NACHT", "NACHT"),
        ("NB-ARC", "NB-ARC"),
        ("AAA_16", "AAA"),
        ("AAA_22", "AAA"),
        ("HeLo", "HeLo"),
        ("HeLo-like", "HeLo-like"),  # order matters: must not classify as HeLo
        ("HELL", "HeLo-like"),
        ("HET", "HET"),
        ("HET-s", "HET-s"),
        ("PF00168", "C2"),
        ("C2 domain", "C2"),
        ("NLR07", "PP"),
        ("PP-motif", "PP"),
        ("NLR05", "Basidio_ASM"),
        ("NLR17", "Lineage_ASM"),
        ("PUASM", "PUASM"),
        ("WD40__class__Agaricomycetes", "WD40"),  # taxonomic suffix stripped
        ("Leucine-rich repeat", "LRR"),
        ("Zinc finger", "Zinc_finger"),
        ("Coiled_coil", "CC"),
        ("some_unrelated_thing", "Other"),
    ],
)
def test_classify_domain_type(token: str, expected: str) -> None:
    assert classify_domain_type(token) == expected


def test_domain_category_map() -> None:
    assert plots.domain_category("NACHT") == "NBD"
    assert plots.domain_category("WD40") == "Sensor"
    assert plots.domain_category("HeLo") == "Effector"
    assert plots.domain_category("HRAM") == "ASM"
    assert plots.domain_category("Other") == "Other"


# ---------------------------------------------------------------------------
# Full artifact set on the hand-computed fixture
# ---------------------------------------------------------------------------

def test_generate_all_artifacts(tmp_path: Path, log: logging.Logger) -> None:
    hits_path, tiered_path = _main_fixture(tmp_path)
    outdir = tmp_path / "domain_plots"
    produced = generate_domain_plots(hits_path, tiered_path, outdir, log)

    expected = [outdir / name for name in SUMMARY_PNGS if name != "nlr_domains_by_tier.png"]
    expected += [outdir / "individual" / f"{pid}.png" for pid in ("p1", "p2", "p3", "p4", "p5")]
    expected += [outdir / "nlr_domains_by_tier.png"]
    expected += [outdir / "domain_architecture_summary.tsv"]
    assert produced == expected

    for path in produced:
        assert path.exists(), path
        if path.suffix == ".png":
            assert path.read_bytes()[:8] == PNG_MAGIC, path
    # fusions and Tier-4 proteins never plotted
    assert not (outdir / "individual" / "p6.png").exists()
    assert not (outdir / "individual" / "p8.png").exists()
    # no temp files left behind (atomic temp-then-rename writes)
    assert list(outdir.rglob("*.tmp")) == []


def test_summary_tsv_content(tmp_path: Path, log: logging.Logger) -> None:
    hits_path, tiered_path = _main_fixture(tmp_path)
    outdir = tmp_path / "domain_plots"
    generate_domain_plots(hits_path, tiered_path, outdir, log)
    text = (outdir / "domain_architecture_summary.tsv").read_text()
    assert text == EXPECTED_SUMMARY


# ---------------------------------------------------------------------------
# Defensive behavior
# ---------------------------------------------------------------------------

def test_missing_input_files_raise(tmp_path: Path, log: logging.Logger) -> None:
    hits_path, tiered_path = _main_fixture(tmp_path)
    with pytest.raises(FileNotFoundError, match="all_domain_hits"):
        generate_domain_plots(tmp_path / "nope.tsv", tiered_path, tmp_path / "o", log)
    with pytest.raises(FileNotFoundError, match="tiered_candidates"):
        generate_domain_plots(hits_path, tmp_path / "nope.tsv", tmp_path / "o", log)


def test_missing_required_column_raises(tmp_path: Path, log: logging.Logger) -> None:
    hits_path, tiered_path = _main_fixture(tmp_path)
    # strip the 'domain' column from the hits table
    hits = pd.read_csv(hits_path, sep="\t")
    hits.drop(columns=["domain"]).to_csv(hits_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="domain"):
        generate_domain_plots(hits_path, tiered_path, tmp_path / "o", log)


def test_no_tier123_candidates_skips_gracefully(
    tmp_path: Path, log: logging.Logger
) -> None:
    hits_path, tiered_path = _write_fixture(
        tmp_path,
        HITS_ROWS,
        [("p6", "TIER_4A_REPEAT_ONLY_NO_NBD", 90),
         ("px", "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL", 95)],
    )
    outdir = tmp_path / "domain_plots"
    assert generate_domain_plots(hits_path, tiered_path, outdir, log) == []
    assert not list(outdir.rglob("*.png"))


def test_no_asm_domains_still_plots(tmp_path: Path, log: logging.Logger) -> None:
    """asm_presence.png is produced even when no ASM-category hit exists."""
    hits_path, tiered_path = _write_fixture(
        tmp_path,
        [("p1", "NACHT", 100, 300), ("p1", "WD40", 500, 700)],
        [("p1", "TIER_1A_HIGH_CONFIDENCE", 10)],
    )
    outdir = tmp_path / "domain_plots"
    produced = generate_domain_plots(hits_path, tiered_path, outdir, log)
    assert outdir / "asm_presence.png" in produced
    summary = (outdir / "domain_architecture_summary.tsv").read_text()
    assert "NACHT;WD40" in summary


def test_missing_asm_column_in_tiered_ok(tmp_path: Path, log: logging.Logger) -> None:
    """tiered_candidates.tsv without ASM/rescue annotations still plots."""
    hits_path, tiered_path = _write_fixture(
        tmp_path,
        [("p1", "NACHT", 100, 300), ("p1", "HRAM", 500, 700)],
        [("p1", "TIER_1A_HIGH_CONFIDENCE")],
        tiered_columns=["protein_id", "tier"],
    )
    outdir = tmp_path / "domain_plots"
    produced = generate_domain_plots(hits_path, tiered_path, outdir, log)
    # length_vs_priority needs rescue_priority -> skipped with a warning
    names = [p.name for p in produced]
    assert "length_vs_priority.png" not in names
    for name in ("sensor_distribution.png", "nbd_distribution.png",
                 "asm_presence.png", "order_sanity.png",
                 "nlr_domains_by_tier.png", "domain_architecture_summary.tsv"):
        assert name in names


def test_candidates_without_hits_graceful(tmp_path: Path, log: logging.Logger) -> None:
    """Tier 1-3 candidates exist but none have domain hits."""
    hits_path, tiered_path = _write_fixture(
        tmp_path,
        [("other", "NACHT", 100, 300)],
        [("p1", "TIER_1A_HIGH_CONFIDENCE", 10)],
    )
    outdir = tmp_path / "domain_plots"
    produced = generate_domain_plots(hits_path, tiered_path, outdir, log)
    summary = outdir / "domain_architecture_summary.tsv"
    assert summary in produced
    assert summary.read_text() == "tier\ttotal_proteins\tarchitecture\tcount\n"
    assert list((outdir / "individual").glob("*.png")) == []
    for name in SUMMARY_PNGS:
        assert (outdir / name).exists()


def test_actual_sequence_length_is_used_for_scatter(tmp_path, log, monkeypatch):
    """An unannotated C-terminus must not shorten the plotted protein."""
    hits, tiered = _write_fixture(
        tmp_path, [("p1", "NACHT", 100, 300)],
        [("p1", "TIER_1A_HIGH_CONFIDENCE", 10)],
    )
    observed = []
    original = plots._scatter_plot

    def capture(path, features, tier_levels):
        observed.extend(features["protein_length"].tolist())
        return original(path, features, tier_levels)

    monkeypatch.setattr(plots, "_scatter_plot", capture)
    generate_domain_plots(hits, tiered, tmp_path / "plots", log)
    assert observed == [1200.0]


def test_scatter_retains_visible_tier_and_sensor_legends(tmp_path, monkeypatch):
    from matplotlib.legend import Legend
    features = pd.DataFrame([
        {"protein_id": "p1", "tier": "TIER_1A_HIGH_CONFIDENCE",
         "sensor_class": "WD40", "protein_length": 1200, "rescue_priority": 10},
        {"protein_id": "p2", "tier": "TIER_2A_HIGH_PRIORITY_RESCUE",
         "sensor_class": "None", "protein_length": 600, "rescue_priority": 100},
    ])
    seen = {}

    def inspect(fig, path):
        fig.canvas.draw()
        legends = [child for child in fig.axes[0].get_children()
                   if isinstance(child, Legend)]
        seen["titles"] = {legend.get_title().get_text() for legend in legends}
        renderer = fig.canvas.get_renderer()
        for legend in legends:
            box = legend.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0
            assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
        plots.plt.close(fig)
        return path

    monkeypatch.setattr(plots, "_atomic_save_fig", inspect)
    plots._scatter_plot(tmp_path / "scatter.png", features,
                        ["TIER_1A_HIGH_CONFIDENCE", "TIER_2A_HIGH_PRIORITY_RESCUE"])
    assert seen["titles"] == {"Tier", "Sensor class"}


@pytest.mark.parametrize("protein_id", ["../escape", "p/child", "p\\child", "..", "p\nchild"])
def test_plot_filenames_cannot_escape_output(tmp_path, log, protein_id):
    hits, tiered = _write_fixture(tmp_path, [(protein_id, "NACHT", 1, 100)],
                                 [(protein_id, "TIER_1A_HIGH_CONFIDENCE", 10)])
    with pytest.raises(ValueError, match="plot filename"):
        generate_domain_plots(hits, tiered, tmp_path / "plots", log)
    assert not (tmp_path / "plots").exists()


@pytest.mark.parametrize("start,end", [(0, 100), (100, 1), (1, "inf"), ("bad", 100)])
def test_invalid_domain_coordinates_are_rejected(tmp_path, log, start, end):
    hits, tiered = _write_fixture(tmp_path, [("p1", "NACHT", start, end)],
                                 [("p1", "TIER_1A_HIGH_CONFIDENCE", 10)])
    with pytest.raises(ValueError, match="invalid start/end"):
        generate_domain_plots(hits, tiered, tmp_path / "plots", log)


def test_no_length_proxy_when_sequence_length_missing(tmp_path, log):
    hits, tiered = _main_fixture(tmp_path)
    frame = pd.read_csv(tiered, sep="\t")
    frame.drop(columns=["protein_length"]).to_csv(tiered, sep="\t", index=False)
    with pytest.raises(ValueError, match="protein_length"):
        generate_domain_plots(hits, tiered, tmp_path / "plots", log)


def test_stage7_default_renderer_writes_validated_outputs_without_r(tmp_path, log, monkeypatch):
    import json
    from types import SimpleNamespace
    from funlr.config import FunlrConfig
    from funlr.core.context import RunPaths
    from funlr.stages.reporting_v3 import domain_plots

    paths = RunPaths(tmp_path / "run")
    stage7 = paths.stage_dir(7)
    for name in ("logs", "metadata"):
        (stage7 / name).mkdir()
    hits = pd.DataFrame([("p1", "NACHT", 1, 100), ("fusion", "NACHT", 1, 100)],
                        columns=["protein_id", "domain", "start", "end"])
    hits.to_csv(paths.stage_dir(2) / "all_domain_hits.tsv", sep="\t", index=False)
    tiered = pd.DataFrame([("p1", "TIER_1A_HIGH_CONFIDENCE", 1200, 10, False),
                           ("fusion", "TIER_1A_HIGH_CONFIDENCE", 1200, 10, True)],
                         columns=["protein_id", "tier", "protein_length", "rescue_priority", "is_fusion_model"])
    monkeypatch.setenv("MPLCONFIGDIR", str(tmp_path / "inherited-cache"))
    # Deliberately omit a command runner: the default renderer must not call R.
    ctx = SimpleNamespace(config=FunlrConfig(), paths=paths, logger=log)
    products = domain_plots(ctx, tiered, stage7, paths.final_dir)
    status = json.loads((stage7 / "metadata/plot_status.json").read_text())
    assert status["backend"] == "matplotlib"
    assert status["protein_length_source"] == "Stage3 protein_length"
    for directory in (stage7 / "domain_plots", paths.final_dir / "domain_plots"):
        assert (directory / "length_vs_priority.png") in products
        assert (directory / "length_vs_priority.png").read_bytes()[:8] == PNG_MAGIC
        assert (directory / "individual/p1.png").is_file()
        assert not (directory / "individual/fusion.png").exists()
        assert not (directory / "R_sessionInfo.txt").exists()
        provenance = json.loads((directory / "plot_provenance.json").read_text())
        assert provenance["renderer"].lower() == "agg"
        assert provenance["matplotlib_version"] == plots.matplotlib.__version__
        assert provenance["configured_cache"] == str((stage7 / "plot_work/.matplotlib").resolve())
        assert set(provenance["font_files"]) == {"normal", "bold"}
        assert all(len(row["sha256"]) == 64 for row in provenance["font_files"].values())
        assert not (directory / ".matplotlib").exists()
