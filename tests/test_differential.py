"""Fresh-run reproducibility for the updated reference scientific contract.

The earlier prototype uses different discovery/tiering rules and is no longer
an equality oracle. This check runs the current 12-protein fixture twice through
independent CLI processes, output paths, hash seeds, and equivalent nested/flat
configurations. It verifies all selected scientific outputs exactly and proves
that the comparison detects changed classifications. Native reference validation is
recorded separately; the executables here are explicit test doubles.
"""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


smoke = load_module("funlr_synthetic_fixture_helpers", REPO / "tests/integration/test_smoke.py")
comparison = load_module("funlr_comparison_helper", REPO / "scripts/compare_runs.py")


def _run(fixture, out, env, config):
    args = [sys.executable, "-m", "funlr", "run", "--genome", str(fixture / "genome.fa"),
            "--proteins", str(fixture / "proteins.faa"), "--gff3", str(fixture / "annotation.gff3"),
            "--eggnog", str(fixture / "emapper.annotations"), "--pfam", str(fixture / "db/pfam.hmm"),
            "--nbd-hmms", str(fixture / "db/nbd.hmm"), "--custom-hmms", str(fixture / "db/custom.hmm"),
            "--outdir", str(out), "--sample", "differential", "--threads", "2", "--config", str(config)]
    result = subprocess.run(args, cwd=REPO, env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_updated_science_reproduces_across_fresh_runs_and_equivalent_configs(tmp_path):
    fixture = smoke._prepare_fixture(tmp_path / "fixture")
    shim_dir = tmp_path / "shims"
    shutil.copytree(REPO / "tests/integration/shims", shim_dir)
    for path in shim_dir.iterdir():
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
    env = dict(os.environ, FUNLR_TEST_DOUBLES="1", FUNLR_FIXTURE_DIR=str(fixture),
               PYTHONPATH=str(REPO / "src"))
    env["PATH"] = os.pathsep.join((str(shim_dir), str(Path(sys.executable).parent), env.get("PATH", "")))
    nested = tmp_path / "nested.json"
    nested.write_text(json.dumps({"discovery": {"discovery_mode": "BROAD"}, "asm": {"asm_enable": 0},
                                  "fusion": {"enabled": 0}, "reporting": {"plots": 0}}))
    flat = tmp_path / "flat.json"
    flat.write_text(json.dumps({"DISCOVERY_MODE": "BROAD", "ASM_ENABLE": 0, "ENABLE_FUSION": 0, "REPORT_PLOTS": 0}))
    before, after = tmp_path / "first run", tmp_path / "second run"
    _run(fixture, before, dict(env, PYTHONHASHSEED="11"), nested)
    _run(fixture, after, dict(env, PYTHONHASHSEED="97"), flat)
    for root in (before, after):
        with (root / "final_results/nlr_final_report.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert {row["protein_id"]: row["tier"] for row in rows} == smoke.EXPECTED_TIERS
        assert len(rows) == 12
        with (root / "results/stage3/tables/rescue_priority.tsv").open() as handle:
            order = [(int(row["rescue_priority"]), row["protein_id"]) for row in csv.DictReader(handle, delimiter="\t")]
        assert order == sorted(order, key=lambda item: (-item[0], item[1]))
    report = comparison.compare_runs(before, after)
    assert report["status"] == "match", json.dumps(report["differences"], indent=2)
    assert report["checked_files"] >= 100  # includes both rescue tracks and dynamic tier exports
    changed = after / "final_results/nlr_final_report.tsv"
    changed.write_text(changed.read_text().replace("TIER_1A_HIGH_CONFIDENCE", "TIER_4G_SUSPECTED_PSEUDOGENE", 1))
    report = comparison.compare_runs(before, after)
    assert report["status"] == "different"
    assert [(item["path"], item["status"]) for item in report["differences"]] == [("final_results/nlr_final_report.tsv", "different")]
