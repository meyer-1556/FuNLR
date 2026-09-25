"""End-to-end smoke test: full `funlr run` on the synthetic fixture using
shim executables for all external bioinformatics tools.

The fixture (tests/fixtures/synthetic) contains 4 scaffolds and 12 proteins
designed to exercise canonical, exclusion, and rescue tiers, several flags, and the miniprot->exonerate
rescue handoff. The shims (tests/integration/shims) emit canned
domtblout/GFF/FAIDX outputs read from the fixture's canned/ directory.

Run from the repo root:
    PYTHONPATH=src python3 -m pytest tests/integration/test_smoke.py
"""
import os
import json
import hashlib
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "synthetic"
SHIM_SRC = REPO / "tests" / "integration" / "shims"

# Updated reference tiers under explicit BROAD discovery, with ASM/fusions/plots off.
# HMM evidence gives high confidence to the seven strong NACHT hits and medium
# confidence to ORC4's weaker hit. No biology is inferred from the random FASTAs.
EXPECTED_TIERS = {
    "FUN_000001-T1": "TIER_1A_HIGH_CONFIDENCE",          # NACHT+WD40, interior
    "FUN_000002-T1": "TIER_2A_HIGH_PRIORITY_RESCUE",     # NBD-only + CONTIG_END
    "FUN_000003-T1": "TIER_4B_LIKELY_STAND_HOUSEKEEPING",  # ORC4-like
    "FUN_000004-T1": "TIER_4A_REPEAT_ONLY_NO_NBD",       # Kelch-only
    "FUN_000005-T1": "TIER_1B_NEEDS_REVIEW",             # NLR arch at contig end
    "FUN_000006-T1": "TIER_3B_ARCHITECTURAL_VARIANT",     # NBD+AAA_16, no repeat
    "FUN_000007-T1": "TIER_4A_REPEAT_ONLY_NO_NBD",       # short WD40-only
    "FUN_000008-T1": "TIER_4A_REPEAT_ONLY_NO_NBD",       # updated TPR_1/2 grouping
    "FUN_000009-T1": "TIER_2B_RESCUE_CANDIDATE",         # strict NBD-only, isolated interior
    "FUN_000010-T1": "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL",  # HET-keyword quirk
    "FUN_000011-T1": "TIER_4E_ORIENTATION_INVALID",      # sensor before NBD is excluded
    "FUN_000012-T1": "TIER_1A_HIGH_CONFIDENCE",          # FP away from NBD flags/penalizes only
}

EXPECTED_STRICT = {
    "FUN_000001-T1", "FUN_000002-T1", "FUN_000005-T1", "FUN_000006-T1",
    "FUN_000009-T1", "FUN_000012-T1",
}


def _prepare_fixture(destination):
    """Upgrade old canned parser placeholders in an isolated fixture copy."""
    shutil.copytree(FIXTURE, destination)
    # Updated architecture consumes HMM length/bitscore; populate valid numeric
    # domtbl columns instead of the old parser's unused '-' placeholders.
    lengths = {}
    current = None
    for line in (destination / "proteins.faa").read_text().splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            lengths[current] = 0
        else:
            lengths[current] += len(line)
    for path in (destination / "canned").glob("*.domtblout"):
        rows = []
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                rows.append(line)
                continue
            fields = line.split()
            fields[2] = "400"
            fields[5] = str(lengths[fields[3]])
            fields[6], fields[7], fields[8] = fields[12], "100", "0"
            fields[9], fields[10], fields[11] = "1", "1", fields[12]
            fields[13], fields[14], fields[15] = "100", "0", "1"
            fields[16] = str(int(fields[18]) - int(fields[17]) + 1)
            fields[19], fields[20], fields[21] = fields[17], fields[18], "0.95"
            rows.append(" ".join(fields))
        path.write_text("\n".join(rows)+("\n" if rows else ""))
    return destination


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    """Run the full pipeline once per module; share results across tests."""
    run = tmp_path_factory.mktemp("funlr_smoke")
    fixture = _prepare_fixture(run / "fixture")

    # Stage shims into a writable dir (the repo mount may drop exec bits).
    shim_dir = run / "shims"
    shutil.copytree(SHIM_SRC, shim_dir)
    for f in shim_dir.iterdir():
        f.chmod(f.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    cfg = {
        "sample": {"sample_id": "smoke"},
        "inputs": {
            "genome": str(fixture / "genome.fa"),
            "proteins": str(fixture / "proteins.faa"),
            "gff3": str(fixture / "annotation.gff3"),
            "eggnog": str(fixture / "emapper.annotations"),
        },
        "databases": {
            "pfam": str(fixture / "db" / "pfam.hmm"),
            "nbd_hmms": str(fixture / "db" / "nbd.hmm"),
            "custom_hmms": str(fixture / "db" / "custom.hmm"),
        },
        "execution": {"output_dir": str(run / "out"), "cpu_threads": 2},
        "discovery": {"discovery_mode": "BROAD"},
        "asm": {"asm_enable": 0},
        "fusion": {"enabled": 0},
        "reporting": {"plots": 0},
    }
    cfg_path = run / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    env = dict(os.environ)
    env["PATH"] = f"{shim_dir}{os.pathsep}{Path(sys.executable).parent}{os.pathsep}{env['PATH']}"
    env["FUNLR_TEST_DOUBLES"] = "1"
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["FUNLR_FIXTURE_DIR"] = str(fixture)
    proc = subprocess.run(
        [sys.executable, "-m", "funlr", "run", "--config", str(cfg_path)],
        capture_output=True, text=True, env=env, cwd=REPO,
    )
    (run / "stdout.log").write_text(proc.stdout)
    (run / "stderr.log").write_text(proc.stderr)
    assert proc.returncode == 0, (
        f"funlr run failed ({proc.returncode}):\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )
    return run / "out"


def _tsv(run_dir, relpath):
    p = run_dir / relpath
    assert p.exists(), f"missing {p}"
    return pd.read_csv(p, sep="\t")


def test_stage_dirs_created(run_dir):
    for n in range(8):
        assert (run_dir / "results" / f"stage{n}").is_dir()


def test_tier_assignments(run_dir):
    tiers = _tsv(run_dir, "results/stage3/tables/tiered_candidates.tsv")
    got = dict(zip(tiers["protein_id"], tiers["tier"]))
    assert got == EXPECTED_TIERS


def test_expected_candidate_set(run_dir):
    tiers = _tsv(run_dir, "results/stage3/tables/tiered_candidates.tsv")
    # All 12 fixture proteins enter the union (FUN_000010 via the legacy
    # HET-keyword quirk matching "synthetase").
    assert set(tiers["protein_id"]) == set(EXPECTED_TIERS)


def test_flags_present(run_dir):
    tiers = _tsv(run_dir, "results/stage3/tables/tiered_candidates.tsv").set_index("protein_id")

    def flags(pid):
        v = tiers.loc[pid, "flags"]
        return set() if pd.isna(v) or not v else set(str(v).split(";"))

    assert "HARD_END" in flags("FUN_000005-T1")
    assert "HARD_END" in flags("FUN_000002-T1")
    assert "NBD_ONLY_UNKNOWN_SENSOR" in flags("FUN_000002-T1")
    assert "NON_NLR_ANNOT" in flags("FUN_000003-T1")
    assert "ORIENTATION_INVALID" in flags("FUN_000011-T1")
    assert "FP_DOMAIN_PRESENT" in flags("FUN_000012-T1")
    assert "VERY_SHORT_PROTEIN" in flags("FUN_000007-T1")


def test_rescue_score_bounds(run_dir):
    tiers = _tsv(run_dir, "results/stage3/tables/tiered_candidates.tsv")
    assert tiers["rescue_priority"].between(0, 100).all()


def test_stage4_exports(run_dir):
    s4 = run_dir / "results" / "stage4"
    fa = (s4 / "TIER_1A_HIGH_CONFIDENCE.faa").read_text()
    assert ">FUN_000001-T1" in fa
    assert ">FUN_000011-T1" not in fa  # 011 is excluded for reversed orientation
    fa2a = (s4 / "TIER_2A_HIGH_PRIORITY_RESCUE.faa").read_text()
    assert ">FUN_000002-T1" in fa2a
    assert (s4 / "tiered_candidates.faa").read_text().count(">") == 12
    assert (s4 / "fasta_export_manifest.tsv").exists()


def test_miniprot_outputs(run_dir):
    s5 = run_dir / "results" / "stage5"
    assert (s5 / "gff" / "miniprot.gff3").exists()
    summ = _tsv(run_dir, "results/stage5/tables/miniprot_summary.tsv")
    # canned miniprot resolves the 2A and one 1B query
    assert set(summ["query"]) == {"FUN_000002-T1", "FUN_000005-T1"}
    comprehensive = _tsv(run_dir, "results/stage5/comprehensive/tables/miniprot_summary.tsv")
    assert set(comprehensive["query"]) == {"FUN_000002-T1", "FUN_000005-T1"}
    # unresolved queries land on the exonerate refine list
    refine = (s5 / "queries" / "exonerate_refine_ids.txt").read_text().split()
    assert refine == ["FUN_000009-T1"]


def test_exonerate_refine(run_dir):
    summ = _tsv(run_dir, "results/stage6/tables/exonerate_summary.tsv")
    assert set(summ["protein_id"]) == {"FUN_000009-T1"}
    comprehensive = _tsv(run_dir, "results/stage6/comprehensive/tables/exonerate_summary.tsv")
    assert set(comprehensive["protein_id"]) == {"FUN_000009-T1"}


def test_final_report_and_strict(run_dir):
    s7 = run_dir / "results" / "stage7" / "tables"
    report = pd.read_csv(s7 / "nlr_final_report.tsv", sep="\t")
    assert set(report["protein_id"]) == set(EXPECTED_TIERS)
    strict = pd.read_csv(s7 / "nlr_strict_candidates.tsv", sep="\t")
    assert set(strict["protein_id"]) == EXPECTED_STRICT
    priority = pd.read_csv(s7 / "nlr_strict_candidates.priority.tsv", sep="\t")
    assert set(priority["protein_id"]) == {"FUN_000002-T1", "FUN_000005-T1", "FUN_000009-T1"}
    comprehensive = pd.read_csv(s7 / "nlr_strict_candidates.comprehensive.tsv", sep="\t")
    assert set(comprehensive["protein_id"]) == {"FUN_000002-T1", "FUN_000005-T1"}
    assert report["is_fusion_model"].eq(False).all()
    assert {"NLR_architecture", "nbd_confidence", "comprehensive_mp_scaffold"}.issubset(report)


def test_final_results_mirror(run_dir):
    fr = run_dir / "final_results"
    for name in ("nlr_final_report.tsv", "nlr_strict_candidates.tsv",
                 "tier_summary.tsv", "qc_summary.tsv", "nlr_candidates.bed"):
        assert (fr / name).exists(), f"missing final_results/{name}"


def test_provenance_and_state(run_dir):
    assert (run_dir / "logs" / "commands.jsonl").exists()
    lines = (run_dir / "logs" / "commands.jsonl").read_text().strip().splitlines()
    assert len(lines) >= 6  # faidx, hmmpress/hmmscan xN, seqkit xN, miniprot, exonerate
    manifest = json.loads((run_dir / "manifest.json").read_text())
    state = json.loads((run_dir / "state.json").read_text())
    assert manifest["status"] == "complete"
    assert set(state["stages"]) == set(map(str, range(8)))
    assert all(stage["status"] == "complete" for stage in state["stages"].values())
    assert manifest["code"] and manifest["tools"] and manifest["runtime"]
    assert manifest["signature"] == state["signature"]


def _repeat_run(run_dir, *args):
    env = dict(os.environ)
    env.update(FUNLR_TEST_DOUBLES="1", FUNLR_FIXTURE_DIR=str(run_dir.parent / "fixture"),
               PYTHONPATH=str(REPO / "src"),
               PATH=os.pathsep.join((str(run_dir.parent / "shims"),
                                    str(Path(sys.executable).parent), env["PATH"])))
    return subprocess.run(
        [sys.executable, "-m", "funlr", "run", "--config", str(run_dir.parent / "config.yaml"), *args],
        capture_output=True, text=True, env=env, cwd=REPO,
    )


def test_resume_verifies_every_stage_without_running_commands(run_dir):
    state_before = (run_dir / "state.json").read_bytes()
    commands_before = (run_dir / "logs" / "commands.jsonl").read_bytes()
    proc = _repeat_run(run_dir, "--resume")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (run_dir / "state.json").read_bytes() == state_before
    assert (run_dir / "logs" / "commands.jsonl").read_bytes() == commands_before


@pytest.mark.parametrize("flag", ["--dry-run", "--dry_run"])
def test_dry_run_writes_nothing_and_preserves_completed_run(run_dir, tmp_path, flag):
    def snapshot():
        return {str(p.relative_to(run_dir)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in run_dir.rglob("*") if p.is_file()}

    before = snapshot()
    proc = _repeat_run(run_dir, flag)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert snapshot() == before
    fresh = tmp_path / "planned-output"
    proc = _repeat_run(run_dir, flag, "--outdir", str(fresh))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not fresh.exists()
