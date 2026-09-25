"""Public CLI boundaries for templates and read-only inspection."""
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]


def command(*args, cwd):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), PYTHONDONTWRITEBYTECODE="1")
    return subprocess.run([sys.executable, "-m", "funlr", *args], cwd=cwd, env=env,
                          text=True, capture_output=True, timeout=30)


def test_template_commands_write_only_stdout(tmp_path):
    result = command("init-samples", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    header = result.stdout.strip().split("\t")
    assert {"SAMPLE_ID", "SPECIES_ID", "GENOME_PATH", "PROTEINS_PATH", "GFF3_PATH"} <= set(header)
    assert len(header) == len(set(header))
    result = command("init-config", "--minimal", "--profile", "HIFI", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    config = yaml.safe_load(result.stdout)
    assert config["context"]["profile"] == "HIFI"
    assert "thresholds" not in config
    assert not list(tmp_path.iterdir())


def test_missing_run_is_a_concise_cli_error_without_writes(tmp_path):
    result = command("status", "--run-dir", str(tmp_path / "missing"), "--verify", cwd=tmp_path)
    assert result.returncode == 1
    assert "Traceback" not in result.stderr
    assert not list(tmp_path.iterdir())


def test_unsupported_original_readme_columns_fail_before_output_creation(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("SAMPLE_ID\tSPECIES_ID\tFUNLR_PROFILE\nexample\tFungus_species\tFungi\n")
    outdir = tmp_path / "cohort"
    result = command("batch", "--samples", str(sheet), "--outdir", str(outdir), "--dry-run", cwd=tmp_path)
    assert result.returncode == 1
    assert "FUNLR_PROFILE" in result.stderr
    assert not outdir.exists()
