"""Configuration merge regressions: readable YAML must not silently change a run."""
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from funlr.cli import parser, resolve_args
from funlr.config import DEFAULTS, FunlrConfig, write_example_yaml


def test_starter_yaml_preserves_original_scientific_defaults():
    config = FunlrConfig(yaml.safe_load(write_example_yaml()))
    assert config.scientific_settings() == DEFAULTS


def test_relative_yaml_paths_and_explicit_cli_override(tmp_path):
    config = tmp_path / "settings.yaml"
    config.write_text("inputs:\n  genome: data/genome.fa\nexecution:\n  output_dir: results\n  cpu_threads: 7\n  dry_run: true\nthresholds:\n  eval_nbd: 1e-5\n")
    args = parser().parse_args(["run", "--config", str(config), "--threads", "3", "--proteins", "cli.faa"])
    resolved = resolve_args(args)
    assert args.genome == str(tmp_path / "data/genome.fa")
    assert args.outdir == str(tmp_path / "results")
    assert args.proteins == "cli.faa"
    assert args.threads == 3
    assert args.dry_run is True
    assert float(resolved["EVAL_NBD"]) == 1e-5


@pytest.mark.parametrize("data", [
    {"thresholds": {"eval_ndb": 1e-4}},
    {"execution": {"cpu_threads": True}},
    {"execution": {"cpu_threads": 0}},
    {"execution": {"dry_run": "false"}},
    {"tools": {"hmmscan": None}},
    {"inputs": ["genome.fa"]},
    {"inputs": {"genome": "first\nsecond"}},
])
def test_invalid_nested_config_rejected(data):
    with pytest.raises(ValueError):
        FunlrConfig(data)


def test_yaml_object_constructors_are_rejected(tmp_path):
    path = tmp_path / "invalid.yaml"
    path.write_text("!!python/object/apply:builtins.str [unsafe]\n")
    with pytest.raises(ValueError, match="Invalid YAML"):
        FunlrConfig.from_yaml(path)


def test_hifi_starter_and_explicit_override_order():
    cfg = FunlrConfig(yaml.safe_load(write_example_yaml('HIFI')))
    values = cfg.scientific_settings()
    assert values['HARD_END_BP'] == values['CONTIG_END_BP'] == 1000
    assert values['NBD_CONF_COV_HIGH'] == 0.50
    cfg = FunlrConfig({'context': {'profile': 'HIFI', 'hard_end_bp': 2000}})
    assert cfg.scientific_settings()['CONTIG_END_BP'] == 2000


@pytest.mark.parametrize('old,new,value', [('CONTIG_END_BP', 'HARD_END_BP', 500), ('INCLUDE_TIER3', 'INCLUDE_TIER3_IN_RESCUE', 1)])
def test_legacy_alias_updates_active_setting_and_rejects_conflicts(old, new, value):
    settings = FunlrConfig({old: value}).scientific_settings()
    assert settings[new] == settings[old] == value
    with pytest.raises(ValueError, match='Conflicting'):
        FunlrConfig({old: value, new: DEFAULTS[new]})


def test_unused_settings_cannot_silently_ignore_changes():
    from funlr.updated_defaults import UNUSED_LEGACY
    for key in UNUSED_LEGACY:
        value = DEFAULTS[key]
        changed = value + '_CHANGED' if isinstance(value, str) else (1 - value if value in (0, 1) else value + 1)
        with pytest.raises(ValueError, match='not implemented by the analysis workflow'):
            FunlrConfig({key: changed})


def test_unused_exonerate_retry_is_only_recorded_profile_default():
    assert FunlrConfig().scientific_settings()['EXON_MAX_INTRON_RETRY'] == 0
    assert FunlrConfig({'context': {'profile': 'HIFI'}}).scientific_settings()['EXON_MAX_INTRON_RETRY'] == 10000
    with pytest.raises(ValueError, match='not implemented by the analysis workflow'):
        FunlrConfig({'rescue': {'exon_max_intron_retry': 5000}})
