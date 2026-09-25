"""Installed-resource diagnostics and explicit CLI settings remain observable."""
from pathlib import Path

import pytest
import yaml

from funlr import cli, doctor, public_demo


def healthy(monkeypatch):
    calls = []
    def check(**kwargs):
        calls.append(kwargs)
        return {"ok": True, "errors": [], "tools": {}}
    monkeypatch.setattr(doctor, "check_installation", check)
    return calls


def test_self_test_verifies_real_package_resources_without_writes(tmp_path, monkeypatch):
    calls = healthy(monkeypatch)
    monkeypatch.chdir(tmp_path)
    result = doctor.check_self_test()
    assert result["ok"]
    assert result["public_fixture"]["status"] == "VERIFIED"
    assert result["public_fixture"]["files_checked"] == 10
    assert "db/pfam_mini.hmm" in result["public_fixture"]["sha256"]
    assert result["production_models"]["status"] == "NOT_REQUESTED"
    assert calls[0]["include_demo"] is True
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("error", [FileNotFoundError("fixture missing"),
                                    ValueError("checksum mismatch"), KeyError("sha256")])
def test_self_test_fails_for_missing_or_damaged_packaged_data(monkeypatch, error):
    healthy(monkeypatch)
    def broken():
        raise error
    monkeypatch.setattr(public_demo, "verified_assets", broken)
    result = doctor.check_self_test()
    assert not result["ok"]
    assert result["public_fixture"]["status"] == "FAILED"
    assert any("Packaged public fixture" in item for item in result["errors"])


def test_self_test_missing_requested_models_fails_without_creating_them(tmp_path, monkeypatch):
    healthy(monkeypatch)
    result = doctor.check_self_test(models_dir=tmp_path / "missing")
    assert not result["ok"]
    assert result["public_fixture"]["status"] == "VERIFIED"
    assert result["production_models"]["status"] == "FAILED"
    assert not list(tmp_path.iterdir())


def test_self_test_does_not_hide_runtime_failure(monkeypatch):
    monkeypatch.setattr(doctor, "check_installation", lambda **kwargs:
                        {"ok": False, "errors": ["exonerate: cannot execute"], "tools": {}})
    result = doctor.check_self_test()
    assert result["public_fixture"]["status"] == "VERIFIED"
    assert not result["ok"]
    assert result["errors"] == ["exonerate: cannot execute"]


@pytest.mark.parametrize("command", ["doctor", "validate"])
def test_installation_entrypoints_honor_config_tool_paths(tmp_path, monkeypatch, command, capsys):
    calls = healthy(monkeypatch)
    config = tmp_path / "settings.yaml"
    config.write_text(yaml.safe_dump({"tools": {"miniprot": "./local-miniprot"},
                                     "reporting": {"plot_backend": "r"}}))
    before = set(tmp_path.iterdir())
    assert cli.main([command, "--config", str(config), "--self-test"]) == 0
    assert calls[0]["configured_tools"]["miniprot"] == str(tmp_path / "local-miniprot")
    assert calls[0]["plot_backend"] == "r"
    assert '"self_test": true' in capsys.readouterr().out
    assert set(tmp_path.iterdir()) == before


@pytest.mark.parametrize("command", ["doctor", "validate"])
def test_model_check_requires_explicit_self_test(command, tmp_path, capsys):
    assert cli.main([command, "--models-dir", str(tmp_path)]) == 1
    assert "--models-dir requires --self-test" in capsys.readouterr().err


def test_self_test_does_not_silently_ignore_scientific_cli_inputs(capsys):
    assert cli.main(["validate", "--self-test", "--genome", "genome.fa"]) == 1
    assert "use plain validate" in capsys.readouterr().err


def test_doctor_explicit_backend_overrides_config(tmp_path, monkeypatch):
    calls = healthy(monkeypatch)
    config = tmp_path / "config.yaml"
    config.write_text("reporting:\n  plot_backend: r\n")
    assert cli.main(["doctor", "--config", str(config), "--plot-backend", "matplotlib"]) == 0
    assert calls[0]["plot_backend"] == "matplotlib"


def test_single_run_explicit_mode_and_profile_alias_preserve_written_threshold(tmp_path):
    config = tmp_path / "settings.yaml"
    config.write_text("context:\n  profile: ILLUMINA\ndiscovery:\n  discovery_mode: BROAD\nthresholds:\n  eval_pfam: 1e-8\n")
    args = cli.parser().parse_args(["run", "--config", str(config), "--nlr_profile", "HIFI",
                                   "--discovery_mode", "STRICT"])
    values = cli.resolve_args(args)
    assert values["NLR_PROFILE"] == "HIFI"
    assert values["DISCOVERY_MODE"] == "STRICT"
    assert float(values["EVAL_PFAM"]) == 1e-8


def test_fetch_dest_alias_uses_same_destination():
    args = cli.parser().parse_args(["fetch-pfam", "--dest", "db"])
    assert args.outdir == "db"


def test_ambient_variables_do_not_override_explicit_configuration(tmp_path, monkeypatch):
    monkeypatch.setenv("PFAM_DB", "/unexpected/pfam.hmm")
    monkeypatch.setenv("FUNLR_DB_DIR", "/unexpected/models")
    monkeypatch.setenv("FUNLR_DRY_RUN", "1")
    config = tmp_path / "settings.yaml"
    config.write_text("databases:\n  pfam: explicit.hmm\n")
    args = cli.parser().parse_args(["run", "--config", str(config)])
    cli.resolve_args(args)
    assert args.pfam == str(tmp_path / "explicit.hmm")
    assert args.nbd_hmms is None
    assert args.dry_run is False
