"""Installation checks and demo failure boundaries; real tools are opt-in."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from funlr import demo, doctor


def test_doctor_missing_tools_is_read_only_and_lists_every_tool(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PATH", "")
    report = doctor.check_installation(include_demo=True)
    assert not report["ok"]
    assert set(report["tools"]) == {*doctor.TOOLS, "hmmbuild"}
    assert all(not item["ok"] for item in report["tools"].values())
    assert all("Executable not found" in item["error"] for item in report["tools"].values())
    assert list(tmp_path.iterdir()) == []
    assert report["python"]["executable"] == sys.executable


def test_doctor_reports_failed_executable_instead_of_calling_it_installed(tmp_path):
    broken = tmp_path / "broken-tool"
    broken.write_text(f"#!{sys.executable}\nimport sys\nprint('incompatible executable')\nsys.exit(9)\n")
    broken.chmod(0o755)
    report = doctor.check_installation(configured_tools={key: str(broken) for key in doctor.TOOLS})
    assert not report["ok"]
    assert all("exited 9" in entry["error"] for entry in report["tools"].values())
    assert "hmmbuild" not in report["tools"]


def test_doctor_captures_python_dependency_probe_failure(monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("cannot launch interpreter")
    monkeypatch.setattr(doctor.subprocess, "run", fail)
    result = doctor._python_dependencies()
    assert set(result) == {"pandas", "numpy", "PyYAML", "matplotlib"}
    assert all(not item["ok"] and "cannot launch interpreter" in item["error"] for item in result.values())


def test_demo_missing_dependencies_fail_before_creating_output(tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "check_installation", lambda **kwargs: {
        "ok": False, "errors": ["hmmbuild: executable not found"]})
    root = tmp_path / "new-demo"
    with pytest.raises(ValueError, match="hmmbuild"):
        demo.run_demo(root)
    assert not root.exists()


def test_demo_preserves_existing_directory_without_probing_or_writing(tmp_path, monkeypatch):
    marker = tmp_path / "keep.txt"
    marker.write_text("original")
    def unexpected(**kwargs):
        pytest.fail("existing output must be rejected before dependency checks")
    monkeypatch.setattr(demo, "check_installation", unexpected)
    with pytest.raises(ValueError, match="new directory"):
        demo.run_demo(tmp_path)
    assert list(tmp_path.iterdir()) == [marker]
    assert marker.read_text() == "original"


def test_demo_rejects_broken_symlink_output(tmp_path, monkeypatch):
    output = tmp_path / "broken"
    output.symlink_to(tmp_path / "nonexistent")
    with pytest.raises(ValueError, match="new directory"):
        demo.run_demo(output)
    assert output.is_symlink()
    assert not (tmp_path / "nonexistent").exists()


def test_demo_build_failure_records_failure_evidence(tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "check_installation", lambda **kwargs: {
        "ok": True, "errors": [], "tools": {"hmmbuild": {"path": "/missing/hmmbuild"}}})
    def fail(command, logfile, env, cwd, command_log):
        raise RuntimeError(f"hmmbuild failed; see {logfile}")
    monkeypatch.setattr(demo, "invoke", fail)
    root = tmp_path / "failed-demo"
    with pytest.raises(RuntimeError, match="hmmbuild failed"):
        demo.run_demo(root)
    evidence = json.loads((root / "evidence.json").read_text())
    assert evidence["status"] == "failed"
    assert "hmmbuild failed" in evidence["error"]
    assert "NOT biological validation" in evidence["purpose"]
    assert "completed_stages" not in evidence


def test_bundled_demo_inputs_are_deterministic_and_self_contained(tmp_path):
    first = tmp_path / "a"
    second = tmp_path / "b"
    demo._write_inputs(first)
    demo._write_inputs(second)
    files = {path.name for path in first.iterdir()}
    assert files == {"genome.fa", "proteins.faa", "genes.gff3", "NACHT.afa", "WD40.afa"}
    assert all((first / name).read_bytes() == (second / name).read_bytes() for name in files)
    assert (first / "proteins.faa").read_text().count(">") == 3
    assert "ID=synthetic_rescue;Parent=gene1" in (first / "genes.gff3").read_text()


def test_demo_cli_help_works_outside_repository(tmp_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(demo.__file__).resolve().parent.parent)
    result = subprocess.run([sys.executable, "-m", "funlr", "demo", "--help"],
                            cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--outdir" in result.stdout
    assert not list(tmp_path.iterdir())


@pytest.mark.skipif(os.environ.get("FUNLR_TEST_REAL_TOOLS") != "1", reason="opt-in real bioinformatics executables")
def test_real_bundled_demo_runs_all_stages_and_verifies_resume(tmp_path):
    evidence = demo.run_demo(tmp_path / "installation demo with spaces")
    assert evidence["status"] == "passed"
    assert evidence["completed_stages"] == 8
    assert evidence["resume_reused_stages"] == 8
    assert any(row["protein_id"] == "synthetic_rescue" for row in evidence["positive_exonerate_rows"])
    assert Path(evidence["evidence_file"]).is_file()


def _mock_healthy_installation(monkeypatch, dependency_versions=None):
    versions = dependency_versions or doctor._EXPECTED_VERSIONS
    monkeypatch.setattr(doctor, "tool_inventory", lambda configured, scientific=None: {
        name: {"path": f"/test/{name}", "version_output": "test version"} for name in doctor.TOOLS})
    payload = {name: {"ok": True, "version": version} for name, version in versions.items()}
    monkeypatch.setattr(doctor.subprocess, "run", lambda *args, **kwargs:
                        subprocess.CompletedProcess(args[0], 0, json.dumps(payload), ""))


@pytest.mark.parametrize("version,expected_ok", [((3, 10, 14), False), ((3, 11, 10), True),
                                               ((3, 12, 14), True), ((3, 13, 0), False)])
def test_doctor_enforces_supported_python_range(monkeypatch, version, expected_ok):
    _mock_healthy_installation(monkeypatch)
    monkeypatch.setattr(doctor.sys, "version_info", version)
    report = doctor.check_installation()
    assert report["python"]["ok"] is expected_ok
    assert report["ok"] is expected_ok
    assert report["python"]["requires"] == ">=3.11,<3.13"
    if not expected_ok:
        assert "Python >=3.11,<3.13 is required" in report["errors"]


@pytest.mark.parametrize("package,wrong_version", [("pandas", "2.3.1"), ("numpy", "2.0.0"), ("PyYAML", "6.0.2"), ("matplotlib", "3.9.0")])
def test_doctor_rejects_importable_but_wrong_dependency_version(monkeypatch, package, wrong_version):
    versions = {**doctor._EXPECTED_VERSIONS, package: wrong_version}
    _mock_healthy_installation(monkeypatch, versions)
    report = doctor.check_installation()
    assert not report["ok"]
    entry = report["python_dependencies"][package]
    assert entry["version"] == wrong_version
    assert entry["expected_version"] == doctor._EXPECTED_VERSIONS[package]
    assert entry["ok"] is False
    assert f"found {wrong_version}" in entry["error"]
    assert len(report["errors"]) == 1


def test_doctor_dependency_contract_matches_package_metadata():
    import tomllib
    project = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())["project"]
    assert project["requires-python"] == ">=3.11,<3.13"
    assert set(project["dependencies"]) == {f"{package}=={version}" for package, version in doctor._EXPECTED_VERSIONS.items()}
