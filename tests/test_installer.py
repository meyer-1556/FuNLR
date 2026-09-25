"""Installer boundaries: no network in dry-run, integrity before execution."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("funlr_installer", ROOT / "install.py")
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


def test_dry_run_does_not_write_or_download(tmp_path):
    destination = tmp_path / "new env"
    before = list(tmp_path.rglob("*"))
    with patch.object(installer, "host_platform", return_value="linux-64"), \
         patch.object(installer, "urlopen", side_effect=AssertionError("network")), \
         patch.object(installer.subprocess, "run", side_effect=AssertionError("execution")):
        result = installer.install(destination, demo=True, dry_run=True)
    assert result["status"] == "PLANNED"
    assert result["packages"] == 103
    assert list(tmp_path.rglob("*")) == before


def test_existing_prefix_never_overwritten(tmp_path):
    sentinel = tmp_path / "keep"
    sentinel.write_text("original")
    with patch.object(installer, "host_platform", return_value="linux-64"):
        with pytest.raises(ValueError, match="new prefix"):
            installer.install(tmp_path, dry_run=True)
    assert sentinel.read_text() == "original"


def test_tampered_cached_download_is_rejected(tmp_path):
    cached = tmp_path / "download.tar.bz2"
    cached.write_bytes(b"altered")
    with patch.object(installer, "urlopen", side_effect=AssertionError("network")):
        with pytest.raises(ValueError, match="checksum mismatch"):
            installer.download_verified("https://example.invalid", cached, hashlib.sha256(b"original").hexdigest())


def test_bootstrap_extracts_only_expected_regular_member(tmp_path):
    archive, executable = tmp_path / "archive.tar.bz2", tmp_path / "micromamba"
    with tarfile.open(archive, "w:bz2") as out:
        for name, data in [("../../escaped", b"bad"), ("bin/micromamba", b"known executable")]:
            member = tarfile.TarInfo(name); member.size = len(data)
            out.addfile(member, io.BytesIO(data))
    installer.extract_micromamba(archive, executable)
    assert executable.read_bytes() == b"known executable"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["archive.tar.bz2", "micromamba"]


def test_bootstrap_rejects_link_instead_of_executable(tmp_path):
    archive = tmp_path / "archive.tar.bz2"
    with tarfile.open(archive, "w:bz2") as out:
        member = tarfile.TarInfo("bin/micromamba"); member.type = tarfile.SYMTYPE; member.linkname = "/bin/sh"
        out.addfile(member)
    with pytest.raises(ValueError, match="not a regular"):
        installer.extract_micromamba(archive, tmp_path / "micromamba")


def test_both_locks_agree_with_package_manifest():
    metadata = json.loads((ROOT / "environments/locks.json").read_text())
    for platform_name in metadata["platforms"]:
        packages = installer.validate_lock(ROOT / "environments" / f"conda-{platform_name}.lock", platform_name)
        assert len({p["name"] for p in packages}) == len(packages)
        assert all(len(p["sha256"]) == 64 and p["url"].startswith("https://") for p in packages)
        assert {"hmmer", "samtools", "seqkit", "miniprot", "exonerate", "python", "pyyaml"} <= {p["name"] for p in packages}


def test_lock_tampering_is_reported(tmp_path):
    lock = tmp_path / "changed.lock"
    lock.write_text("@EXPLICIT\nhttps://example.invalid/changed.tar.bz2#0123\n")
    with pytest.raises(ValueError, match="disagree"):
        installer.validate_lock(lock, "linux-64")


def test_source_nested_prefix_is_rejected_without_writes():
    prefix = ROOT / "src" / "must-not-be-created"
    with patch.object(installer, "host_platform", return_value="linux-64"):
        with pytest.raises(ValueError, match="outside the source"):
            installer.install(prefix, dry_run=True)
    assert not prefix.exists()


@pytest.mark.parametrize("prefix", ["line\nbreak", "null\0byte", "carriage\rreturn"])
def test_control_character_prefix_is_rejected(prefix):
    with patch.object(installer, "host_platform", return_value="linux-64"):
        with pytest.raises(ValueError, match="control characters"):
            installer.install(prefix, dry_run=True)


def test_inherited_pip_settings_cannot_redirect_install(tmp_path):
    with patch.dict(installer.os.environ, {"PIP_TARGET": "/unrelated", "PIP_PREFIX": "/elsewhere", "PIP_USER": "1", "CONDA_SUBDIR": "wrong-arch", "PYTHONPATH": "/foreign"}):
        env = installer.installation_environment(tmp_path)
    assert not {"PIP_TARGET", "PIP_PREFIX", "PIP_USER", "CONDA_SUBDIR", "PYTHONPATH"}.intersection(env)
    assert env["PIP_CONFIG_FILE"] == installer.os.devnull


@pytest.mark.parametrize("system,machine,libc", [("Windows", "AMD64", ("", "")), ("Linux", "aarch64", ("glibc", "2.35")), ("Linux", "x86_64", ("musl", "1.2")), ("Linux", "x86_64", ("glibc", "2.16"))])
def test_unsupported_hosts_fail_early(system, machine, libc):
    with patch.object(installer.platform, "system", return_value=system), \
         patch.object(installer.platform, "machine", return_value=machine), \
         patch.object(installer.platform, "libc_ver", return_value=libc):
        with pytest.raises(ValueError):
            installer.host_platform()
