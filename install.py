#!/usr/bin/env python3
"""Install a locked FuNLR environment without modifying shell startup files.

Bootstrap requires Python 3.9+ and internet access. The installed application uses
its own pinned Python. Run `python3 install.py --prefix /path/to/new/env --demo`.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
VERSION = "0.5.0a1"
MAMBA_VERSION = "2.9.0"
BOOTSTRAP = {
    "linux-64": {
        "url": "https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.9.0-0.tar.bz2",
        "sha256": "8761c382127e6363bd9e0a2451aa3ef90d071a79133f736e2f759a3bf13040dd",
    },
    "osx-arm64": {
        "url": "https://conda.anaconda.org/conda-forge/osx-arm64/micromamba-2.9.0-0.tar.bz2",
        "sha256": "500f5074feb8d02c4296ef9921c3650ed2874171805a9fbb8fbb53896433646b",
    },
}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def host_platform():
    system, machine = platform.system(), platform.machine().lower()
    if system == "Linux" and machine in ("x86_64", "amd64"):
        libc, version = platform.libc_ver()
        if libc != "glibc" or not version:
            raise ValueError("The Linux lock requires glibc >=2.28. Use the Linux container on other Linux distributions.")
        if tuple(int(v) for v in version.split(".")[:2]) < (2, 28):
            raise ValueError("This environment requires glibc >=2.28; use a supported node or the Linux container.")
        return "linux-64"
    if system == "Darwin" and machine in ("arm64", "aarch64"):
        return "osx-arm64"
    raise ValueError(f"No tested dependency lock for {system}/{machine}. Supported locks: Linux x86_64 and macOS Apple Silicon. Windows users can use Linux x86_64 under WSL2; see docs/INSTALLATION.md.")


def validate_lock(lock, platform_name):
    metadata = json.loads((ROOT / "environments/locks.json").read_text())
    packages = metadata["platforms"][platform_name]["packages"]
    lines = [line.strip() for line in lock.read_text().splitlines() if line.strip() and not line.startswith("#")]
    expected = [p["url"] + "#" + p["md5"] for p in packages]
    if lines != ["@EXPLICIT", *expected]:
        raise ValueError(f"Lock file and package manifest disagree: {lock}")
    return packages


def download_verified(url, target, expected_sha256):
    if target.exists():
        if target.is_symlink() or digest(target) != expected_sha256:
            raise ValueError(f"Cached bootstrap checksum mismatch: {target}. Remove that cache file before retrying.")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as out:
            temporary = Path(out.name)
            with urlopen(url, timeout=90) as response:
                if not response.geturl().startswith("https://"):
                    raise ValueError("Bootstrap redirected away from HTTPS")
                size = 0
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    size += len(chunk)
                    if size > 100 * 1024 * 1024:
                        raise ValueError("Bootstrap archive exceeds expected size limit")
                    out.write(chunk)
        if digest(temporary) != expected_sha256:
            raise ValueError("Downloaded micromamba checksum mismatch; executable was not installed")
        temporary.replace(target)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def extract_micromamba(archive, executable):
    """Copy only the known regular executable; never extract archive paths."""
    with tarfile.open(archive, "r:bz2") as package:
        member = package.getmember("bin/micromamba")
        if not member.isfile():
            raise ValueError("Bootstrap executable is not a regular archive member")
        source = package.extractfile(member)
        if source is None:
            raise ValueError("Bootstrap executable is missing")
        with source, tempfile.NamedTemporaryFile(dir=executable.parent, delete=False) as out:
            temporary = Path(out.name)
            shutil.copyfileobj(source, out)
    try:
        temporary.chmod(0o755)
        temporary.replace(executable)
    finally:
        if temporary.exists():
            temporary.unlink()


def command(args, *, env, cwd=None, capture=False):
    print("+ " + shlex.join([str(arg) for arg in args]), flush=True)
    result = subprocess.run([str(arg) for arg in args], env=env, cwd=cwd, text=True,
                            stdout=subprocess.PIPE if capture else None, check=False)
    if result.returncode and capture and result.stdout:
        print(result.stdout, file=sys.stderr)
    result.check_returncode()
    return result.stdout if capture else None


def installation_environment(cache):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("PIP_", "CONDA_", "MAMBA_"))}
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env.update(MAMBA_ROOT_PREFIX=str(cache / "mamba"), CONDA_PKGS_DIRS=str(cache / "mamba/pkgs"),
               XDG_CACHE_HOME=str(cache / "cache"), PIP_CONFIG_FILE=os.devnull, PIP_DISABLE_PIP_VERSION_CHECK="1", PYTHONNOUSERSITE="1")
    return env


def install(prefix, *, micromamba=None, demo=False, dry_run=False):
    platform_name = host_platform()
    if any(ord(char) < 32 or ord(char) == 127 for char in str(prefix)):
        raise ValueError("Installation prefix must not contain control characters")
    original = Path(prefix).expanduser()
    if original.is_symlink() or original.exists():
        raise ValueError(f"Installation requires a new prefix: {original}. Existing environments are never overwritten.")
    prefix = original.resolve()
    if prefix.is_relative_to((ROOT / "src").resolve()):
        raise ValueError("Installation prefix must be outside the source package directory (src)")
    lock = ROOT / "environments" / f"conda-{platform_name}.lock"
    packages = validate_lock(lock, platform_name)
    cache = prefix.parent / ".funlr-bootstrap" / platform_name
    demo_path = prefix.with_name(prefix.name + "-demo")
    if demo and (demo_path.exists() or demo_path.is_symlink()):
        raise ValueError(f"Demo needs a fresh output directory: {demo_path}")
    plan = {"status": "PLANNED", "version": VERSION, "platform": platform_name,
            "prefix": str(prefix), "lock": str(lock), "lock_sha256": digest(lock),
            "packages": len(packages), "package_download_bytes": sum(p["size"] for p in packages),
            "bootstrap": BOOTSTRAP[platform_name] if micromamba is None else str(micromamba),
            "demo_outdir": str(demo_path) if demo else None}
    print(json.dumps(plan, indent=2), flush=True)
    if dry_run:
        return plan
    for required in ("pyproject.toml", "README.md", "LICENSE", "NOTICE", "src/funlr/cli.py"):
        if not (ROOT / required).is_file():
            raise ValueError(f"Incomplete source archive: missing {required}")
    cache.mkdir(parents=True, exist_ok=True)
    env = installation_environment(cache)
    if micromamba is None:
        bootstrap = BOOTSTRAP[platform_name]
        archive = cache / f"micromamba-{MAMBA_VERSION}.tar.bz2"
        download_verified(bootstrap["url"], archive, bootstrap["sha256"])
        executable = cache / "micromamba"
        extract_micromamba(archive, executable)
    else:
        executable = Path(micromamba).expanduser().resolve()
    version = command([executable, "--version"], env=env, capture=True).strip()
    if version != MAMBA_VERSION:
        raise ValueError(f"Expected micromamba {MAMBA_VERSION}, found {version!r}")
    command([executable, "--no-rc", "create", "--yes", "--prefix", prefix,
             "--root-prefix", cache / "mamba", "--file", lock], env=env)
    env["PATH"] = str(prefix / "bin") + os.pathsep + env.get("PATH", "")
    # Build from a small temporary copy, so source can remain read-only and its
    # local caches/test data do not enter the distribution.
    with tempfile.TemporaryDirectory(prefix="funlr-source-", dir=cache) as staging:
        stage = Path(staging)
        for name in ("pyproject.toml", "README.md", "LICENSE", "NOTICE"):
            shutil.copy2(ROOT / name, stage / name)
        shutil.copytree(ROOT / "src", stage / "src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info"))
        command([prefix / "bin/python", "-m", "pip", "install", "--no-index", "--no-deps", "--no-build-isolation", "."], env=env, cwd=stage)
    report = json.loads(command([prefix / "bin/funlr", "doctor", "--demo"], env=env, cwd=prefix, capture=True))
    if not report["ok"]:
        raise RuntimeError("Installation completed but doctor failed; see the report above")
    if demo:
        command([prefix / "bin/funlr", "demo", "--outdir", demo_path], env=env, cwd=prefix)
    result = dict(plan, status="INSTALLED", doctor=report, bootstrap_version=version)
    (prefix / "funlr-installation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"\nInstalled FuNLR {VERSION}. For this terminal session:\n  export PATH={shlex.quote(str(prefix / 'bin'))}:\"$PATH\"\nThen run: funlr doctor\nNo shell startup files were changed.")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="funlr-env", help="New environment path (default: ./funlr-env)")
    parser.add_argument("--micromamba", help="Use an existing micromamba 2.9.0 executable")
    parser.add_argument("--demo", action="store_true", help="Run the real-tool synthetic installation demo")
    parser.add_argument("--dry-run", action="store_true", help="Print the installation plan without downloads or writes")
    args = parser.parse_args(argv)
    try:
        install(args.prefix, micromamba=args.micromamba, demo=args.demo, dry_run=args.dry_run)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError, tarfile.TarError) as exc:
        print(f"ERROR: {exc}\nA partial prefix is preserved for diagnosis; choose a new prefix when retrying.", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Installation interrupted; partial files are preserved.", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
