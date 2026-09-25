#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Build a separately installed, serial Exonerate 2.4.0 from verified source."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tarfile
import urllib.request


SOURCE_URL = (
    "https://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/"
    "exonerate-2.4.0.tar.gz"
)
SOURCE_SHA256 = "f849261dc7c97ef1f15f222e955b0d3daf994ec13c9db7766f1ac7e77baa4042"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(value: str) -> str:
    result = shutil.which(value)
    if result is None:
        raise SystemExit(f"Required build tool not found: {value}")
    return str(Path(result).resolve())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glib-prefix", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--install-prefix", required=True, type=Path)
    parser.add_argument("--pkg-config", default="pkg-config")
    parser.add_argument("--cc", default="cc")
    args = parser.parse_args()
    glib = args.glib_prefix.expanduser().resolve()
    work = args.work_dir.expanduser().resolve()
    prefix = args.install_prefix.expanduser().resolve()
    if work.exists() or prefix.exists():
        raise SystemExit("Use fresh, nonexistent work and installation directories.")
    if work == prefix or work in prefix.parents or prefix in work.parents:
        raise SystemExit("Work and installation directories must be separate.")
    if not (glib / "lib/pkgconfig/glib-2.0.pc").is_file():
        raise SystemExit(f"Missing GLib development metadata under {glib}")
    for path in (glib, work, prefix):
        if any(character.isspace() for character in str(path)):
            raise SystemExit("This upstream build system requires paths without whitespace.")
    cc = executable(args.cc)
    make = executable("make")
    patch_tool = executable("patch")
    pkg_config = executable(args.pkg_config)
    patch = Path(__file__).with_name("disable-pthreads-build.patch").resolve()
    if not patch.is_file():
        raise SystemExit(f"Missing companion patch: {patch}")

    work.mkdir(parents=True)
    archive = work / "exonerate-2.4.0.tar.gz"
    with urllib.request.urlopen(SOURCE_URL, timeout=120) as response:
        with archive.open("wb") as destination:
            shutil.copyfileobj(response, destination)
    if digest(archive) != SOURCE_SHA256:
        raise SystemExit("Downloaded source SHA256 mismatch; refusing to build.")
    with tarfile.open(archive, "r:gz") as source_archive:
        for member in source_archive.getmembers():
            destination = (work / member.name).resolve()
            if (
                not destination.is_relative_to(work)
                or member.issym()
                or member.islnk()
                or not (member.isfile() or member.isdir())
            ):
                raise SystemExit(f"Unexpected archive member: {member.name}")
        if hasattr(tarfile, "data_filter"):
            source_archive.extractall(work, filter="data")
        else:
            source_archive.extractall(work)
    source = work / "exonerate-2.4.0"
    wrapper_dir = work / "configure-bin"
    wrapper_dir.mkdir()
    wrapper = wrapper_dir / "pkg-config"
    wrapper.write_text(
        "#!/bin/sh\nexec "
        + shlex.quote(pkg_config)
        + " "
        + shlex.quote(f"--define-variable=prefix={glib}")
        + ' "$@"\n'
    )
    wrapper.chmod(0o755)
    env = os.environ.copy()
    env.update(
        PATH=str(wrapper_dir) + os.pathsep + os.environ.get("PATH", ""),
        PKG_CONFIG_LIBDIR=os.pathsep.join(
            (str(glib / "lib/pkgconfig"), str(glib / "share/pkgconfig"))
        ),
        PKG_CONFIG_PATH="",
        CC=cc,
        CFLAGS="-O2 -g -fno-strict-aliasing",
        LDFLAGS=f"-Wl,-rpath,{glib / 'lib'}",
    )
    manifest = {
        "source_url": SOURCE_URL,
        "source_sha256": SOURCE_SHA256,
        "patch_sha256": digest(patch),
        "patch_file": str(patch),
        "compiler": cc,
        "pkg_config": pkg_config,
        "glib_prefix": str(glib),
        "install_prefix": str(prefix),
        "environment": {key: env[key] for key in ("CC", "CFLAGS", "LDFLAGS")},
        "commands": [],
        "status": "building",
    }
    manifest_path = work / "build-provenance.json"

    def run(label: str, argv: list[str]) -> None:
        log_path = work / f"{label}.log"
        command = {"argv": argv, "log": str(log_path)}
        manifest["commands"].append(command)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"{label}: {shlex.join(argv)}", flush=True)
        with log_path.open("wb") as log:
            result = subprocess.run(argv, cwd=source, env=env, stdout=log, stderr=log)
        command["returncode"] = result.returncode
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        if result.returncode:
            manifest["status"] = "failed"
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
            raise SystemExit(f"{label} failed ({result.returncode}); inspect {log_path}")

    run("patch", [patch_tool, "--batch", "--forward", "-p1", "-i", str(patch)])
    run(
        "configure",
        [
            "./configure", f"--prefix={prefix}", "--disable-pthreads",
            "--disable-assert", "--disable-utilities",
        ],
    )
    # Upstream generated-header dependencies are unsafe with parallel make.
    run("make", [make, "-j1"])
    run("install", [make, "install"])
    binary = prefix / "bin/exonerate"
    manifest.update(status="built", binary=str(binary), binary_sha256=digest(binary))
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Built {binary}\nProvenance: {manifest_path}")
    print("Build success is not a scientific validation; run doctor and both demos.")


if __name__ == "__main__":
    main()
