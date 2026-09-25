"""Read-only profile inventory, independent of discovery and HMMER execution."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from importlib.resources import files
from pathlib import Path


def _metadata() -> dict:
    return json.loads(files("funlr").joinpath("data/model_sources.json").read_text())


def _profile(lines: list[bytes], ordinal: int, library: Path) -> dict:
    headers = {}
    for line in lines[1:]:
        if line.startswith(b"HMM ") or line.startswith(b"HMM\t"):
            break
        fields = line.split(None, 1)
        if fields and fields[0] in {b"NAME", b"ACC", b"LENG", b"NSEQ"}:
            key = fields[0].decode("ascii")
            if key in headers or len(fields) != 2:
                raise ValueError(f"Invalid {key} header in {library.name}, profile {ordinal}")
            try:
                headers[key] = fields[1].decode("utf-8").strip()
            except UnicodeDecodeError as exc:
                raise ValueError(f"Invalid text header in {library.name}, profile {ordinal}") from exc
    if not headers.get("NAME") or not headers.get("LENG"):
        raise ValueError(f"Missing NAME/LENG in {library.name}, profile {ordinal}")
    for key in ("LENG", "NSEQ"):
        if key in headers:
            try:
                value = int(headers[key])
            except ValueError as exc:
                raise ValueError(f"Invalid {key} in {library.name}, profile {ordinal}") from exc
            if value < 1:
                raise ValueError(f"Invalid {key} in {library.name}, profile {ordinal}")
            headers[key] = value
    body = b"".join(lines)
    without_name = b"".join(line for line in lines if line.split(None, 1)[:1] != [b"NAME"])
    return {
        "ordinal": ordinal,
        "name": headers["NAME"],
        "accession": headers.get("ACC"),
        "length": headers["LENG"],
        "nseq": headers.get("NSEQ"),
        "bytes": len(body),
        "sha256": hashlib.sha256(body).hexdigest(),
        "sha256_without_name": hashlib.sha256(without_name).hexdigest(),
        # A library-level source audit does not identify each publication profile.
        "publication_profile_identity": None,
        "construction_provenance": None,
    }


def _library(path: Path, metadata: dict) -> dict:
    digest = hashlib.sha256()
    size = 0
    records = []
    active: list[bytes] = []
    with path.open("rb") as stream:
        for lineno, line in enumerate(stream, 1):
            digest.update(line)
            size += len(line)
            if line.startswith(b"HMMER3/"):
                if active:
                    raise ValueError(f"Unterminated HMM profile: {path.name}:{lineno}")
                active = [line]
            elif active:
                active.append(line)
                if line.strip() == b"//":
                    records.append(_profile(active, len(records) + 1, path))
                    active = []
            elif line.strip():
                raise ValueError(f"Unexpected text outside HMM profile: {path.name}:{lineno}")
    if active or not records:
        raise ValueError(f"Empty or unterminated HMM library: {path.name}")
    source = next((row for row in metadata["libraries"]
                   if row["file"] == path.name and row["sha256"] == digest.hexdigest()), None)
    groups = defaultdict(list)
    for row in records:
        groups[row["sha256_without_name"]].append({"ordinal": row["ordinal"], "name": row["name"]})
    duplicates = [{"sha256_without_name": key, "profiles": value}
                  for key, value in sorted(groups.items()) if len(value) > 1]
    return {
        "file": path.name, "bytes": size, "sha256": digest.hexdigest(),
        "profile_count": len(records), "profiles": records,
        "snapshot_identity": "RECOGNIZED" if source else "UNRECOGNIZED",
        "role": source["role"] if source else None,
        "source_audit": source["source_audit"] if source else None,
        "citations": source["citations"] if source else None,
        "construction_provenance": None,
        "redistribution_permission": None,
        "identical_except_name_groups": duplicates,
        "identical_except_name_extra_records": sum(len(row["profiles"]) - 1 for row in duplicates),
    }


def inventory_models(directory, verify_snapshot: bool = False) -> dict:
    """Inventory direct *.hmm files without modifying models or running tools.

    Header/record structure is checked, not HMM parameter validity. Source
    annotations attach only to known filename-and-SHA256 pairs. Duplicate
    profiles remain in order; this command never deduplicates a database.
    """
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"Model directory does not exist: {directory}")
    paths = sorted(directory.glob("*.hmm"), key=lambda path: path.name)
    if not paths:
        raise ValueError(f"No *.hmm libraries found: {directory}")
    for path in paths:
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Model library must be a regular, non-symlink file: {path.name}")
    snapshot = None
    if verify_snapshot:
        from .databases import verify_models
        snapshot = verify_models(directory)
    metadata = _metadata()
    libraries = [_library(path, metadata) for path in paths]
    if snapshot is not None:
        by_name = {row["file"]: row for row in libraries}
        for checked in snapshot["libraries"]:
            current = by_name.get(checked["file"])
            if current is None or any(current[key] != checked[key] for key in ("sha256", "bytes")) or current["profile_count"] != checked["profiles"]:
                raise ValueError(f"Model library changed during snapshot inventory: {checked['file']}")
    return {
        "schema_version": 1, "status": "INVENTORIED",
        "source_metadata_version": metadata["metadata_version"],
        "directory": str(directory),
        "snapshot_verification": snapshot,
        "scope": "HMM header and byte inventory; library source audits do not prove per-profile training history, biological specificity, or redistribution permission.",
        "libraries": libraries,
    }
