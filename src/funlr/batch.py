"""Sequential cohorts using the ordinary, independently verified sample runner."""
from argparse import Namespace
from collections import Counter
import copy
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile

import yaml

from . import __version__
from .config import FunlrConfig
from .runner import now, run as run_sample
from .validation import INPUT_FLAGS, validate_inputs

PATH_COLUMNS = {
    "GENOME_PATH": ("inputs", "genome"),
    "PROTEINS_PATH": ("inputs", "proteins"),
    "GFF3_PATH": ("inputs", "gff3"),
    "EGGNOG_PATH": ("inputs", "eggnog"),
    "ANNOTATION_TABLE_PATH": ("inputs", "annotations"),
    "PFAM_PATH": ("databases", "pfam"),
    "NBD_HMMS": ("databases", "nbd_hmms"),
    "CUSTOM_HMMS": ("databases", "custom_hmms"),
    "ASM_HMMS": ("databases", "asm_hmms"),
    "EFFECTOR_HMMS": ("databases", "effector_hmms"),
    "SENSOR_HMMS": ("databases", "sensor_hmms"),
}
METADATA_COLUMNS = ("ORGANISM_KINGDOM", "ORGANISM_KARYOTE", "PLOIDY")
OTHER_COLUMNS = ("SAMPLE_ID", "SPECIES_ID", "ASSEMBLY_PROFILE", "NLR_PROFILE", "DISCOVERY_MODE", "ASSEMBLY_VERSION")
COLUMNS = (*OTHER_COLUMNS[:2], *PATH_COLUMNS, "MASKED_GENOME_PATH", *OTHER_COLUMNS[2:], *METADATA_COLUMNS)
COUNT_COLUMNS = ("n_candidates", "n_report_rows", "n_strict_candidates", "n_fusion_models")
SUMMARY_COLUMNS = ("sample_id", "species_id", "output_dir", "status", *METADATA_COLUMNS, *COUNT_COLUMNS, "error")
ABSENT = {"", "none", "na", "null"}


def add_batch_arguments(parser):
    parser.add_argument("--samples", "--input-tsv", "--input_tsv", required=True, dest="samples")
    parser.add_argument("--config", help="Shared YAML/JSON settings; sample cells override nonempty values")
    parser.add_argument("-o", "--outdir", "--output-dir", "--output_dir", required=True, dest="outdir")
    parser.add_argument("-t", "--threads", "--cpu-threads", "--cpu_threads", type=int, dest="threads", help="Threads for the currently active sample; samples run sequentially")
    parser.add_argument("-r", "--ram-gb", "--ram_gb", type=int, dest="ram_gb")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", "--dry_run", action="store_true", dest="dry_run")
    parser.add_argument("--validate-only", action="store_true", dest="validate_only")
    failures = parser.add_mutually_exclusive_group()
    failures.add_argument("--keep-going", action="store_true", dest="keep_going")
    failures.add_argument("--stop-on-error", "--stop_on_error", action="store_false", dest="keep_going",
                          help="Stop after the first sample failure (the default)")
    parser.set_defaults(keep_going=False)


def sample_template():
    return "\t".join(COLUMNS) + "\n"


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _read_sheet(path):
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("Sample sheet must be UTF-8 TSV") from exc
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True)
    try:
        header = next(reader, [])
        if reader.line_num != 1 or not header:
            raise ValueError("Sample sheet needs one single-line header")
        if any(any(ord(char) < 32 or 127 <= ord(char) < 160 for char in cell) for cell in header):
            raise ValueError("Control characters are unsupported in sample-sheet headers")
        header = [cell.strip().upper() for cell in header]
        if len(set(header)) != len(header):
            raise ValueError("Duplicate sample-sheet header")
        unknown = sorted(set(header) - set(COLUMNS))
        if unknown:
            raise ValueError("Unsupported sample-sheet columns: " + ", ".join(unknown))
        missing = {"SAMPLE_ID", "SPECIES_ID"} - set(header)
        if missing:
            raise ValueError("Missing sample-sheet columns: " + ", ".join(sorted(missing)))
        rows, seen = [], set()
        previous_line = reader.line_num
        for cells in reader:
            line = previous_line + 1
            if reader.line_num != line:
                raise ValueError(f"Multiline sample-sheet cells are unsupported at line {line}")
            previous_line = reader.line_num
            if len(cells) != len(header):
                raise ValueError(f"Sample-sheet line {line} has {len(cells)} cells; expected {len(header)}")
            if any(any(ord(char) < 32 or 127 <= ord(char) < 160 for char in cell) for cell in cells):
                raise ValueError(f"Control characters are unsupported in sample-sheet cells at line {line}")
            row = {key: (None if cell.strip().casefold() in ABSENT else cell.strip()) for key, cell in zip(header, cells)}
            sample = row.get("SAMPLE_ID")
            if not sample or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or len(sample) > 128:
                raise ValueError(f"Unsafe SAMPLE_ID at line {line}; use 1–128 letters/digits, dots, underscores or hyphens, starting with a letter/digit")
            if sample.casefold() in seen:
                raise ValueError(f"Duplicate or case-colliding SAMPLE_ID: {sample}")
            seen.add(sample.casefold())
            if not row.get("SPECIES_ID"):
                raise ValueError(f"SPECIES_ID is required for sample {sample}")
            for key in ("ASSEMBLY_PROFILE", "NLR_PROFILE", "DISCOVERY_MODE"):
                if row.get(key):
                    row[key] = row[key].upper()
            if row.get("ASSEMBLY_PROFILE") and row.get("NLR_PROFILE") and row["ASSEMBLY_PROFILE"] != row["NLR_PROFILE"]:
                raise ValueError(f"Sample {sample}: conflicting ASSEMBLY_PROFILE and NLR_PROFILE")
            for key, supported in (("ORGANISM_KINGDOM", "funga"), ("ORGANISM_KARYOTE", "eukaryote")):
                if row.get(key) and row[key].casefold() != supported:
                    raise ValueError(f"Sample {sample}: {key} must be {supported} when provided; this workflow supports annotated fungi")
            rows.append(row)
    except csv.Error as exc:
        raise ValueError(f"Invalid sample-sheet TSV: {exc}") from exc
    if not rows:
        raise ValueError("Sample sheet contains no samples")
    return raw, rows


def _sample_args(config, resume=False):
    data = config.to_dict()
    values = {key: data[section][key] for key in INPUT_FLAGS for section in ("inputs", "databases") if key in data[section]}
    values.update(outdir=data["execution"]["output_dir"], threads=data["execution"]["cpu_threads"],
                  ram_gb=data["execution"]["ram_gb"], sample=data["sample"]["sample_id"],
                  species_id=data["sample"]["species_id"], assembly_version=data["sample"]["assembly_version"],
                  resume=resume, dry_run=False, tools=data["tools"], scientific=config.scientific_settings())
    return Namespace(**values)


def _read_json(path):
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read batch provenance: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Invalid batch JSON object: {path}")
    return value


def _owned_paths(root):
    if root.is_symlink():
        raise ValueError(f"Batch output root cannot be a symlink: {root}")
    if root.exists() and not root.is_dir():
        raise ValueError(f"Batch output root is not a directory: {root}")
    if root.exists():
        for name in ("configs", "samples"):
            if (root / name).exists() and not (root / name).is_dir():
                raise ValueError(f"Batch-owned directory is not a directory: {root / name}")
        for name in ("batch_manifest.json", "batch_state.json", "batch_summary.json", "batch_summary.tsv", "samples.tsv"):
            if (root / name).exists() and not (root / name).is_file():
                raise ValueError(f"Batch-owned file is not a regular file: {root / name}")
        for path in root.rglob("*"):
            if path.is_symlink():
                raise ValueError(f"Symlinks in batch-owned output paths are unsupported: {path}")
            if not path.is_file() and not path.is_dir():
                raise ValueError(f"Unsupported special file in batch output: {path}")
            if not path.resolve().is_relative_to(root.resolve()):
                raise ValueError(f"Batch output path escapes its root: {path}")


def _identity(sheet_sha, samples):
    payload = {"schema_version": 1, "funlr_version": __version__, "source_sheet_sha256": sheet_sha,
               "samples": [{key: sample[key] for key in ("sample_id", "species_id", "output_dir", "metadata", "resolved_config")} for sample in samples]}
    return _digest(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode())


def _prepare(args):
    sheet = Path(args.samples).expanduser().resolve()
    raw, rows = _read_sheet(sheet)
    requested_root = Path(args.outdir).expanduser().absolute()
    if requested_root.is_symlink():
        raise ValueError("Batch output root cannot be a symlink")
    root = requested_root.resolve()
    if root.exists() and not root.is_dir():
        raise ValueError("Batch output root must be a directory")
    if root.exists() and any(root.iterdir()) and not args.resume:
        raise ValueError("Batch output directory is not empty; use --resume for the identical batch or a new directory")
    shared_path = Path(args.config).expanduser().resolve() if args.config else None
    shared = FunlrConfig.from_yaml(shared_path) if shared_path else FunlrConfig()
    samples, prepared = [], []
    all_inputs = {sheet, *([shared_path] if shared_path else [])}
    for row in rows:
        sample = row["SAMPLE_ID"]
        config = copy.deepcopy(shared)  # Retain explicit-setting provenance for profile defaults.
        overrides = {"sample": {"sample_id": sample, "species_id": row["SPECIES_ID"], "assembly_version": row.get("ASSEMBLY_VERSION")},
                     "context": {"profile": row.get("ASSEMBLY_PROFILE") or row.get("NLR_PROFILE")},
                     "discovery": {"discovery_mode": row.get("DISCOVERY_MODE")},
                     "execution": {"output_dir": str(root / "samples" / sample), "cpu_threads": args.threads, "ram_gb": args.ram_gb, "resume": False, "dry_run": False}}
        for column, (section, key) in PATH_COLUMNS.items():
            if row.get(column):
                overrides.setdefault(section, {})[key] = str((sheet.parent / Path(row[column]).expanduser()).resolve())
        row_warnings = []
        if row.get("MASKED_GENOME_PATH"):
            if row.get("GENOME_PATH"):
                row_warnings.append("Both GENOME_PATH and MASKED_GENOME_PATH are set; using GENOME_PATH.")
            else:
                overrides.setdefault("inputs", {})["genome"] = str((sheet.parent / Path(row["MASKED_GENOME_PATH"]).expanduser()).resolve())
        try:
            config.apply_overrides(overrides)
            single_args = _sample_args(config, resume=args.resume)
            inputs, warnings, counts = validate_inputs(single_args)
            warnings = row_warnings + warnings
        except (ValueError, OSError) as exc:
            raise ValueError(f"Sample {sample}: {exc}") from exc
        all_inputs.update(inputs.values())
        entry = {"sample_id": sample, "species_id": row["SPECIES_ID"], "output_dir": f"samples/{sample}",
                 "metadata": {key: row[key] for key in METADATA_COLUMNS if row.get(key)},
                 "config_path": f"configs/{sample}.yaml", "resolved_config": config.to_dict(),
                 "input_counts": counts, "warnings": warnings}
        samples.append(entry)
        prepared.append((single_args, inputs, config.scientific_settings(), warnings, counts))
    for path in all_inputs:
        if path.is_relative_to(root) or root.is_relative_to(path):
            raise ValueError(f"Batch inputs and output root must not contain each other: {path}")
    identity = _identity(_digest(raw), samples)
    manifest = {"schema_version": 1, "funlr_version": __version__, "created": now(), "identity_sha256": identity,
                "source_sheet": {"file": sheet.name, "sha256": _digest(raw)},
                "source_config": {"file": shared_path.name, "sha256": _digest(shared_path.read_bytes())} if shared_path else None,
                "samples": samples}
    if args.resume:
        manifest_path = root / "batch_manifest.json"
        if manifest_path.is_symlink():
            raise ValueError("Batch manifest cannot be a symlink")
        previous = _read_json(manifest_path)
        if previous.get("schema_version") != 1 or previous.get("funlr_version") != __version__ or previous.get("identity_sha256") != identity:
            raise ValueError("Cannot resume: batch sample sheet, version, or effective configurations changed")
        if previous.get("samples") != samples or previous.get("source_sheet", {}).get("sha256") != _digest(raw):
            raise ValueError("Cannot resume: frozen batch provenance differs from the requested samples")
        manifest = previous
    _owned_paths(root)
    if (root / ".batch.lock").exists():
        raise ValueError("Batch lock exists; verify the recorded process has stopped before removing a stale lock")
    if args.resume:
        hashed_inputs = {}
        for sample, parameters in zip(samples, prepared):
            child_manifest = root / sample["output_dir"] / "manifest.json"
            if not child_manifest.exists():
                continue
            previous_inputs = _read_json(child_manifest).get("inputs")
            if not isinstance(previous_inputs, dict) or set(previous_inputs) != set(parameters[1]):
                raise ValueError(f"Cannot resume: child input provenance is invalid for {sample['sample_id']}")
            for key, input_path in parameters[1].items():
                expected = previous_inputs[key]
                if not isinstance(expected, dict) or expected.get("path") != str(input_path):
                    raise ValueError(f"Cannot resume: child input path changed for {sample['sample_id']}: {key}")
                if input_path not in hashed_inputs:
                    digest = hashlib.sha256()
                    with input_path.open("rb") as handle:
                        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                            digest.update(chunk)
                    hashed_inputs[input_path] = digest.hexdigest()
                if expected.get("sha256") != hashed_inputs[input_path]:
                    raise ValueError(f"Cannot resume: input contents changed for {sample['sample_id']}: {key}")
        saved_sheet = root / "samples.tsv"
        if saved_sheet.exists() and (not saved_sheet.is_file() or saved_sheet.read_bytes() != raw):
            raise ValueError("Cannot resume: saved samples.tsv differs from the source sheet")
        for sample in samples:
            path = root / sample["config_path"]
            if path.exists():
                try:
                    saved = yaml.safe_load(path.read_text())
                except (OSError, yaml.YAMLError) as exc:
                    raise ValueError(f"Cannot read saved sample configuration: {path}") from exc
                if saved != sample["resolved_config"]:
                    raise ValueError(f"Cannot resume: saved configuration changed for {sample['sample_id']}")
    return root, raw, manifest, prepared


def _atomic_bytes(path, data):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".batch-write-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(data)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _write_json(path, data):
    _atomic_bytes(path, (json.dumps(data, indent=2, sort_keys=True) + "\n").encode())


def _counts(root):
    final = root / "final_results"
    records = []
    for name in ("nlr_final_report.tsv", "nlr_strict_candidates.tsv"):
        with (final / name).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"protein_id", "tier"}.issubset(reader.fieldnames or []):
                raise RuntimeError(f"Invalid completed sample report: {final / name}")
            records.append(list(reader))
    fusions = sum(str(row.get("is_fusion_model", "")).lower() in {"true", "1", "yes"} for row in records[0])
    return dict(zip(COUNT_COLUMNS, (len(records[0]) - fusions, len(records[0]), len(records[1]), fusions)))


def _publish_progress(root, manifest, state):
    state["updated"] = now()
    _write_json(root / "batch_state.json", state)
    rows = []
    for sample, progress in zip(manifest["samples"], state["samples"]):
        rows.append({"sample_id": sample["sample_id"], "species_id": sample["species_id"], "output_dir": sample["output_dir"],
                     "status": progress["status"], **{key: sample["metadata"].get(key, "") for key in METADATA_COLUMNS},
                     **{key: progress.get("counts", {}).get(key) for key in COUNT_COLUMNS}, "error": progress.get("error", "")})
    summary = {"schema_version": 1, "status": state["status"], "sample_status_counts": dict(Counter(row["status"] for row in rows)), "samples": rows}
    _write_json(root / "batch_summary.json", summary)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    _atomic_bytes(root / "batch_summary.tsv", buffer.getvalue().encode())
    return summary


def run_batch(args):
    """Validate every row before writes/tools, then run each sample sequentially."""
    root, raw, manifest, prepared = _prepare(args)
    if args.dry_run or args.validate_only:
        return {"status": "PLANNED" if args.dry_run else "VALIDATED", "outdir": str(root),
                "sample_count": len(manifest["samples"]), "samples": manifest["samples"]}
    root.mkdir(parents=True, exist_ok=True)
    lock = root / ".batch.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise ValueError("Another process owns the batch lock") from None
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(json.dumps({"pid": os.getpid(), "started": now()}) + "\n")
        _owned_paths(root)
        for name in ("configs", "samples"):
            (root / name).mkdir(exist_ok=True)
        if not args.resume:
            _write_json(root / "batch_manifest.json", manifest)
        if not (root / "samples.tsv").exists():
            _atomic_bytes(root / "samples.tsv", raw)
        for sample in manifest["samples"]:
            path = root / sample["config_path"]
            if not path.exists():
                _atomic_bytes(path, yaml.safe_dump(sample["resolved_config"], sort_keys=False).encode())
        state = {"schema_version": 1, "identity_sha256": manifest["identity_sha256"], "status": "RUNNING", "started": now(),
                 "samples": [{"sample_id": sample["sample_id"], "species_id": sample["species_id"], "output_dir": sample["output_dir"], "status": "NOT_RUN"} for sample in manifest["samples"]]}
        _publish_progress(root, manifest, state)
        for sample, progress, parameters in zip(manifest["samples"], state["samples"], prepared):
            progress.update(status="RUNNING", started=now())
            _publish_progress(root, manifest, state)
            try:
                run_sample(*parameters)
                progress.update(status="COMPLETE", counts=_counts(root / sample["output_dir"]))
            except KeyboardInterrupt:
                progress.update(status="INTERRUPTED", finished=now(), error="Interrupted")
                state["status"] = "INTERRUPTED"
                _publish_progress(root, manifest, state)
                raise
            except Exception as exc:
                progress.update(status="FAILED", error=f"{type(exc).__name__}: {exc}")
                if not args.keep_going:
                    progress["finished"] = now()
                    break
            progress["finished"] = now()
            _publish_progress(root, manifest, state)
        state["status"] = "FAILED" if any(row["status"] != "COMPLETE" for row in state["samples"]) else "COMPLETE"
        state["finished"] = now()
        return {**_publish_progress(root, manifest, state), "outdir": str(root)}
    finally:
        lock.unlink(missing_ok=True)
