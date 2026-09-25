"""Read recorded run results without running tools or modifying pipeline files.

The inspector reports recorded state, not process liveness or biological
validation. Verification checks recorded output bytes against recorded hashes;
it does not authenticate the manifest or revalidate the original analysis.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime
import hashlib
from html import escape
import io
import json
import os
from pathlib import Path, PurePosixPath
import re

from . import __version__


SCHEMA_VERSION = 1
STAGE_NAMES = (
    "Validate inputs", "Discover candidates", "Domain architecture",
    "Tier and flag", "Export FASTAs", "Miniprot rescue",
    "Exonerate rescue", "Final reports",
)
REVIEW_LIMIT = 200
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_FINAL_REPORT = "final_results/nlr_final_report.tsv"
_STRICT_REPORT = "final_results/nlr_strict_candidates.tsv"
_PROTECTED = {"results", "final_results", "work", "logs", "state", "configs", "samples"}


def _hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def _hash_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_path(root, relative):
    """Accept only a relative path whose existing components are not symlinks."""
    if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
        raise ValueError("Expected a nonempty relative POSIX path")
    path = PurePosixPath(relative)
    if path.is_absolute() or ".." in path.parts or str(path) in {"", "."}:
        raise ValueError("Path must remain within the run directory")
    current = root
    for component in path.parts:
        current /= component
        if current.is_symlink():
            raise ValueError("Symlinked run paths are not supported")
    if not current.resolve(strict=False).is_relative_to(root):
        raise ValueError("Path escapes the run directory")
    return current


def _read_json(root, relative, errors, snapshots):
    try:
        path = _safe_path(root, relative)
        if not path.exists():
            return None
        data = path.read_bytes()
        snapshots[relative] = _hash_bytes(data)
        value = json.loads(data)
        if not isinstance(value, dict):
            raise ValueError("Expected a JSON object")
        return value
    except (OSError, ValueError, UnicodeError) as exc:
        errors.append(f"Cannot read {relative}: {exc}")
        return None


def _lock_info(root, name, errors, snapshots):
    path = root / name
    present = path.exists() or path.is_symlink()
    metadata = _read_json(root, name, errors, snapshots) if present else None
    return {"path": name, "present": present, "metadata": metadata,
            "process_liveness": "NOT_CHECKED"}


def _new_integrity(verify):
    return {
        "status": "FAILED" if verify else "NOT_CHECKED",
        "scope": "Recorded stage and final output files; original inputs, tools, and manifest authenticity are not revalidated.",
        "expected_files": 0, "checked_files": 0, "matched_files": 0,
        "changed_count": 0, "missing_count": 0, "unsafe_count": 0,
        "unreadable_count": 0, "invalid_record_count": 0, "issues": [],
    }


def _empty_counts():
    return {"n_report_rows": None, "n_candidates": None,
            "n_unique_candidate_ids": None, "n_unique_report_ids": None,
            "n_strict_candidates": None, "n_fusion_models": None,
            "n_fusion_review_required": None, "n_flagged_candidates": None,
            "n_review_candidates": None}


def _read_table(root, relative, required, errors, observed_hashes, missing_is_error):
    try:
        path = _safe_path(root, relative)
        if not path.exists():
            if missing_is_error:
                errors.append(f"Missing final table: {relative}")
            return None
        raw = path.read_bytes()
        observed_hashes[relative] = _hash_bytes(raw)
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("Missing required table columns: " + ", ".join(sorted(required)))
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Duplicate table column names")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("A table row has a different number of fields from its header")
            rows.append(row)
        return rows
    except (OSError, ValueError, UnicodeError, csv.Error) as exc:
        errors.append(f"Cannot read {relative}: {exc}")
        return None


def _table_summary(root, report, complete, observed_hashes):
    rows = _read_table(root, _FINAL_REPORT, {"protein_id", "tier", "flags"},
                       report["errors"], observed_hashes, complete)
    strict = _read_table(root, _STRICT_REPORT, {"protein_id"},
                         report["errors"], observed_hashes, complete)
    if strict is not None:
        report["counts"]["n_strict_candidates"] = len(strict)
    if rows is None:
        return
    tiers, flags = Counter(), Counter()
    fusion_count = fusion_review = flagged_count = review_count = 0
    ids, candidate_ids = set(), set()
    review_rows = []
    for row in rows:
        ids.add(row["protein_id"])
        tiers[row["tier"]] += 1
        row_flags = {part.strip() for part in row["flags"].split(";")
                     if part.strip() and part.strip() not in {"NA", "nan"}}
        flags.update(row_flags)
        is_flagged = bool(row_flags - {"PASS"})
        is_fusion = (row.get("is_fusion_model", "").strip().lower() in {"true", "yes", "1", "t"}
                     or row["tier"].startswith("FUSION"))
        if not is_fusion:
            candidate_ids.add(row["protein_id"])
        needs_review = "REVIEW_REQUIRED" in row_flags
        fusion_count += is_fusion
        fusion_review += is_fusion and needs_review
        flagged_count += is_flagged
        if is_flagged or is_fusion or "NEEDS_REVIEW" in row["tier"] or row["tier"].startswith("TIER_2"):
            review_count += 1
            if len(review_rows) < REVIEW_LIMIT:
                review_rows.append({key: row.get(key, "") for key in (
                    "protein_id", "tier", "rescue_priority", "flags", "nbd_confidence",
                    "fusion_members", "fusion_evidence_rule",
                )})
    report["counts"].update(
        n_report_rows=len(rows), n_candidates=len(rows) - fusion_count,
        n_unique_candidate_ids=len(candidate_ids), n_unique_report_ids=len(ids),
        n_fusion_models=fusion_count, n_fusion_review_required=fusion_review,
        n_flagged_candidates=flagged_count, n_review_candidates=review_count,
    )
    report["tier_counts"] = dict(sorted(tiers.items()))
    report["flag_counts"] = dict(sorted(flags.items()))
    report["review_rows"] = review_rows
    report["review_selection"] = {
        "rule": "Rows with non-PASS flags, fusion hypotheses, NEEDS_REVIEW tiers, or Tier 2 rescue tiers.",
        "limit": REVIEW_LIMIT, "shown": len(review_rows), "total": review_count,
        "order": "Original final-report row order",
    }
    if len(ids) != len(rows):
        report["warnings"].append("The final report contains repeated protein IDs; counts report rows as well as unique IDs.")


def _records(manifest, state, integrity):
    records = {}
    sources = [("manifest.final_outputs", manifest.get("final_outputs"))]
    stages = state.get("stages", {})
    if isinstance(stages, dict):
        for number, stage in sorted(stages.items(), key=lambda item: str(item[0])):
            if isinstance(stage, dict) and stage.get("status") == "complete":
                sources.append((f"state.stages.{number}.outputs", stage.get("outputs")))
    for label, values in sources:
        if not isinstance(values, dict) or not values:
            integrity["issues"].append({"kind": "invalid_record", "path": label,
                                         "message": "Missing or empty output checksum mapping"})
            continue
        for name, digest in values.items():
            if not isinstance(name, str) or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
                integrity["issues"].append({"kind": "invalid_record", "path": str(name),
                                             "message": f"Invalid output SHA256 in {label}"})
                continue
            digest = digest.lower()
            if name in records and records[name] != digest:
                integrity["issues"].append({"kind": "invalid_record", "path": name,
                                             "message": "Conflicting recorded output checksums"})
                continue
            records[name] = digest
    return records


def _verify(root, manifest, state, integrity, observed_hashes):
    records = _records(manifest, state, integrity)
    integrity["expected_files"] = len(records)
    for relative, expected in sorted(records.items()):
        try:
            path = _safe_path(root, relative)
        except ValueError as exc:
            integrity["issues"].append({"kind": "unsafe", "path": relative, "message": str(exc)})
            continue
        if not path.exists():
            integrity["issues"].append({"kind": "missing", "path": relative})
            continue
        if not path.is_file():
            integrity["issues"].append({"kind": "unreadable", "path": relative,
                                         "message": "Recorded output is not a regular file"})
            continue
        try:
            actual = observed_hashes.get(relative) or _hash_file(path)
            integrity["checked_files"] += 1
        except OSError as exc:
            integrity["issues"].append({"kind": "unreadable", "path": relative, "message": str(exc)})
            continue
        if actual != expected:
            integrity["issues"].append({"kind": "changed", "path": relative,
                                         "expected_sha256": expected, "observed_sha256": actual})
        else:
            integrity["matched_files"] += 1
    final_outputs = manifest.get("final_outputs")
    final_outputs = final_outputs if isinstance(final_outputs, dict) else {}
    for relative in (_FINAL_REPORT, _STRICT_REPORT):
        if relative not in final_outputs:
            integrity["issues"].append({"kind": "invalid_record", "path": relative,
                                         "message": "Required final table lacks a manifest checksum"})
    for kind in ("changed", "missing", "unsafe", "unreadable", "invalid_record"):
        integrity[kind + "_count"] = sum(item["kind"] == kind for item in integrity["issues"])


def _check_snapshot(root, snapshots, errors):
    for relative, expected in snapshots.items():
        try:
            if _hash_file(_safe_path(root, relative)) != expected:
                errors.append(f"Metadata changed during inspection: {relative}")
        except (OSError, ValueError):
            errors.append(f"Metadata became unavailable during inspection: {relative}")


def _elapsed(started, finished):
    """Elapsed recorded interval, never a CPU-time or process-liveness estimate."""
    try:
        seconds = (datetime.fromisoformat(finished.replace("Z", "+00:00"))
                   - datetime.fromisoformat(started.replace("Z", "+00:00"))).total_seconds()
        return seconds if seconds >= 0 else None
    except (AttributeError, TypeError, ValueError):
        return None


def _supplemental_summary(root, report, manifest, state, snapshots, observed_hashes):
    records = {}
    stages = state.get("stages", {})
    for stage in stages.values() if isinstance(stages, dict) else []:
        if isinstance(stage, dict) and isinstance(stage.get("outputs"), dict):
            records.update(stage["outputs"])
    if isinstance(manifest.get("final_outputs"), dict):
        records.update(manifest["final_outputs"])
    for key, relative in (("input_qc", "results/stage0/input_qc.json"),
                           ("evidence_summary", "final_results/evidence/evidence_summary.json")):
        report[key] = None
        # A stray JSON added after a run cannot become apparently verified QC.
        if relative not in records:
            continue
        value = _read_json(root, relative, report["errors"], snapshots)
        if value is not None:
            if value.get("schema_version") != 1:
                report["warnings"].append(f"Unsupported summary schema: {relative}")
                continue
            report[key] = value
            observed_hashes[relative] = snapshots[relative]
    byte_count = present = unavailable = 0
    for relative in records:
        try:
            path = _safe_path(root, relative)
            if not path.is_file():
                raise ValueError("Not a regular file")
            byte_count += path.stat().st_size
            present += 1
        except (OSError, ValueError, TypeError):
            unavailable += 1
    report["resource_observations"] = {
        "latest_invocation_elapsed_seconds": _elapsed(manifest.get("started"), manifest.get("finished")),
        "recorded_output_paths": len(records), "observed_output_paths": present,
        "unavailable_output_paths": unavailable, "observed_recorded_output_bytes": byte_count,
        "peak_memory_bytes": None, "cpu_seconds": None,
        "scope": "Elapsed time uses saved timestamps, not CPU time. Output sizes are current logical bytes of unique recorded paths, including mirrored products; unrecorded work files and original inputs are excluded. Peak memory was not measured. A resumed invocation is not the full analysis duration.",
    }


def _single_run(root, verify):
    report = {
        "schema_version": SCHEMA_VERSION, "reporter_version": __version__, "kind": "run", "run_dir": str(root),
        "status": "UNKNOWN", "status_detail": "No recognized run metadata.",
        "recorded_status": None, "errors": [], "warnings": [],
        "integrity": _new_integrity(verify), "counts": _empty_counts(),
        "tier_counts": {}, "flag_counts": {}, "review_rows": [], "review_selection": None,
        "provenance": {}, "stage_progress": {"complete": 0, "total": 8, "stages": []},
    }
    snapshots = {}
    manifest = _read_json(root, "manifest.json", report["errors"], snapshots)
    state = _read_json(root, "state.json", report["errors"], snapshots)
    resolved = _read_json(root, "resolved_config.json", report["errors"], snapshots)
    lock = _lock_info(root, ".run.lock", report["errors"], snapshots)
    report["lock"] = lock
    manifest, state = manifest or {}, state or {}
    recorded = str(manifest.get("status", "")).lower()
    report["recorded_status"] = recorded or None
    source_warnings = manifest.get("warnings", [])
    if isinstance(source_warnings, list):
        report["warnings"].extend(str(item) for item in source_warnings)
    elif source_warnings:
        report["warnings"].append(str(source_warnings))
    if manifest.get("error"):
        report["errors"].append("Recorded run error: " + str(manifest["error"]))
    report["provenance"] = {key: manifest.get(key) for key in (
        "version", "sample", "species_id", "assembly_version", "signature", "baseline",
        "started", "finished", "threads", "ram_gb_label", "runtime", "inputs", "tools",
        "code", "input_counts", "upstream_commit",
    )}
    report["provenance"]["resolved_settings"] = resolved if resolved is not None else manifest.get("config", {})
    if resolved is not None and manifest.get("config") is not None and resolved != manifest["config"]:
        report["errors"].append("resolved_config.json differs from the manifest configuration.")
    stages = state.get("stages", {})
    if not isinstance(stages, dict):
        report["errors"].append("state.json has an invalid stages mapping.")
        stages = {}
    for number, name in enumerate(STAGE_NAMES):
        data = stages.get(str(number), {})
        if not isinstance(data, dict):
            report["errors"].append(f"Invalid stage {number} record.")
            data = {}
        status = data.get("status", "not_recorded")
        outputs = data.get("outputs", {})
        valid_complete = status == "complete" and isinstance(outputs, dict) and bool(outputs)
        report["stage_progress"]["complete"] += valid_complete
        report["stage_progress"]["stages"].append({
            "number": number, "name": name, "recorded_status": status,
            "started": data.get("started"), "finished": data.get("finished"),
            "elapsed_seconds": _elapsed(data.get("started"), data.get("finished")),
            "recorded_output_count": len(outputs) if isinstance(outputs, dict) else 0,
        })
    signatures_match = (isinstance(manifest.get("signature"), str)
                        and bool(manifest["signature"])
                        and manifest["signature"] == state.get("signature"))
    outputs_recorded = isinstance(manifest.get("final_outputs"), dict) and bool(manifest["final_outputs"])
    if recorded == "failed":
        report.update(status="FAILED", status_detail="The run manifest records a failed run.")
    elif recorded == "running":
        if lock["present"]:
            report.update(status="RUNNING", status_detail="Running metadata and a run lock are present; process liveness was not checked.")
        else:
            report.update(status="INCOMPLETE", status_detail="Running metadata has no run lock; the run may have been interrupted.")
    elif recorded == "complete":
        if (report["stage_progress"]["complete"] == 8 and signatures_match
                and outputs_recorded and not lock["present"] and not report["errors"]):
            report.update(status="RECORDED_COMPLETE", status_detail="Completion is recorded; output integrity has not been checked.")
        else:
            report.update(status="INCOMPLETE", status_detail="The completion record is inconsistent with stage state, signature, lock, or final-output metadata.")
    elif manifest or state or lock["present"]:
        report.update(status="INCOMPLETE", status_detail="Run files are present, but no complete run is recorded.")
    observed_hashes = {}
    _table_summary(root, report, recorded == "complete", observed_hashes)
    _supplemental_summary(root, report, manifest, state, snapshots, observed_hashes)
    if report["errors"] and report["status"] == "RECORDED_COMPLETE":
        report.update(status="INCOMPLETE", status_detail="Recorded completion has missing or unreadable metadata or final tables.")
    if verify:
        _verify(root, manifest, state, report["integrity"], observed_hashes)
    _check_snapshot(root, snapshots, report["errors"])
    if (root / ".run.lock").exists() != lock["present"]:
        report["errors"].append("The run lock changed during inspection.")
    if report["errors"] and report["status"] == "RECORDED_COMPLETE":
        report.update(status="INCOMPLETE", status_detail="Run metadata changed or could not be read consistently.")
    if verify:
        integrity = report["integrity"]
        if report["status"] == "RECORDED_COMPLETE" and not integrity["issues"] and not report["errors"]:
            integrity["status"] = "VERIFIED"
            report["status_detail"] = "Recorded completion and all recorded stage/final output checksums were verified."
        else:
            integrity["status"] = "FAILED"
            if integrity["issues"] and recorded == "complete":
                report.update(status="FAILED", status_detail="The run was recorded complete, but output verification failed.")
    return report


def _batch_run(root, verify):
    result = {"schema_version": SCHEMA_VERSION, "reporter_version": __version__, "kind": "batch", "run_dir": str(root),
              "status": "UNKNOWN", "status_detail": "No valid batch manifest.",
              "errors": [], "warnings": [], "integrity": _new_integrity(verify), "samples": []}
    snapshots = {}
    manifest = _read_json(root, "batch_manifest.json", result["errors"], snapshots) or {}
    state = _read_json(root, "batch_state.json", result["errors"], snapshots) or {}
    result["lock"] = _lock_info(root, ".batch.lock", result["errors"], snapshots)
    result["recorded_status"] = state.get("status")
    result["provenance"] = {key: manifest.get(key) for key in (
        "funlr_version", "created", "source_sheet", "identity_sha256",
    )}
    samples = manifest.get("samples")
    if manifest.get("schema_version") != 1 or not isinstance(samples, list) or not samples:
        result["errors"].append("Missing, empty, or unsupported batch manifest schema.")
        return result
    identity = manifest.get("identity_sha256")
    if not isinstance(identity, str) or not _SHA256.fullmatch(identity):
        result["errors"].append("Invalid batch manifest identity SHA256.")
    if state.get("schema_version") != 1:
        result["errors"].append("Missing or unsupported batch state schema.")
    if state.get("identity_sha256") != identity:
        result["errors"].append("Batch state identity does not match the manifest.")
    states = state.get("samples", [])
    by_id, state_casefold_ids = {}, set()
    if not isinstance(states, list):
        result["errors"].append("Invalid batch state sample list.")
        states = []
    for row in states:
        if not isinstance(row, dict) or not isinstance(row.get("sample_id"), str):
            result["errors"].append("Invalid batch state sample record.")
            continue
        sid = row["sample_id"]
        if sid.casefold() in state_casefold_ids:
            result["errors"].append("Batch state has repeated or case-conflicting sample IDs.")
        state_casefold_ids.add(sid.casefold())
        by_id[sid] = row
    seen = set()
    for sample in samples:
        try:
            if not isinstance(sample, dict) or not isinstance(sample.get("sample_id"), str):
                raise ValueError("Invalid batch sample record")
            sid = sample["sample_id"]
            if sid.casefold() in seen:
                raise ValueError("Repeated or case-conflicting batch sample ID")
            seen.add(sid.casefold())
            relative = sample.get("output_dir")
            child = _safe_path(root, relative)
            if PurePosixPath(relative).parts != ("samples", sid):
                raise ValueError("Batch sample output must be samples/SAMPLE_ID")
            if sid in by_id and by_id[sid].get("output_dir") != relative:
                result["errors"].append(f"Batch state output directory differs for sample {sid}.")
            # A child is always a single run. Never recurse into untrusted
            # nested batch manifests or follow arbitrary manifest links.
            inspected = _single_run(child, verify)
            result["samples"].append({
                "sample_id": sid, "species_id": sample.get("species_id"),
                "metadata": sample.get("metadata", {}), "output_dir": relative,
                "recorded_batch_status": by_id.get(sid, {}).get("status"),
                **{key: inspected[key] for key in ("status", "status_detail", "counts", "stage_progress", "integrity", "errors", "warnings", "lock")},
            })
        except (OSError, ValueError, TypeError) as exc:
            result["errors"].append(f"Invalid batch sample: {exc}")
    expected_ids = {row["sample_id"] for row in samples
                    if isinstance(row, dict) and isinstance(row.get("sample_id"), str)}
    if set(by_id) != expected_ids:
        result["errors"].append("Batch state sample IDs do not exactly match the manifest.")
    statuses = Counter(sample["status"] for sample in result["samples"])
    result["sample_counts"] = {"total": len(samples), **dict(sorted(statuses.items()))}
    recorded = str(state.get("status", "")).upper()
    if statuses["FAILED"] or recorded == "FAILED":
        result.update(status="FAILED", status_detail="A batch or child run records a failure.")
    elif recorded == "RUNNING" and result["lock"]["present"]:
        result.update(status="RUNNING", status_detail="Batch running metadata and lock are present; process liveness was not checked.")
    elif (recorded == "COMPLETE" and statuses["RECORDED_COMPLETE"] == len(samples)
          and not result["lock"]["present"] and not result["errors"]):
        result.update(status="RECORDED_COMPLETE", status_detail="All child runs have recorded completion.")
    else:
        result.update(status="INCOMPLETE", status_detail="The batch is incomplete; running metadata without a lock may indicate interruption.")
    _check_snapshot(root, snapshots, result["errors"])
    if (root / ".batch.lock").exists() != result["lock"]["present"]:
        result["errors"].append("The batch lock changed during inspection.")
    if result["errors"] and result["status"] == "RECORDED_COMPLETE":
        result.update(status="INCOMPLETE", status_detail="Batch metadata changed or is inconsistent.")
    if verify:
        integrity = result["integrity"]
        integrity["scope"] = "Recorded stage and final outputs of all manifest-listed child runs; batch-summary tables and original inputs are not revalidated."
        for child in result["samples"]:
            for key in ("expected_files", "checked_files", "matched_files", "changed_count", "missing_count", "unsafe_count", "unreadable_count", "invalid_record_count"):
                integrity[key] += child["integrity"][key]
            integrity["issues"].extend({"sample_id": child["sample_id"], **issue}
                                       for issue in child["integrity"]["issues"])
        if (result["status"] == "RECORDED_COMPLETE"
                and all(row["integrity"]["status"] == "VERIFIED" for row in result["samples"])):
            integrity["status"] = "VERIFIED"
            result["status_detail"] = "All recorded child-run outputs were verified."
    return result


def inspect_run(run_dir, verify=False):
    """Return schema-1 recorded status without creating files or running tools.

    With ``verify=True``, hash every unique recorded stage/final output. A lock
    is only an observed file; this function never infers whether its PID lives.
    """
    candidate = Path(run_dir).expanduser().absolute()
    if candidate.is_symlink():
        raise ValueError("The run directory must not be a symlink")
    root = candidate.resolve()
    if root.exists() and not root.is_dir():
        raise ValueError("The run path is not a directory")
    if (root / "batch_manifest.json").exists() or (root / "batch_manifest.json").is_symlink():
        return _batch_run(root, bool(verify))
    return _single_run(root, bool(verify))


def format_status(report):
    """Format an inspection result as concise plain text."""
    lines = [f"FuNLR {report['kind']} status: {report['status']}",
             f"Directory: {report['run_dir']}", report["status_detail"],
             f"Integrity: {report['integrity']['status']}"]
    if report["kind"] == "batch":
        for sample in report["samples"]:
            lines.append(f"  {sample['sample_id']}: {sample['status']} (integrity {sample['integrity']['status']})")
    else:
        progress = report["stage_progress"]
        lines.append(f"Stages recorded complete: {progress['complete']}/{progress['total']}")
        counts = report["counts"]
        if counts["n_candidates"] is not None:
            lines.append(f"Observed final rows: {counts['n_report_rows']}; original candidates: {counts['n_candidates']}; strict candidates: {counts['n_strict_candidates']}; fusion hypotheses: {counts['n_fusion_models']}")
    if report["integrity"]["status"] != "NOT_CHECKED":
        integrity = report["integrity"]
        lines.append(f"Output checks: {integrity['matched_files']}/{integrity['expected_files']} match; {integrity['changed_count']} changed, {integrity['missing_count']} missing, {integrity['unsafe_count']} unsafe, {integrity['unreadable_count']} unreadable, {integrity['invalid_record_count']} invalid records")
    lines.append("Lock present: " + ("yes" if report.get("lock", {}).get("present") else "no") + "; process liveness not checked")
    lines.extend("Warning: " + item for item in report["warnings"])
    lines.extend("Error: " + item for item in report["errors"])
    return "\n".join(lines) + "\n"


def _html_table(headers, rows):
    head = "".join("<th scope=\"col\">" + escape(str(value)) + "</th>" for value in headers)
    body = "".join("<tr>" + "".join("<td>" + escape(str(value)) + "</td>" for value in row) + "</tr>" for row in rows)
    return "<div class=\"table-wrap\"><table><thead><tr>" + head + "</tr></thead><tbody>" + body + "</tbody></table></div>"


def _supplemental_html(report):
    parts = []
    qc = report.get("input_qc")
    if isinstance(qc, dict):
        observations = []
        for section, metrics in (
            ("genome", ("n_sequences", "total_sequence_characters", "n50_bp", "pct_n_bases", "pct_lowercase_acgt")),
            ("proteins", ("n_sequences", "min_length", "max_length", "n_with_internal_stop", "n_with_terminal_stop", "n_starting_with_m")),
            ("gff3", ("n_gene_features", "n_distinct_gene_ids", "n_transcript_features", "n_distinct_transcript_ids")),
        ):
            values = qc.get(section, {})
            if isinstance(values, dict):
                observations.extend((section, metric.replace("_", " "), values.get(metric) if values.get(metric) is not None else "Not recorded") for metric in metrics)
        parts += ["<h2>Input observations</h2><p>Descriptive measurements of the supplied files. These do not grade assembly quality, establish complete open reading frames, or infer missing genes.</p>",
                  _html_table(["Input", "Observation", "Value"], observations)]
    evidence = report.get("evidence_summary")
    if isinstance(evidence, dict):
        parts.append("<h2>Discovery evidence and rescue accounting</h2><p>Source channels overlap; their counts must not be added. Recorded decisions are unchanged. Full records are in <code>final_results/evidence/</code>.</p>")
        channels = evidence.get("discovery_channels", {})
        if isinstance(channels, dict):
            parts.append(_html_table(["Recorded discovery channel", "Unique protein IDs"], [(key, value if value is not None else "Not recorded") for key, value in channels.items()]))
        rescue = evidence.get("rescue", {})
        if isinstance(rescue, dict):
            rows = []
            for track, counts in rescue.items():
                if isinstance(counts, dict):
                    rows.append([track, *[counts.get(key) if counts.get(key) is not None else "Not recorded" for key in (
                        "miniprot_requested", "exonerate_requested", "exonerate_selected", "exonerate_not_selected", "exonerate_features_observed")]])
            parts.append(_html_table(["Track", "Miniprot queries", "Original queries requesting refinement", "Selected originals", "Not selected", "Queries with Exonerate features (includes fusions)"], rows))
        counts = evidence.get("counts", {})
        if isinstance(counts, dict):
            parts.append(_html_table(["Hit evidence", "Count"], [(key.replace("_", " "), counts.get(key, "Not recorded")) for key in (
                "reported_hits", "accepted_reported_hits", "unassessed_reported_hits", "exported_segments")]))
        scan = evidence.get("scan_status", {})
        if isinstance(scan, dict) and scan:
            parts.append(_html_table(["Optional Pfam-NBD scan", "Recorded status"], [(key.replace("_", " "), scan.get(key) if scan.get(key) is not None else "Not recorded") for key in ("pfam_nbd_requested", "pfam_nbd_completed")]))
        parts.append("<p>Not selected is distinct from no alignment. Aligned evidence FASTA records are individual query spans accepted by the recorded discovery filters, with valid coordinates; they may overlap and are not complete NBD sequences or a multiple sequence alignment. Rejected and unassessed raw hits remain in the mapping table without sequences. An empty optional scan does not establish biological absence. Review and exclusion records explain saved membership and flags, not biological function.</p>")
    resources = report.get("resource_observations", {})
    if resources:
        parts.append("<h2>Recorded duration and output size</h2>")
        parts.append(_html_table(["Observation", "Value"], [(key.replace("_", " "), resources.get(key) if resources.get(key) is not None else "Not measured / unavailable") for key in (
            "latest_invocation_elapsed_seconds", "observed_recorded_output_bytes", "observed_output_paths", "unavailable_output_paths", "peak_memory_bytes", "cpu_seconds")]))
        parts.append("<p>" + escape(resources["scope"]) + "</p>")
    return "".join(parts)


def _html(report):
    p, counts = report["provenance"], report["counts"]
    sample = p.get("sample") or "Unnamed sample"
    labels = {"started": "latest invocation started", "finished": "latest invocation finished"}
    details = _html_table(["Recorded field", "Value"], [(labels.get(key, key.replace("_", " ")), p.get(key) or "Not recorded")
                           for key in ("sample", "species_id", "assembly_version", "version", "started", "finished", "threads")])
    metrics = "".join("<div class=\"metric\"><strong>" + str(counts[key]) + "</strong><span>" + label + "</span></div>"
                      for key, label in (("n_candidates", "original candidates"), ("n_strict_candidates", "strict candidates"),
                                         ("n_fusion_models", "fusion hypotheses"), ("n_fusion_review_required", "fusions explicitly flagged for review")))
    tier_table = _html_table(["Tier", "Report rows"], report["tier_counts"].items())
    flag_table = _html_table(["Flag", "Rows carrying flag"], report["flag_counts"].items())
    stages = _html_table(["Stage", "Recorded status", "Recorded outputs", "Stage started", "Stage finished", "Elapsed seconds"],
                         [(f"{row['number']}. {row['name']}", row["recorded_status"], row["recorded_output_count"], row["started"] or "Not recorded", row["finished"] or "Not recorded", row.get("elapsed_seconds") if row.get("elapsed_seconds") is not None else "Not recorded")
                          for row in report["stage_progress"]["stages"]])
    review = _html_table(["Protein ID", "Tier", "Priority", "Flags", "NBD confidence", "Fusion members"],
                         [[row[key] for key in ("protein_id", "tier", "rescue_priority", "flags", "nbd_confidence", "fusion_members")]
                          for row in report["review_rows"]])
    selection = report["review_selection"]
    warning_html = "<ul>" + "".join("<li>" + escape(item) + "</li>" for item in report["warnings"]) + "</ul>" if report["warnings"] else "<p>No warnings were recorded.</p>"
    encoded = escape(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; img-src 'none'; base-uri 'none'; form-action 'none'">
<title>FuNLR report — {escape(str(sample))}</title>
<style>
:root{{color-scheme:light;font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#203132;background:#f3f6f5;line-height:1.5}}
*{{box-sizing:border-box}}body{{margin:0}}main{{max-width:1150px;margin:auto;padding:36px 24px 64px}}header{{border-bottom:3px solid #22776b;padding-bottom:22px;margin-bottom:24px}}h1{{font-size:2.2rem;line-height:1.2;margin:8px 0}}h2{{font-size:1.35rem;margin:30px 0 12px}}h3{{font-size:1.05rem}}.eyebrow{{color:#22776b;font-weight:700;letter-spacing:.07em;text-transform:uppercase;font-size:.8rem}}.badge{{display:inline-block;background:#d9eee6;color:#175749;border-radius:4px;padding:5px 10px;font-weight:650}}.note{{background:#eaf0ef;border-left:4px solid #22776b;padding:12px 16px}}.metrics{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin:22px 0}}.metric{{background:white;border:1px solid #d7e0dd;border-radius:7px;padding:18px}}.metric strong{{display:block;font-size:2rem;color:#185e54}}.metric span{{font-size:.9rem}}.table-wrap{{overflow-x:auto;background:white;border:1px solid #d7e0dd;border-radius:6px}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}th,td{{text-align:left;vertical-align:top;padding:10px 12px;border-bottom:1px solid #e3e9e7;overflow-wrap:anywhere}}th{{background:#e9f1ee;font-weight:650}}tr:last-child td{{border-bottom:0}}.two-col{{display:grid;grid-template-columns:1fr 1fr;gap:24px}}details{{background:white;border:1px solid #d7e0dd;border-radius:6px;margin-top:18px;padding:14px}}summary{{cursor:pointer;font-weight:650}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.78rem;line-height:1.5}}footer{{margin-top:32px;font-size:.85rem;color:#526664}}@media(max-width:760px){{.metrics{{grid-template-columns:1fr 1fr}}.two-col{{grid-template-columns:1fr}}main{{padding:24px 16px}}}}@media print{{body{{background:white}}main{{max-width:none;padding:0}}.table-wrap{{overflow:visible}}details{{break-inside:avoid}}}}
</style></head><body><main>
<header><div class="eyebrow">FuNLR · recorded analysis</div><h1>{escape(str(sample))}</h1>
<p>Fungal NLR candidate discovery and rescue report</p><span class="badge">Recorded outputs verified</span></header>
<p class="note">These are computational candidates and fusion hypotheses. Alignment support and a strict-set label do not establish a biologically validated NLR or gene fusion.</p>
<div class="metrics">{metrics}</div>
<h2>Run record</h2>{details}
<p>The final report contains {counts['n_report_rows']} rows, including the fusion hypotheses shown above. Verification matched {report['integrity']['matched_files']} unique recorded stage and final output files. This checks agreement with the saved checksums, not the authenticity of the manifest or the scientific validity of the analysis. Original inputs and executable files were not rehashed.</p>
<div class="two-col"><section><h2>Tier counts</h2>{tier_table}</section><section><h2>Flags</h2><p>Each flag is counted once per report row; rows can carry several flags.</p>{flag_table}</section></div>
<h2>Candidates requiring inspection</h2><p>{escape(selection['rule'])} Showing {selection['shown']} of {selection['total']} rows, limited to the first {selection['limit']} in original report order. The complete records remain in <code>final_results/nlr_final_report.tsv</code>.</p>{review}
<h2>Stage completion</h2>{stages}
{_supplemental_html(report)}
<h2>Recorded warnings</h2>{warning_html}
<details><summary>Full settings, database/input hashes, tool versions, and inspection record</summary><pre>{encoded}</pre></details>
<footer>Created by FuNLR reporter {escape(report['reporter_version'])} from the recorded FuNLR {escape(str(p.get('version') or 'unknown'))} run. This self-contained report has no scripts, remote fonts, or external assets. It is a read-only summary of the recorded run; FuNLR's discovery, tiering, and rescue outputs were not modified.</footer>
</main></body></html>'''


def write_html_report(run_dir, output):
    """Write one new local HTML report only after verifying a completed run.

    Existing files and symlinks are never overwritten. Inside the run directory,
    only a direct child HTML file is allowed; stage/product directories remain
    untouched. No output directories are created implicitly.
    """
    target = Path(output).expanduser().absolute()
    if target.suffix.lower() != ".html":
        raise ValueError("Report output must have a .html extension")
    if any(part.is_symlink() for part in (target, *target.parents)):
        raise ValueError("Report output and its parent directories must not be symlinks")
    if target.exists():
        raise FileExistsError(f"Report output already exists: {target}")
    if not target.parent.is_dir():
        raise ValueError("Report output parent directory must already exist")
    root = Path(run_dir).expanduser().resolve()
    resolved_target = target.resolve(strict=False)
    if resolved_target.is_relative_to(root) and resolved_target.parent != root:
        relative = resolved_target.relative_to(root)
        if any(part in _PROTECTED for part in relative.parts[:-1]):
            raise ValueError("Cannot write an HTML report into pipeline-owned output directories")
        raise ValueError("Inside a run directory, write the HTML directly beside manifest.json")
    report = inspect_run(run_dir, verify=True)
    if report["kind"] != "run":
        raise ValueError("HTML reports currently support one sample run; inspect a batch with status")
    if report["status"] != "RECORDED_COMPLETE" or report["integrity"]["status"] != "VERIFIED":
        raise ValueError("HTML report requires a recorded-complete run with verified outputs: " + report["status_detail"])
    for item in (report["provenance"].get("inputs") or {}).values():
        if isinstance(item, dict) and item.get("path"):
            try:
                is_input = Path(item["path"]).expanduser().resolve() == resolved_target
            except (OSError, ValueError, TypeError):
                continue
            if is_input:
                raise ValueError("Report output must not be an input path")
    content = _html(report).encode("utf-8")
    # Exclusive creation is the final operation after all validation. A failed
    # write removes only the new report file this call just created.
    descriptor = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
    except BaseException:
        target.unlink(missing_ok=True)
        raise
    return {"schema_version": SCHEMA_VERSION, "reporter_version": __version__, "kind": "html_report", "path": str(resolved_target),
            "run_dir": report["run_dir"], "status": report["status"],
            "integrity": report["integrity"], "counts": report["counts"],
            "sha256": _hash_bytes(content), "size_bytes": len(content)}
