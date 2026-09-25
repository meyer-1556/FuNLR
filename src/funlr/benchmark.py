"""Compare an explicit called-ID set to a curated reference without rerunning it."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

_LABELS = {"positive", "negative", "uncertain"}
_EVIDENCE = {"experimental", "curated", "computational", "synthetic"}
_OVERLAP = {"yes", "no", "unknown"}


def _read(path) -> tuple[Path, bytes, str]:
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"Not a readable file: {path}")
    data = path.read_bytes()
    return path, data, hashlib.sha256(data).hexdigest()


def _text(value, field: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"{field} must be nonempty single-line text")
    return value


def _identifier(value, field: str) -> str:
    value = _text(value, field)
    if any(c.isspace() for c in value):
        raise ValueError(f"{field} must be one FASTA identifier without whitespace")
    return value


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate reference JSON field: {key}")
        result[key] = value
    return result


def _input_ids(data: bytes) -> set[str]:
    """Require unique FASTA first-token IDs and nonempty sequences."""
    result = set()
    active = None
    length = 0
    for line in data.decode("utf-8").splitlines():
        if not line.strip():
            continue
        if line.startswith(">"):
            if active is not None and not length:
                raise ValueError(f"Empty input FASTA sequence: {active}")
            fields = line[1:].split()
            if not fields:
                raise ValueError("Empty input FASTA identifier")
            active = _text(fields[0], "input identifier")
            if active in result:
                raise ValueError(f"Duplicate input FASTA identifier: {active}")
            result.add(active)
            length = 0
        else:
            if active is None:
                raise ValueError("Sequence before first input FASTA header")
            length += len(line.strip())
    if not result or not length:
        raise ValueError("Input FASTA must contain nonempty sequences")
    return result


def _called_ids(data: bytes, id_column: str) -> set[str]:
    import io
    reader = csv.reader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t", strict=True)
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValueError("Called-ID TSV is empty; include its one-column header") from exc
    if header != [id_column]:
        raise ValueError(f"Called-ID TSV must have exactly one column named {id_column}; export the selected positive IDs, not a full report")
    result = set()
    for row in reader:
        if len(row) != 1:
            raise ValueError("Called-ID TSV must have one ID per row")
        value = _identifier(row[0], "called identifier")
        if value in result:
            raise ValueError(f"Duplicate called identifier: {value}")
        result.add(value)
    return result


def _ratio(numerator: int, denominator: int):
    return numerator / denominator if denominator else None


def compare_reference(reference, predictions, input_path, id_column: str = "protein_id") -> dict:
    """Score one explicit evaluation unit; unlisted calls are never negatives.

    Presence is checked against input FASTA IDs. For locus evaluation, creation
    of the locus-ID FASTA and mapping of calls are curator responsibilities.
    No protein-to-gene fallback or automatic fusion interpretation is performed.
    """
    ref_path, ref_data, ref_sha = _read(reference)
    pred_path, pred_data, pred_sha = _read(predictions)
    fasta_path, fasta_data, fasta_sha = _read(input_path)
    try:
        manifest = json.loads(ref_data, object_pairs_hook=_unique_object)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError("Reference manifest must be valid JSON") from exc
    if not isinstance(manifest, dict):
        raise ValueError("Reference manifest must be a JSON object")
    required = {"schema_version", "benchmark_id", "unit", "prediction_set", "input_sha256",
                "exhaustive_labels", "synthetic", "records"}
    if set(manifest) != required:
        raise ValueError(f"Reference fields must be exactly: {', '.join(sorted(required))}")
    if type(manifest["schema_version"]) is not int or manifest["schema_version"] != 1:
        raise ValueError("Unsupported reference schema_version (expected 1)")
    for field in ("benchmark_id", "prediction_set"):
        _text(manifest[field], field)
    if not isinstance(manifest["unit"], str) or manifest["unit"] not in {"protein", "locus"}:
        raise ValueError("Reference unit must be protein or locus")
    if id_column != manifest["unit"] + "_id":
        raise ValueError(f"Use --id-column {manifest['unit']}_id for this reference unit")
    for field in ("exhaustive_labels", "synthetic"):
        if type(manifest[field]) is not bool:
            raise ValueError(f"{field} must be a JSON boolean")
    if manifest["input_sha256"] != fasta_sha:
        raise ValueError("Reference input_sha256 does not match the supplied evaluation FASTA")
    try:
        universe = _input_ids(fasta_data)
        calls = _called_ids(pred_data, id_column)
    except (UnicodeDecodeError, csv.Error) as exc:
        raise ValueError(f"Invalid input FASTA or called-ID TSV: {exc}") from exc
    outside = calls - universe
    if outside:
        raise ValueError(f"Called IDs absent from evaluation FASTA: {', '.join(sorted(outside))}")
    records = manifest["records"]
    if not isinstance(records, list) or not records:
        raise ValueError("Reference records must be a nonempty list")
    seen = set()
    evaluated = []
    row_fields = {id_column, "expected_call", "evidence_category", "evidence_source", "reference_present", "training_overlap"}
    for row in records:
        if not isinstance(row, dict) or set(row) != row_fields:
            raise ValueError(f"Each reference record must contain exactly: {', '.join(sorted(row_fields))}")
        rid = _identifier(row[id_column], id_column)
        if rid in seen:
            raise ValueError(f"Duplicate reference identifier: {rid}")
        seen.add(rid)
        if not isinstance(row["expected_call"], str) or row["expected_call"] not in _LABELS:
            raise ValueError(f"Invalid expected_call for {rid}")
        if not isinstance(row["evidence_category"], str) or row["evidence_category"] not in _EVIDENCE:
            raise ValueError(f"Invalid evidence_category for {rid}")
        if row["evidence_source"] is not None:
            _text(row["evidence_source"], "evidence_source")
        if row["expected_call"] != "uncertain" and not row["evidence_source"]:
            raise ValueError(f"A labeled record requires evidence_source: {rid}")
        if not isinstance(row["training_overlap"], str) or row["training_overlap"] not in _OVERLAP:
            raise ValueError(f"Invalid training_overlap for {rid}")
        if type(row["reference_present"]) is not bool:
            raise ValueError(f"reference_present must be a JSON boolean: {rid}")
        present = rid in universe
        if row["reference_present"] != present:
            raise ValueError(f"reference_present disagrees with evaluation FASTA: {rid}")
        if row["evidence_category"] == "synthetic" and not manifest["synthetic"]:
            raise ValueError("Synthetic reference records require synthetic=true")
        outcome = "UNAVAILABLE" if not present else {
            "positive": "TP" if rid in calls else "FN",
            "negative": "FP" if rid in calls else "TN",
            "uncertain": "UNRESOLVED_CALLED" if rid in calls else "UNRESOLVED_NOT_CALLED",
        }[row["expected_call"]]
        evaluated.append({**row, "called": rid in calls, "outcome": outcome})
    resolved = {row[id_column] for row in evaluated
                if row["reference_present"] and row["expected_call"] != "uncertain"}
    if manifest["exhaustive_labels"] and resolved != universe:
        raise ValueError("exhaustive_labels=true requires positive/negative labels for every evaluation FASTA ID")
    counts = Counter(row["outcome"] for row in evaluated)
    tp, fn, fp, tn = (counts[name] for name in ("TP", "FN", "FP", "TN"))
    negative_panel = fp + tn
    # A selected panel cannot estimate proteome-wide precision or specificity.
    precision_available = manifest["exhaustive_labels"] and negative_panel > 0
    unknown_calls = sorted(calls - resolved)
    strata = {}
    for field in ("evidence_category", "training_overlap"):
        strata[field] = {}
        for value in sorted({row[field] for row in evaluated}):
            subset = Counter(row["outcome"] for row in evaluated if row[field] == value)
            positive_count = subset["TP"] + subset["FN"]
            strata[field][value] = {
                "positive_present": positive_count, "positive_recovered": subset["TP"],
                "positive_recovery": _ratio(subset["TP"], positive_count),
                "negative_present": subset["FP"] + subset["TN"],
                "negative_called": subset["FP"], "unavailable": subset["UNAVAILABLE"],
            }
    return {
        "schema_version": 1, "status": "COMPARED", "benchmark_id": manifest["benchmark_id"],
        "unit": manifest["unit"], "prediction_set": manifest["prediction_set"],
        "synthetic": manifest["synthetic"], "exhaustive_labels": manifest["exhaustive_labels"],
        "provenance": {
            "reference": {"path": str(ref_path), "sha256": ref_sha},
            "predictions": {"path": str(pred_path), "sha256": pred_sha},
            "input": {"path": str(fasta_path), "sha256": fasta_sha},
        },
        "counts": {"input_units": len(universe), "called_units": len(calls),
                   "positive_present": tp + fn, "positive_recovered": tp,
                   "negative_present": negative_panel, "negative_called": fp,
                   "unavailable_references": counts["UNAVAILABLE"],
                   "unresolved_or_unlisted_calls": len(unknown_calls)},
        "metrics": {"positive_recovery": _ratio(tp, tp + fn),
                    "precision": _ratio(tp, tp + fp) if precision_available else None,
                    "specificity": _ratio(tn, negative_panel) if precision_available else None},
        "precision_scope": "Entire declared evaluation FASTA" if precision_available else None,
        "precision_unavailable_reason": None if precision_available else "Requires exhaustive positive/negative labels for the evaluation FASTA and a negative reference panel",
        "independence": "SYNTHETIC" if manifest["synthetic"] else "NOT_ESTABLISHED_BY_THIS_HARNESS",
        "presence_check": "Exact evaluation FASTA IDs; locus mapping and biological labels are supplied curator assertions",
        "run_integrity": "NOT_CHECKED: verify completed pipeline outputs before exporting the called-ID set",
        "unknown_called_ids": unknown_calls, "strata": strata, "records": evaluated,
    }
