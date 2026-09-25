#!/usr/bin/env python3
"""Strict, read-only comparison of FuNLR scientific output artifacts.

No sorting, numerical tolerance, identifier rewriting, NA substitution, or
sequence rewrapping is performed. Only the explicitly declared metadata below
is excluded. Standard-library-only; run from a downloaded source archive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import sys


TIER_FASTAS = {
    "tier1A": "tier1A_high_confidence", "tier1B": "tier1B_needs_review",
    "tier2A": "tier2A_high_priority_rescue", "tier2B": "tier2B_rescue_candidate",
    "tier3": "tier3_architectural_variant", "tier4A": "tier4A_repeat_only_no_nbd",
    "tier4B": "tier4B_likely_stand_housekeeping", "tier4C": "tier4C_no_nbd_no_repeat_low_signal",
    "all_tiered_candidates": "all_tiered_candidates",
}
REQUIRED_FILES = {
    "results/stage0/contig_lengths.tsv", "results/stage0/proteins_clean.faa",
    "results/stage0/protein_lengths.tsv", "results/stage0/id_coords.tsv",
    "results/stage0/eggnog_light.tsv", "results/stage0/master_table.tsv",
    "results/stage1/nbd_hits.tsv", "results/stage1/nbd_candidate_ids.txt",
    "results/stage1/nbd_parse_report.txt", "results/stage1/eggnog_priority_ids.txt",
    "results/stage1/union_candidate_ids.txt", "results/stage1/union_candidates.faa",
    "results/stage2/stage1_nbd.ids", "results/stage2/all_domain_hits.tsv",
    "results/stage2/architecture_summary.tsv", "results/stage2/parse_filter_report.tsv",
    "results/stage3/tables/tiered_candidates.tsv", "results/stage3/tables/rescue_priority.tsv",
    "results/stage3/beds/nlr_candidates.bed", "results/stage3/summaries/tier_summary.tsv",
    "results/stage3/summaries/flag_counts.tsv", "results/stage4/fasta_export_summary.tsv",
    "results/stage5/queries/rescue_queries.raw.faa", "results/stage5/queries/rescue_queries.faa",
    "results/stage5/queries/exonerate_refine_ids.txt", "results/stage5/tables/miniprot_features.tsv",
    "results/stage5/tables/miniprot_summary.tsv", "results/stage5/tables/miniprot_hit_counts.tsv",
    "results/stage5/summaries/query_counts.tsv", "results/stage6/queries/exonerate_refine_ids.work.txt",
    "results/stage6/queries/exonerate_queries.faa", "results/stage6/queries/exonerate_queries.ids",
    "results/stage6/tables/exonerate_summary.tsv", "results/stage6/summaries/exonerate_counts.tsv",
    "results/stage7/tables/nlr_final_report.tsv", "results/stage7/tables/nlr_strict_candidates.tsv",
    "results/stage7/summaries/tier_summary.tsv", "results/stage7/summaries/qc_summary.tsv",
}
REQUIRED_FILES.update(f"results/stage4/{stem}.faa" for stem in TIER_FASTAS.values())
REQUIRED_FILES.update(f"results/stage4/{label}.ids" for label in TIER_FASTAS)
for stem in ("nlr_candidates", "nlr_miniprot_loci", "nlr_exonerate_loci"):
    REQUIRED_FILES.update({f"results/stage7/beds/{stem}.bed", f"final_results/{stem}.bed"})
REQUIRED_FILES.update(f"final_results/{stem}.tsv" for stem in (
    "nlr_final_report", "nlr_strict_candidates", "tier_summary", "qc_summary",
))

# Historical and v3 layouts are both understood. Cross-version comparisons
# retain their schema/path differences as differences; no conversion occurs.
LEGACY_REQUIRED_FILES = REQUIRED_FILES.copy()
REQUIRED_FILES = {path for path in LEGACY_REQUIRED_FILES
                  if not path.startswith("results/stage4/") and path != "results/stage1/eggnog_priority_ids.txt"}
REQUIRED_FILES.update({
    "results/stage0/id_mapping_report.tsv", "results/stage1/nbd_hits_detailed.tsv",
    "results/stage1/nbd_candidate_ids_strict.txt", "results/stage1/nbd_candidate_ids_relaxed.txt",
    "results/stage1/pfam_nbd_hits.tsv", "results/stage1/pfam_nbd_candidate_ids.txt",
    "results/stage1/nbd_hits_with_master_strict.tsv", "results/stage1/nbd_hits_with_master_relaxed.tsv",
    "results/stage1/desc_priority_ids.txt", "results/stage1/desc_priority_ids.structural_hint.txt",
    "results/stage1/pfam_priority_ids.txt", "results/stage1/pfam_enriched_top50.txt",
    "results/stage1/pfam_enriched_in_nbd.tsv", "results/stage1/union_candidate_ids.raw.txt",
    "results/stage1/union_candidate_ids.missing_from_fasta.txt", "results/stage1/union_candidate_accounting.tsv",
    "results/stage2/stage1_nbd_hits.tsv", "results/stage2/asm_hits_summary.tsv",
    "results/stage3/summaries/universe_summary.tsv", "results/stage4/fasta_export_summary.tsv",
    "results/stage4/tiered_candidates.faa", "results/stage4/tiered_candidates.ids",
    "results/stage4/rescue_priority.faa", "results/stage4/rescue_priority.ids",
    "results/stage5/queries/miniprot_priority_ids.txt", "results/stage5/qc/exonerate_refine_reasons.tsv",
    "results/stage5/comprehensive/queries/miniprot_comprehensive_ids.txt",
    "results/stage5/comprehensive/queries/comprehensive_queries.faa",
    "results/stage5/comprehensive/queries/comprehensive_queries.dedup.faa",
})
for relative in ("queries/exonerate_refine_ids.txt", "tables/miniprot_features.tsv", "tables/miniprot_summary.tsv",
                 "tables/miniprot_hit_counts.tsv", "summaries/query_counts.tsv", "qc/exonerate_refine_reasons.tsv"):
    REQUIRED_FILES.add("results/stage5/comprehensive/" + relative)
for prefix in ("results/stage6/", "results/stage6/comprehensive/"):
    REQUIRED_FILES.update(prefix + relative for relative in (
        "queries/exonerate_refine_ids.work.txt", "queries/exonerate_queries.faa", "queries/exonerate_queries.ids",
        "queries/exonerate_all.faa", "queries/exonerate_all.ids", "queries/exonerate_run.ids",
        "queries/fusion_queries.faa", "queries/fusion_refine_ids.txt", "tables/exonerate_summary.tsv",
        "tables/fusion_manifest.tsv", "tables/fusion_evidence.tsv", "summaries/exonerate_counts.tsv"))
for prefix in ("results/stage7/tables/", "final_results/"):
    REQUIRED_FILES.update(prefix + "nlr_strict_candidates." + track + ".tsv" for track in ("priority", "comprehensive"))
for prefix in ("results/stage7/beds/", "final_results/"):
    REQUIRED_FILES.update(prefix + name + ".bed" for name in ("nlr_comprehensive_miniprot", "nlr_comprehensive_exonerate", "nlr_fusions_exonerate"))
REQUIRED_FILES.add("results/stage7/tables/fusion_evidence.tsv")


# These are path-only manifests, not classification or sequence data.
EXCLUDED_FILES = {
    "results/stage0/input_manifest.tsv": "input paths/sizes; root manifests retain identities outside this scientific comparison",
    "results/stage0/software_manifest.tsv": "tool/database provenance is outside this scientific comparison",
    "results/stage5/comprehensive/comprehensive_manifest.tsv": "artifact paths; scientific artifacts compared directly",
    "results/stage6/comprehensive/stage6_manifest.tsv": "artifact paths; scientific artifacts compared directly",
    "results/stage4/fasta_export_manifest.tsv": "export paths; actual FASTAs and ID files are compared",
    **{f"results/stage{n}/stage{n}_manifest.tsv": "artifact paths; scientific artifacts are compared directly"
       for n in (5, 6, 7)},
}
# Only the values in this named column are masked; its header and column order,
# all other cells, quoting outside the masked cells, and line endings stay exact.
EXCLUDED_COLUMNS = {
    "results/stage5/summaries/query_counts.tsv": ["file"],
    "results/stage5/comprehensive/summaries/query_counts.tsv": ["file"],
    "results/stage1/priority_sources_summary.tsv": ["file"],
    "final_results/integration/export_status.tsv": ["path"],
}
EXCLUDED_DIRECTORIES = {"results/stage7/plot_work"}
SCIENTIFIC_SUFFIXES = {".tsv", ".faa", ".fa", ".fasta", ".bed", ".ids"}
EXCLUSION_RULES = [
    "Only results/stage0 through results/stage7 and final_results are scanned.",
    "Within those directories, TSV/FASTA/BED/IDS and the required named TXT files are compared.",
    "Run metadata, timestamps, command logs, work indexes, raw HMMER/miniprot/Exonerate output and other file types are excluded.",
    "Parsed domain/alignment TSVs and rescue BEDs are compared; this does not audit raw alignment output.",
    "Metadata directories and Stage 7 plot_work duplicates are excluded; exported domain summaries are compared.",
    "V3 dynamic tier TSVs imply required matching Stage 4 ID/FASTA exports. Optional integration outputs are compared when present in either run.",
    "Tagged integration GFF3/eggNOG copies and gene ID TXT lists are compared as exact bytes; images, R session logs and raw native alignments are excluded.",
]
MIN_COLUMNS = {
    "master_table.tsv": {"protein_id", "scaffold", "start", "end", "protein_length"},
    "protein_lengths.tsv": {"protein_id", "protein_length"},
    "id_coords.tsv": {"protein_id", "scaffold", "start", "end", "strand"},
    "eggnog_light.tsv": {"protein_id"},
    "nbd_hits.tsv": {"protein_id", "nbd_hmm", "dom_i_evalue", "ali_start", "ali_end", "ali_len"},
    "architecture_summary.tsv": {"protein_id", "has_nbd"},
    "all_domain_hits.tsv": {"protein_id", "domain", "dom_i_evalue", "start", "end"},
    "tiered_candidates.tsv": {"protein_id", "tier", "flags", "rescue_priority"},
    "rescue_priority.tsv": {"protein_id", "tier", "flags", "rescue_priority"},
    "nlr_final_report.tsv": {"protein_id", "tier", "flags"},
    "nlr_strict_candidates.tsv": {"protein_id", "tier", "flags"},
    "miniprot_summary.tsv": {"scaffold", "start", "end"},
    "exonerate_summary.tsv": {"protein_id", "scaffold", "start", "end"},
    "fasta_export_summary.tsv": {"file", "bytes", "seqs"},
    "query_counts.tsv": {"label", "file", "seqs"},
    "qc_summary.tsv": {"metric", "value"},
    "tier_summary.tsv": {"tier", "count"},
}


def required_for(root: Path) -> set[str]:
    architecture = root / "results/stage2/architecture_summary.tsv"
    header = []
    if architecture.is_file():
        with architecture.open() as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
    legacy = "has_repeat" in header and "has_sensor" not in header
    required = (LEGACY_REQUIRED_FILES if legacy else REQUIRED_FILES).copy()
    if not legacy:
        for table in (root / "results/stage3/tables").glob("TIER_*.tsv"):
            required.update({f"results/stage4/{table.stem}.ids", f"results/stage4/{table.stem}.faa"})
        tiered = root / "results/stage3/tables/tiered_candidates.tsv"
        if tiered.is_file():
            with tiered.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    tier = row.get("tier", "") or ""
                    if re.fullmatch(r"TIER_[A-Za-z0-9_]+", tier):
                        required.update({f"results/stage3/tables/{tier}.tsv",
                                         f"results/stage4/{tier}.ids", f"results/stage4/{tier}.faa"})
    return required


def selected(relative: str, path: Path) -> bool:
    if relative in EXCLUDED_FILES or "metadata" in relative.split("/") or any(relative.startswith(folder + "/") for folder in EXCLUDED_DIRECTORIES):
        return False
    if "/integration/" in relative and path.suffix.lower() in {".gff3", ".annotations", ".txt"}:
        return True
    return path.suffix.lower() in SCIENTIFIC_SUFFIXES or relative in (REQUIRED_FILES | LEGACY_REQUIRED_FILES)


def discover(root: Path) -> set[str]:
    found = set()
    for folder in [root / "results" / f"stage{n}" for n in range(8)] + [root / "final_results"]:
        for path in folder.rglob("*"):
            relative = path.relative_to(root).as_posix()
            if path.is_file() and selected(relative, path):
                found.add(relative)
    return found


def mask_metadata_cells(text: str, column_indices: set[int]) -> str:
    """Mask excluded TSV field spans without rewriting any other bytes."""
    pieces, start, row, col, quoted, pos = [], 0, 0, 0, False, 0
    while pos < len(text):
        char = text[pos]
        if char == '"':
            if quoted and pos + 1 < len(text) and text[pos + 1] == '"':
                pos += 2
                continue
            if quoted or pos == start:
                quoted = not quoted
        if not quoted and char in "\t\r\n":
            pieces.append("<EXCLUDED_PATH>" if row and col in column_indices else text[start:pos])
            if char == "\r" and pos + 1 < len(text) and text[pos + 1] == "\n":
                pieces.append("\r\n")
                pos += 1
            else:
                pieces.append(char)
            if char == "\t":
                col += 1
            else:
                row, col = row + 1, 0
            start = pos + 1
        pos += 1
    pieces.append("<EXCLUDED_PATH>" if start < len(text) and row and col in column_indices else text[start:])
    return "".join(pieces)


def inspect_file(path: Path, relative: str) -> dict:
    raw = path.read_bytes()
    text = raw.decode("utf-8")
    result = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw), "compared": raw}
    suffix = path.suffix.lower()
    if suffix == ".tsv":
        rows = list(csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True))
        if not rows:
            raise ValueError("TSV is empty; headed empty tables must retain a schema")
        headerless = relative == "results/stage0/contig_lengths.tsv"
        if not headerless and (not rows[0] or not any(rows[0]) or len(set(rows[0])) != len(rows[0])):
            raise ValueError("TSV header is blank or has duplicate column names")
        width = 2 if headerless else len(rows[0])
        if any(len(row) != width for row in rows):
            raise ValueError("TSV rows have inconsistent column counts")
        if not headerless:
            missing = MIN_COLUMNS.get(path.name, set()) - set(rows[0])
            if missing:
                raise ValueError("TSV lacks required columns: " + ", ".join(sorted(missing)))
            alternatives = {
                "architecture_summary.tsv": ({"has_repeat", "order_ok"}, {"has_sensor", "order_invalid"}),
                "miniprot_summary.tsv": ({"protein_id"}, {"query"}),
            }.get(path.name)
            if alternatives and not any(option <= set(rows[0]) for option in alternatives):
                raise ValueError("TSV does not match a recognized legacy or v3 schema")
            result["header"] = rows[0]
        result["records"] = rows if headerless else rows[1:]
        columns = EXCLUDED_COLUMNS.get(relative, [])
        if columns:
            indices = {rows[0].index(name) for name in columns if name in rows[0]}
            result["compared"] = mask_metadata_cells(text, indices).encode("utf-8")
            result["records"] = [["<EXCLUDED_PATH>" if i in indices else value for i, value in enumerate(row)]
                                 for row in rows[1:]]
    elif suffix in {".faa", ".fa", ".fasta"}:
        headers, seen = [], False
        for line in text.splitlines():
            if line.startswith(">"):
                if not line[1:].strip():
                    raise ValueError("FASTA has an empty header")
                headers.append(line)
                seen = True
            elif line and not seen:
                raise ValueError("FASTA contains sequence data before its first header")
        result["headers"] = headers
    elif suffix == ".bed":
        rows = [line.split("\t") for line in text.splitlines()]
        for row in rows:
            if len(row) != 6:
                raise ValueError("FuNLR BED must have exactly six columns")
            if int(row[1]) < 0 or int(row[2]) < int(row[1]):
                raise ValueError("BED has invalid start/end coordinates")
        result["records"] = rows
    return result


def difference_details(before: dict, after: dict) -> dict:
    left, right = before["compared"], after["compared"]
    offset = next((i for i, (a, b) in enumerate(zip(left, right)) if a != b), min(len(left), len(right)))
    details = {"first_differing_byte": offset, "reference_line": left[:offset].count(b"\n") + 1}
    for label in ("header", "headers"):
        if label in before and before[label] != after.get(label):
            details[label + "_changed"] = True
    if "records" in before and "records" in after:
        a, b = before["records"], after["records"]
        details["reference_rows"], details["candidate_rows"] = len(a), len(b)
        details["first_differing_row"] = next(
            (i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y),
            min(len(a), len(b)) + 1 if len(a) != len(b) else None,
        )
        details["format_only"] = a == b and before.get("header") == after.get("header")
    return details


def compare_runs(reference: Path, candidate: Path, *, required_files=None) -> dict:
    reference, candidate = Path(reference).resolve(), Path(candidate).resolve()
    if not reference.is_dir() or not candidate.is_dir():
        raise ValueError("Both arguments must be existing FuNLR run directories")
    if reference == candidate:
        raise ValueError("Reference and candidate must be different run directories")
    required = (required_for(reference) | required_for(candidate)) if required_files is None else set(required_files)
    paths = sorted(required | discover(reference) | discover(candidate))
    if not paths:
        raise ValueError("No scientific artifacts were selected for comparison")
    differences, matches = [], 0
    for relative in paths:
        first, second = reference / relative, candidate / relative
        if not first.is_file() or not second.is_file():
            differences.append({"path": relative, "status": "missing",
                                "missing_from": [label for label, path in (("reference", first), ("candidate", second)) if not path.is_file()]})
            continue
        details = {}
        for label, path in (("reference", first), ("candidate", second)):
            try:
                details[label] = inspect_file(path, relative)
            except (OSError, UnicodeError, ValueError, csv.Error) as exc:
                details[label] = {"error": str(exc)}
        if any("error" in entry for entry in details.values()):
            differences.append({"path": relative, "status": "invalid", **{
                label + "_error": item["error"] for label, item in details.items() if "error" in item}})
        elif details["reference"]["compared"] != details["candidate"]["compared"]:
            differences.append({"path": relative, "status": "different", **{
                label + "_sha256": item["sha256"] for label, item in details.items()},
                **difference_details(details["reference"], details["candidate"])})
        else:
            matches += 1
    return {
        "schema_version": 1, "status": "different" if differences else "match",
        "reference": str(reference), "candidate": str(candidate), "checked_files": len(paths),
        "matching_files": matches, "differences": differences,
        "comparison_policy": "Exact bytes; only declared metadata cells/files are excluded. No scientific values or ordering are normalized.",
        "exclusions": {"rules": EXCLUSION_RULES, "files": EXCLUDED_FILES, "columns": EXCLUDED_COLUMNS, "directories": sorted(EXCLUDED_DIRECTORIES)},
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="completed baseline run directory containing results/ and final_results/")
    parser.add_argument("candidate", type=Path, help="completed run directory to compare")
    parser.add_argument("--json", type=Path, dest="report", help="save the JSON report here; default prints JSON to stdout")
    args = parser.parse_args(argv)
    try:
        if args.report:
            report_path = args.report.resolve()
            for root in (args.reference.resolve(), args.candidate.resolve()):
                if report_path.is_relative_to(root):
                    raise ValueError("Write the comparison report outside both run directories to preserve their artifacts and provenance")
        report = compare_runs(args.reference, args.candidate)
        result = 0 if report["status"] == "match" else 1
    except (OSError, ValueError) as exc:
        report, result = {"schema_version": 1, "status": "error", "error": str(exc)}, 2
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        if result == 2:
            print(rendered, file=sys.stderr, end="")
            return result
        try:
            args.report.write_text(rendered)
        except OSError as exc:
            print(f"Cannot write comparison report: {exc}", file=sys.stderr)
            return 2
        print(f"{report['status']}: {args.report}", file=sys.stderr)
    else:
        print(rendered, end="")
    return result


if __name__ == "__main__":
    raise SystemExit(main())
