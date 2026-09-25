"""Additive, read-only explanations of saved discovery and rescue decisions.

These exports never select candidates or change their classifications. HMM
statistics come from the original Stage 1 hmmscan tables. Exported sequences
are individual query-alignment spans, not reconstructed or complete domains.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import nullcontext
import csv
import hashlib
import json
import math
from pathlib import Path
import re

from .parsers.fasta import fasta_iter, write_fasta

SCHEMA_VERSION = 1
CHANNELS = {
    "nbd_strict": "nbd_candidate_ids_strict.txt",
    "nbd_relaxed": "nbd_candidate_ids_relaxed.txt",
    "pfam_nbd": "pfam_nbd_candidate_ids.txt",
    "pfam_enriched": "pfam_priority_ids.txt",
    "description": "desc_priority_ids.txt",
    "description_structural_hint": "desc_priority_ids.structural_hint.txt",
}
HIT_COLUMNS = [
    "hit_id", "protein_id", "source", "model_name", "model_accession",
    "model_length", "query_length", "full_evalue", "full_bitscore",
    "domain_i_evalue", "domain_bitscore", "domain_bias", "hmm_start", "hmm_end",
    "model_coverage", "ali_start", "ali_end", "aligned_length", "env_start", "env_end",
    "passes_strict_recorded", "passes_relaxed_recorded", "accepted_pfam_recorded", "accepted_by_recorded_filters",
    "raw_path", "raw_line", "sequence_source", "segment_status", "segment_sha256",
]
CANDIDATE_COLUMNS = [
    "protein_id", "gene_id", "transcript_id", "scaffold", "start", "end", "strand",
    "record_type", "present_in_scanned_proteome", "discovery_channels",
    *CHANNELS, "in_raw_union", "in_candidate_union", "in_final_report", "in_strict_set",
    "tier", "nbd_confidence", "classification_source", "flags", "rescue_priority", "NLR_architecture",
    "has_nbd", "has_nbd_stage1", "has_nbd_pfam", "has_stand_like", "nbd_ok",
    "has_effector", "has_sensor", "has_any_repeat_region", "order_invalid", "non_nlr_annot",
    "fused_into", "fusion_members", "decision_notes", "nbd_library_reported_hits",
    "pfam_nbd_reported_hits", "best_reported_hit_id", "best_model_name",
    "best_domain_i_evalue", "best_domain_bitscore", "best_model_coverage",
    "best_ali_start", "best_ali_end",
    *[f"{track}_{field}" for track in ("priority", "comprehensive") for field in (
        "miniprot_requested", "miniprot_mrna_observed", "exonerate_requested",
        "exonerate_selected", "exonerate_not_selected", "exonerate_run_recorded",
        "exonerate_features_observed")],
]


def _sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "t"}


def _flags(value):
    return sorted({v.strip() for v in re.split(r"[;|,]", str(value or ""))
                   if v.strip() and v.strip().upper() not in {"PASS", "NA", "NAN"}})


class _Saved:
    def __init__(self, root):
        self.root = Path(root)
        self.sources = {}
        self.warnings = []

    def path(self, relative):
        path = self.root / relative
        if relative not in self.sources:
            self.sources[relative] = {"status": "PRESENT" if path.is_file() else "MISSING"}
            if path.is_file():
                self.sources[relative].update(sha256=_sha(path), size_bytes=path.stat().st_size)
        return path

    def ids(self, relative):
        path = self.path(relative)
        return set(path.read_text().splitlines()) - {""} if path.is_file() else None

    def table(self, relative):
        path = self.path(relative)
        if not path.is_file():
            return None
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))


def _key(row, model_key):
    return (row["protein_id"], row[model_key], int(row["ali_start"]),
            int(row["ali_end"]), float(row["dom_i_evalue"]))


def _hit_evidence(saved, proteins):
    detailed = saved.table("results/stage1/nbd_hits_detailed.tsv")
    pfam = saved.table("results/stage1/pfam_nbd_hits.tsv")
    original = {_key(row, "nbd_hmm"): row for row in detailed or []}
    accepted_pfam = {_key(row, "pfam_nbd_hmm") for row in pfam or []}
    rows, segments = [], []
    source_seq = "results/stage0/proteins_clean.faa"
    for source, filename in (("nbd_library", "nbd_whole.domtblout"), ("pfam_nbd", "pfam_nbd.domtblout")):
        relative = "results/stage1/" + filename
        path = saved.path(relative)
        if not path.is_file():
            continue
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                try:
                    if len(fields) < 22:
                        raise ValueError("fewer than 22 columns")
                    ml, ql = int(fields[2]), int(fields[5])
                    hs, he, start, end, es, ee = map(int, fields[15:21])
                    ev, score, bias = map(float, (fields[12], fields[13], fields[14]))
                    full_ev, full_score = float(fields[6]), float(fields[7])
                    if not all(math.isfinite(v) for v in (ev, score, bias, full_ev, full_score)):
                        raise ValueError("nonfinite statistic")
                    if ev < 0 or full_ev < 0:
                        raise ValueError("negative E-value")
                except (ValueError, IndexError) as exc:
                    saved.warnings.append(f"Skipped unreadable HMMER evidence at {relative}:{line_number}: {exc}")
                    continue
                pid, model, accession = fields[3], fields[0], fields[1]
                key = (pid, model, start, end, ev)
                prior = original.get(key) if source == "nbd_library" else None
                hit_id = f"{source}_{line_number:08d}"
                seq = proteins.get(pid)
                accepted = ((_bool(prior["passes_strict"]) or _bool(prior["passes_relaxed"]))
                            if prior is not None else None)
                if source == "pfam_nbd":
                    accepted = key in accepted_pfam if pfam is not None else None
                status = "EXPORTED"
                if accepted is None:
                    status = "NOT_ASSESSED"
                elif not accepted:
                    status = "NOT_ACCEPTED_BY_RECORDED_FILTERS"
                elif seq is None:
                    status = "MISSING_SCANNED_PROTEIN"
                elif len(seq) != ql:
                    status = "QUERY_LENGTH_MISMATCH"
                elif not 1 <= start <= end <= len(seq):
                    status = "INVALID_QUERY_COORDINATES"
                elif not 1 <= hs <= he <= ml:
                    status = "INVALID_MODEL_COORDINATES"
                segment = seq[start - 1:end] if status == "EXPORTED" else None
                if segment is not None:
                    segments.append((hit_id, segment))
                row = dict(zip(HIT_COLUMNS[:21], [
                    hit_id, pid, source, model, "" if accession == "-" else accession,
                    ml, ql, full_ev, full_score, ev, score, bias, hs, he,
                    (he - hs + 1) / ml if 1 <= hs <= he <= ml else "",
                    start, end, abs(end - start) + 1, es, ee,
                    _bool(prior["passes_strict"]) if prior is not None else "",
                ]))
                row.update(
                    passes_relaxed_recorded=_bool(prior["passes_relaxed"]) if prior is not None else "",
                    accepted_pfam_recorded=(key in accepted_pfam) if source == "pfam_nbd" and pfam is not None else "",
                    accepted_by_recorded_filters=accepted if accepted is not None else "",
                    raw_path=relative, raw_line=line_number, sequence_source=source_seq,
                    segment_status=status,
                    segment_sha256=hashlib.sha256(segment.encode()).hexdigest() if segment is not None else "",
                )
                rows.append(row)
    return rows, segments


def _rescue(saved):
    data, summary = {}, {}
    for track in ("priority", "comprehensive"):
        tail = "/comprehensive" if track == "comprehensive" else ""
        five, six = f"results/stage5{tail}", f"results/stage6{tail}"
        values = {
            "miniprot_requested": saved.ids(five + f"/queries/miniprot_{track}_ids.txt"),
            "exonerate_requested": saved.ids(five + "/queries/exonerate_refine_ids.txt"),
            "exonerate_selected": saved.ids(six + "/queries/exonerate_queries.ids"),
            "exonerate_run_recorded": saved.ids(six + "/queries/exonerate_run.ids"),
        }
        features = saved.table(five + "/tables/miniprot_features.tsv")
        values["miniprot_mrna_observed"] = ({r["query"] for r in features if r.get("feature") == "mRNA"}
                                             if features is not None else None)
        features = saved.table(six + "/tables/exonerate_summary.tsv")
        values["exonerate_features_observed"] = ({r.get("query", r.get("protein_id", "")) for r in features}
                                                  if features is not None else None)
        requested, selected = values["exonerate_requested"], values["exonerate_selected"]
        values["exonerate_not_selected"] = (requested - selected if requested is not None and selected is not None else None)
        data[track] = values
        summary[track] = {key: len(ids) if ids is not None else None for key, ids in values.items()}
        summary[track]["selection_note"] = "Requested original queries absent from the saved selection were not attempted; this includes the configured per-track cap. Run IDs may also include fusion queries."
    return data, summary


def _membership(pid, values):
    return pid in values if values is not None else ""


def build_evidence(root, settings=None):
    """Read saved stage products and return deterministic tables and a summary.

    Missing source files remain unknown, never zero evidence. Stage-1 pass
    flags are copied from original calls, not recomputed from current defaults.
    """
    saved = _Saved(root)
    path = saved.path("results/stage0/proteins_clean.faa")
    with path.open() if path.is_file() else nullcontext([]) as handle:
        proteins = dict(fasta_iter(handle))
    channels = {key: saved.ids("results/stage1/" + name) for key, name in CHANNELS.items()}
    raw_union = saved.ids("results/stage1/union_candidate_ids.raw.txt")
    union = saved.ids("results/stage1/union_candidate_ids.txt")
    final = saved.table("final_results/nlr_final_report.tsv")
    strict = saved.table("final_results/nlr_strict_candidates.tsv")
    tiered = saved.table("results/stage3/tables/tiered_candidates.tsv")
    tiered_by_id = {row["protein_id"]: row for row in tiered or []}
    final_by_id = {row["protein_id"]: row for row in final or []}
    strict_ids = {row["protein_id"] for row in strict or []} if strict is not None else None
    hits, segments = _hit_evidence(saved, proteins)
    by_protein = defaultdict(list)
    for hit in hits:
        by_protein[hit["protein_id"]].append(hit)
    rescue, rescue_summary = _rescue(saved)
    observed = set(final_by_id) | set(by_protein) | set(raw_union or []) | set(union or [])
    for ids in channels.values():
        observed.update(ids or [])
    rows, review = [], []
    for pid in sorted(observed):
        call = final_by_id.get(pid, {})
        fusion = _bool(call.get("is_fusion_model")) or str(call.get("tier", "")).startswith("FUSION")
        record_type = "FUSION_HYPOTHESIS" if fusion else "ORIGINAL_CANDIDATE" if call else "DISCOVERY_EVIDENCE_ONLY"
        notes = []
        if not call and union is not None and pid not in union:
            notes.append("NOT_IN_CANDIDATE_UNION")
        if path.is_file() and pid not in proteins and not fusion:
            notes.append("ABSENT_FROM_SCANNED_PROTEOME")
        if call and strict_ids is not None and pid not in strict_ids:
            notes.append("NOT_IN_STRICT_SET")
        if call.get("fused_into"):
            notes.append("MEMBER_OF_FUSION_HYPOTHESIS")
        if _flags(call.get("flags")):
            notes.append("FLAGS_REQUIRE_INSPECTION")
        if fusion:
            notes.append("FUSION_REQUIRES_REVIEW")
        elif "NEEDS_REVIEW" in call.get("tier", "") or call.get("tier", "").startswith("TIER_2"):
            notes.append("REVIEW_OR_RESCUE_TIER")
        classified = {**tiered_by_id.get(pid, {}), **call}
        row = {column: classified.get(column, "") for column in CANDIDATE_COLUMNS}
        row.update(protein_id=pid, record_type=record_type,
                   classification_source=("final_results/nlr_final_report.tsv" if fusion else
                                          "results/stage3/tables/tiered_candidates.tsv" if pid in tiered_by_id else ""),
                   present_in_scanned_proteome=pid in proteins if path.is_file() else "",
                   discovery_channels=";".join(key for key, ids in channels.items() if pid in (ids or set())),
                   in_raw_union=_membership(pid, raw_union), in_candidate_union=_membership(pid, union),
                   in_final_report=bool(call) if final is not None else "", in_strict_set=_membership(pid, strict_ids),
                   decision_notes=";".join(notes))
        row.update({key: _membership(pid, ids) for key, ids in channels.items()})
        candidates = by_protein[pid]
        for source, filename in (("nbd_library", "nbd_whole"), ("pfam_nbd", "pfam_nbd")):
            available = saved.sources[f"results/stage1/{filename}.domtblout"]["status"] == "PRESENT"
            row[source + "_reported_hits"] = sum(h["source"] == source for h in candidates) if available else ""
        if candidates:
            best = min(candidates, key=lambda h: (h["domain_i_evalue"], -h["domain_bitscore"], -h["aligned_length"], h["hit_id"]))
            for target, source in (("hit_id", "hit_id"), ("model_name", "model_name"),
                                   ("domain_i_evalue", "domain_i_evalue"), ("domain_bitscore", "domain_bitscore"),
                                   ("model_coverage", "model_coverage"), ("ali_start", "ali_start"), ("ali_end", "ali_end")):
                row["best_reported_hit_id" if target == "hit_id" else "best_" + target] = best[source]
        for track, values in rescue.items():
            row.update({f"{track}_{key}": _membership(pid, ids) for key, ids in values.items()})
        rows.append(row)
        if notes:
            review.append(row)
    original_final = [r for r in final or [] if not _bool(r.get("is_fusion_model")) and not r.get("tier", "").startswith("FUSION")]
    genes = Counter(r.get("gene_id") for r in original_final if r.get("gene_id"))
    locus_keys = ("scaffold", "start", "end", "strand")
    annotated_loci = {tuple(r.get(k) for k in locus_keys) for r in original_final if all(r.get(k) for k in locus_keys)}
    qc = saved.table("final_results/qc_summary.tsv")
    qc_map = {r["metric"]: r["value"] for r in qc or []}
    reconciliation = []
    for metric, observed_count in (("n_total_candidates", len(original_final)),
                                    ("n_strict_candidates_combined", len(strict or []))):
        value = qc_map.get(metric)
        try:
            matches = float(value) == observed_count if value is not None else None
        except ValueError:
            matches = False
        reconciliation.append({"metric": metric, "recorded_value": value,
                               "observed_value": observed_count if final is not None and strict is not None else None,
                               "matches": matches})
    for values in rescue_summary.values():
        values["configured_original_query_cap"] = (settings or {}).get("EXON_MAXN")
    model_inputs = {}
    manifest_path = saved.root / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for name in ("NBD_HMMS", "PFAM_DB"):
            record = manifest.get("inputs", {}).get(name, {})
            model_inputs[name] = {"sha256": record.get("sha256"), "source": "recorded run manifest; not rehashed"}
    saved.path("results/stage1/pfam_nbd.hmm")
    run_info_path = saved.path("results/stage1/run_info.txt")
    run_info = dict(line.split(": ", 1) for line in run_info_path.read_text().splitlines() if ": " in line) if run_info_path.is_file() else {}
    summary = {
        "schema_version": SCHEMA_VERSION,
        "scope": "Saved Stage 1 discovery evidence, original final/strict membership, and recorded rescue queries; no biological decisions are recomputed.",
        "sequence_scope": "One unaligned sequence segment per Stage 1 hit with recorded strict, relaxed, or Pfam acceptance and valid coordinates, using 1-based inclusive query alignment coordinates on proteins_clean.faa. Segments may overlap or include competing-family matches and are not complete NBDs, fusion proteins, or a multiple sequence alignment. Rejected or unassessed hits remain in the mapping table without sequences.",
        "candidate_table_scope": "Final-report IDs plus IDs observed in discovery channels or reported Stage 1 hits; not all unselected proteome proteins.",
        "best_hit_scope": "Lowest reported domain i-Evalue, then highest domain bitscore and aligned length; this display choice is not a candidate-selection or confidence rule. E-values from different database searches have different search spaces.",
        "counts": {"evidence_ids": len(rows), "review_or_excluded_ids": len(review),
                   "final_report_rows": len(final) if final is not None else None,
                   "original_final_candidates": len(original_final) if final is not None else None,
                   "strict_report_rows": len(strict) if strict is not None else None,
                   "reported_hits": len(hits), "accepted_reported_hits": sum(h["accepted_by_recorded_filters"] is True for h in hits),
                   "unassessed_reported_hits": sum(h["accepted_by_recorded_filters"] == "" for h in hits),
                   "exported_segments": len(segments)},
        "scan_status": {"pfam_nbd_requested": _bool(run_info["pfam_nbd_scan"]) if "pfam_nbd_scan" in run_info else None,
                        "pfam_nbd_completed": _bool(run_info["pfam_nbd_scan_completed"]) if "pfam_nbd_scan_completed" in run_info else None,
                        "note": "Recorded Stage 1 execution status. An empty hit table can reflect a disabled or uncompleted optional scan; it does not establish biological absence."},
        "discovery_channels": {key: len(ids) if ids is not None else None for key, ids in channels.items()},
        "channel_overlaps": [{"left": left, "right": right,
                              "shared_ids": len(channels[left] & channels[right]) if channels[left] is not None and channels[right] is not None else None}
                             for i, left in enumerate(CHANNELS) for right in list(CHANNELS)[i + 1:]],
        "annotation_units": {"original_candidate_rows": len(original_final) if final is not None else None,
                             "distinct_nonempty_gene_ids": len(genes) if final is not None else None,
                             "gene_ids_with_multiple_candidate_rows": sum(n > 1 for n in genes.values()) if final is not None else None,
                             "distinct_recorded_coordinate_tuples": len(annotated_loci) if final is not None else None,
                             "note": "Existing annotation identifiers and exact coordinate tuples; not orthology, allele, or biological locus inference."},
        "union": {"raw": len(raw_union) if raw_union is not None else None,
                  "selected": len(union) if union is not None else None,
                  "raw_not_selected": len(raw_union - union) if raw_union is not None and union is not None else None},
        "final_tier_counts": dict(sorted(Counter(r.get("tier", "") for r in final or []).items())),
        "segment_status_counts": dict(sorted(Counter(r["segment_status"] for r in hits).items())),
        "rescue": rescue_summary, "qc_reconciliation": reconciliation, "model_inputs": model_inputs,
        "sources": saved.sources,
        "discovery_settings": settings or {}, "warnings": saved.warnings,
    }
    missing = [name for name, entry in saved.sources.items() if entry["status"] == "MISSING" and name not in {"manifest.json", "results/stage1/pfam_nbd.hmm"}]
    summary["availability"] = "COMPLETE" if not missing else "PARTIAL"
    summary["missing_sources"] = sorted(missing)
    return rows, review, hits, segments, summary


def export_evidence(root, destinations, settings=None):
    """Write the same additive products to explicitly supplied destinations."""
    rows, review, hits, segments, summary = build_evidence(root, settings)
    outputs = []
    for destination in destinations:
        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=True)
        for name, columns, values in (("candidate_evidence.tsv", CANDIDATE_COLUMNS, rows),
                                      ("review_and_exclusions.tsv", CANDIDATE_COLUMNS, review),
                                      ("nbd_hit_evidence.tsv", HIT_COLUMNS, hits)):
            path = destination / name
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(values)
            outputs.append(path)
        path = destination / "nbd_aligned_segments.faa"
        with path.open("w") as handle:
            write_fasta(iter(segments), handle)
        outputs.append(path)
        path = destination / "evidence_summary.json"
        path.write_text(json.dumps(summary, sort_keys=True, indent=2, allow_nan=False) + "\n")
        outputs.append(path)
    return outputs
