"""Evidence exports must explain saved calls without changing their meaning."""
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from funlr.evidence import CHANNELS, build_evidence, export_evidence
from funlr.parsers.fasta import fasta_iter


def text(root, name, value):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value)
    return path


def table(root, name, rows, columns):
    path = text(root, name, "")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def raw(pid, model="NACHT", start=3, end=12, ev="1e-9", qlen=40, hs=2, he=10):
    return f"{model} PF00001.2 20 {pid} - {qlen} 1e-8 42.0 0.1 1 1 1e-10 {ev} 38.5 0.2 {hs} {he} {start} {end} 1 15 0.9 model description\n"


def fixture(root, empty=False):
    seq = "ACDEFGHIKLMNPQRSTVWY" * 2
    text(root, "results/stage0/proteins_clean.faa", ">p1\n" + seq + "\n>p2\n" + seq[::-1] + "\n>excluded\n" + seq + "\n")
    for key, filename in CHANNELS.items():
        text(root, "results/stage1/" + filename, "" if empty else ("p1\np2\n" if key in {"nbd_relaxed", "pfam_enriched"} else "p1\n" if key in {"nbd_strict", "pfam_nbd"} else "excluded\n"))
    for name in ("union_candidate_ids.raw.txt", "union_candidate_ids.txt"):
        text(root, "results/stage1/" + name, "" if empty else "p1\np2\n")
    custom = [] if empty else [raw("p1"), raw("p1", start=23, end=32), raw("p2", ev="1e-2"), raw("excluded", ev="0.5")]
    text(root, "results/stage1/nbd_whole.domtblout", "# original hmmscan\n" + "".join(custom))
    text(root, "results/stage1/pfam_nbd.domtblout", "" if empty else raw("p1"))
    text(root, "results/stage1/run_info.txt", "pfam_nbd_scan: 1\npfam_nbd_scan_completed: True\n")
    detailed = [] if empty else [dict(protein_id=pid, nbd_hmm="NACHT", dom_i_evalue=ev, ali_start=start, ali_end=end,
                                    passes_strict=pid == "p1", passes_relaxed=pid in {"p1", "p2"})
                               for pid, start, end, ev in (("p1",3,12,"1e-9"),("p1",23,32,"1e-9"),("p2",3,12,"1e-2"),("excluded",3,12,"0.5"))]
    table(root, "results/stage1/nbd_hits_detailed.tsv", detailed,
          ["protein_id","nbd_hmm","dom_i_evalue","ali_start","ali_end","passes_strict","passes_relaxed"])
    table(root, "results/stage1/pfam_nbd_hits.tsv", [] if empty else [dict(protein_id="p1",pfam_nbd_hmm="NACHT",dom_i_evalue="1e-9",ali_start=3,ali_end=12)],
          ["protein_id","pfam_nbd_hmm","dom_i_evalue","ali_start","ali_end"])
    calls = [] if empty else [
        dict(protein_id="p1",gene_id="g1",transcript_id="t1",scaffold="c",start="1",end="100",strand="+",tier="TIER_1A_HIGH_CONFIDENCE",flags="PASS",nbd_confidence="HIGH",is_fusion_model=False,fused_into=""),
        dict(protein_id="p2",gene_id="g1",transcript_id="t2",scaffold="c",start="1",end="100",strand="+",tier="TIER_2C_LOW_PRIORITY_FRAGMENT",flags="NBD_ONLY",nbd_confidence="LOW",is_fusion_model=False,fused_into="FUSION_a"),
        dict(protein_id="FUSION_a",gene_id="",transcript_id="",scaffold="c",start="1",end="200",strand="+",tier="FUSION_RESCUE",flags="REVIEW_REQUIRED",nbd_confidence="HIGH",is_fusion_model=True,fused_into=""),
    ]
    columns = list(calls[0]) if calls else ["protein_id","tier","flags","nbd_confidence","is_fusion_model"]
    table(root, "final_results/nlr_final_report.tsv", calls, columns)
    table(root, "final_results/nlr_strict_candidates.tsv", calls[::2], columns)
    table(root, "results/stage3/tables/tiered_candidates.tsv", calls[:2], columns)
    table(root, "final_results/qc_summary.tsv", [dict(metric="n_total_candidates",value=0 if empty else 2),dict(metric="n_strict_candidates_combined",value=0 if empty else 2)], ["metric","value"])
    for track in ("priority", "comprehensive"):
        tail = "/comprehensive" if track == "comprehensive" else ""
        five, six = "results/stage5" + tail, "results/stage6" + tail
        for name in (f"queries/miniprot_{track}_ids.txt", "queries/exonerate_refine_ids.txt"):
            text(root, five + "/" + name, "" if empty else "p1\np2\n")
        for name in ("queries/exonerate_queries.ids", "queries/exonerate_run.ids"):
            text(root, six + "/" + name, "" if empty else "p1\n")
        table(root, five + "/tables/miniprot_features.tsv", [] if empty else [dict(query="p1",feature="mRNA"),dict(query="p2",feature="CDS")],["query","feature"])
        table(root, six + "/tables/exonerate_summary.tsv", [] if empty else [dict(query="p1")],["query"])
    return seq


def snapshot(root):
    return {str(p.relative_to(root)):(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns)
            for p in root.rglob("*") if p.is_file()}


def test_original_scores_flags_coordinates_and_multiple_segments(tmp_path):
    root = tmp_path / "run"
    seq = fixture(root)
    # Conflicting architecture placeholders must never become measurements.
    text(root, "results/stage2/all_domain_hits.tsv", "protein_id\tdomain\tdom_i_evalue\tstart\tend\np1\tNBD\t0\t1\t40\n")
    before = snapshot(root)
    candidates, review, hits, segments, summary = build_evidence(root, {"EVAL_NBD": "1e-99", "EXON_MAXN": 1})
    indexed = {r["protein_id"]:r for r in candidates}
    assert indexed["p1"]["best_domain_i_evalue"] == 1e-9
    assert indexed["p1"]["best_domain_bitscore"] == 38.5
    assert indexed["p1"]["best_model_coverage"] == 9 / 20
    assert indexed["p1"]["nbd_confidence"] == "HIGH"
    assert all(h["passes_strict_recorded"] for h in hits if h["source"] == "nbd_library" and h["protein_id"] == "p1")
    pfam = next(h for h in hits if h["source"] == "pfam_nbd")
    assert pfam["passes_strict_recorded"] == ""  # no cross-source join
    assert pfam["accepted_pfam_recorded"] is True
    exported = dict(segments)
    first = next(h for h in hits if h["source"] == "nbd_library" and h["protein_id"] == "p1")
    assert exported[first["hit_id"]] == seq[2:12]
    assert all(len(s) == 10 for _,s in segments)  # never the 3–32 bounding interval
    assert len(exported) == 4
    assert len(hits) == 5
    assert summary["counts"]["accepted_reported_hits"] == 4
    assert summary["counts"]["unassessed_reported_hits"] == 0
    assert all(h["model_accession"] == "PF00001.2" for h in hits)
    assert summary["counts"]["original_final_candidates"] == 2
    assert summary["counts"]["final_report_rows"] == 3
    assert summary["annotation_units"]["distinct_nonempty_gene_ids"] == 1
    assert summary["annotation_units"]["gene_ids_with_multiple_candidate_rows"] == 1
    assert all(r["matches"] for r in summary["qc_reconciliation"])
    assert snapshot(root) == before


def test_review_exclusion_and_rescue_accounting_follow_recorded_membership(tmp_path):
    fixture(tmp_path)
    rows, review, _, _, summary = build_evidence(tmp_path, {"EXON_MAXN":1})
    by_id = {r["protein_id"]:r for r in rows}
    assert "NOT_IN_CANDIDATE_UNION" in by_id["excluded"]["decision_notes"]
    assert "NOT_IN_STRICT_SET" in by_id["p2"]["decision_notes"]
    assert "MEMBER_OF_FUSION_HYPOTHESIS" in by_id["p2"]["decision_notes"]
    assert by_id["FUSION_a"]["record_type"] == "FUSION_HYPOTHESIS"
    assert "ABSENT_FROM_SCANNED_PROTEOME" not in by_id["FUSION_a"]["decision_notes"]
    assert {r["protein_id"] for r in review} == {"p2","excluded","FUSION_a"}
    assert by_id["p2"]["priority_exonerate_not_selected"] is True
    assert by_id["p2"]["priority_exonerate_run_recorded"] is False
    assert by_id["p2"]["priority_miniprot_mrna_observed"] is False
    assert summary["rescue"]["priority"]["exonerate_not_selected"] == 1
    shared = next(r for r in summary["channel_overlaps"] if (r["left"],r["right"]) == ("nbd_relaxed","pfam_enriched"))
    assert shared["shared_ids"] == 2


@pytest.mark.parametrize("line,expected", [
    (raw("p1",end=99),"INVALID_QUERY_COORDINATES"),
    (raw("p1",start=0),"INVALID_QUERY_COORDINATES"),
    (raw("p1",start=15,end=3),"INVALID_QUERY_COORDINATES"),
    (raw("p1",qlen=41),"QUERY_LENGTH_MISMATCH"),
    (raw("p1",he=99),"INVALID_MODEL_COORDINATES"),
    (raw("not_present"),"MISSING_SCANNED_PROTEIN"),
])
def test_invalid_spans_are_recorded_and_never_silently_clamped(tmp_path,line,expected):
    fixture(tmp_path)
    text(tmp_path,"results/stage1/nbd_whole.domtblout",line)
    fields = line.split()
    table(tmp_path, "results/stage1/nbd_hits_detailed.tsv", [dict(
        protein_id=fields[3], nbd_hmm=fields[0], dom_i_evalue=fields[12], ali_start=fields[17],
        ali_end=fields[18], passes_strict=True, passes_relaxed=True)],
        ["protein_id", "nbd_hmm", "dom_i_evalue", "ali_start", "ali_end", "passes_strict", "passes_relaxed"])
    _,_,hits,segments,_ = build_evidence(tmp_path)
    custom = next(h for h in hits if h["source"] == "nbd_library")
    assert custom["segment_status"] == expected
    assert custom["hit_id"] not in dict(segments)
    assert custom["segment_sha256"] == ""


def test_malformed_and_nonfinite_raw_lines_have_explicit_audit(tmp_path):
    fixture(tmp_path)
    text(tmp_path,"results/stage1/nbd_whole.domtblout","bad line\n" + raw("p1",ev="nan") + raw("p1",ev="-1"))
    _,_,hits,_,summary = build_evidence(tmp_path)
    assert len(hits) == 1  # the independent Pfam source is preserved
    assert len(summary["warnings"]) == 3
    json.dumps(summary,allow_nan=False)


def test_missing_sources_are_unknown_not_negative(tmp_path):
    fixture(tmp_path)
    (tmp_path/"results/stage1/nbd_whole.domtblout").unlink()
    (tmp_path/"results/stage1/nbd_candidate_ids_strict.txt").unlink()
    (tmp_path/"results/stage6/queries/exonerate_queries.ids").unlink()
    rows,_,_,_,summary = build_evidence(tmp_path)
    p1 = next(r for r in rows if r["protein_id"] == "p1")
    assert p1["nbd_library_reported_hits"] == p1["nbd_strict"] == ""
    assert p1["priority_exonerate_not_selected"] == ""
    assert summary["discovery_channels"]["nbd_strict"] is None
    assert summary["availability"] == "PARTIAL"


def test_rejected_raw_hits_remain_auditable_without_sequence(tmp_path):
    fixture(tmp_path)
    _, _, hits, segments, summary = build_evidence(tmp_path)
    rejected = next(h for h in hits if h["protein_id"] == "excluded")
    assert rejected["accepted_by_recorded_filters"] is False
    assert rejected["segment_status"] == "NOT_ACCEPTED_BY_RECORDED_FILTERS"
    assert rejected["segment_sha256"] == ""
    assert rejected["hit_id"] not in dict(segments)
    assert summary["counts"]["reported_hits"] == 5
    assert summary["counts"]["exported_segments"] == 4


@pytest.mark.parametrize("filename,source", [
    ("nbd_hits_detailed.tsv", "nbd_library"), ("pfam_nbd_hits.tsv", "pfam_nbd"),
])
def test_missing_acceptance_table_is_unknown_and_never_exports(tmp_path, filename, source):
    fixture(tmp_path)
    (tmp_path / "results/stage1" / filename).unlink()
    _, _, hits, segments, summary = build_evidence(tmp_path)
    affected = [h for h in hits if h["source"] == source]
    assert affected
    assert all(h["accepted_by_recorded_filters"] == "" and h["segment_status"] == "NOT_ASSESSED" for h in affected)
    assert not {h["hit_id"] for h in affected} & dict(segments).keys()
    assert summary["counts"]["unassessed_reported_hits"] == len(affected)


def test_uncompleted_optional_scan_is_distinct_from_no_hits(tmp_path):
    fixture(tmp_path, empty=True)
    text(tmp_path, "results/stage1/run_info.txt", "pfam_nbd_scan: 1\npfam_nbd_scan_completed: False\n")
    *_, summary = build_evidence(tmp_path)
    assert summary["scan_status"]["pfam_nbd_requested"] is True
    assert summary["scan_status"]["pfam_nbd_completed"] is False
    assert summary["discovery_channels"]["pfam_nbd"] == 0


def test_empty_candidates_write_header_only_tables_and_empty_valid_fasta(tmp_path):
    root = tmp_path/"run"
    fixture(root,empty=True)
    destinations = [tmp_path/"first",tmp_path/"second"]
    export_evidence(root,destinations)
    first,second = destinations
    assert snapshot(root)  # populated input with no selected candidates
    assert (first/"nbd_aligned_segments.faa").read_bytes() == b""
    for name in ("candidate_evidence.tsv","review_and_exclusions.tsv","nbd_hit_evidence.tsv"):
        assert len((first/name).read_text().splitlines()) == 1
    summary = json.loads((first/"evidence_summary.json").read_text())
    assert summary["counts"]["evidence_ids"] == 0
    assert summary["counts"]["reported_hits"] == 0
    assert summary["availability"] == "COMPLETE"
    assert {p.name:p.read_bytes() for p in first.iterdir()} == {p.name:p.read_bytes() for p in second.iterdir()}


def test_exports_repeat_exactly_and_fasta_mapping_hashes_match(tmp_path):
    root = tmp_path/"run"
    fixture(root)
    target = tmp_path/"exports"
    export_evidence(root,[target])
    before = {p.name:p.read_bytes() for p in target.iterdir()}
    export_evidence(root,[target])
    assert before == {p.name:p.read_bytes() for p in target.iterdir()}
    with (target/"nbd_aligned_segments.faa").open() as handle:
        sequences = dict(fasta_iter(handle))
    with (target/"nbd_hit_evidence.tsv").open() as handle:
        hits = list(csv.DictReader(handle,delimiter="\t"))
    assert len(sequences) == sum(hit["segment_status"] == "EXPORTED" for hit in hits)
    for hit in hits:
        if hit["segment_status"] == "EXPORTED":
            assert hashlib.sha256(sequences[hit["hit_id"]].encode()).hexdigest() == hit["segment_sha256"]
        else:
            assert hit["hit_id"] not in sequences and hit["segment_sha256"] == ""


def test_zero_candidate_finalization_records_all_empty_evidence_artifacts(tmp_path):
    from funlr.stages.stage7_finalize import run
    fixture(tmp_path, empty=True)
    config = SimpleNamespace(data={"reporting": {"plots": False, "integration": False}},
                             get=lambda section, key, default=None: False,
                             scientific_settings=lambda: {})
    paths = SimpleNamespace(output_dir=tmp_path, final_dir=tmp_path / "final_results",
                            stage_dir=lambda stage: tmp_path / "results" / f"stage{stage}")
    outputs = run(SimpleNamespace(dry_run=False, config=config, paths=paths))
    expected = {"candidate_evidence.tsv", "review_and_exclusions.tsv", "nbd_hit_evidence.tsv",
                "nbd_aligned_segments.faa", "evidence_summary.json"}
    assert {p.name for p in paths.final_dir.joinpath("evidence").iterdir()} == expected
    with (paths.stage_dir(7) / "stage7_manifest.tsv").open() as handle:
        labels = {r["label"] for r in csv.DictReader(handle, delimiter="\t")}
    assert {"evidence/" + name for name in expected} <= labels
    assert all(paths.stage_dir(7).joinpath("evidence", name) in outputs for name in expected)
    assert paths.final_dir.joinpath("evidence/nbd_aligned_segments.faa").read_bytes() == b""
    assert len(paths.final_dir.joinpath("evidence/candidate_evidence.tsv").read_text().splitlines()) == 1
