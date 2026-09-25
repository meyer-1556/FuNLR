"""Read-only status and local HTML reports over recorded run metadata."""
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path

import pytest

from funlr import __version__
from funlr.report import inspect_run, format_status, write_html_report


def write_json(path, data):
    path.write_text(json.dumps(data, sort_keys=True) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(root):
    return {str(p.relative_to(root)): (digest(p), p.stat().st_mtime_ns)
            for p in root.rglob("*") if p.is_file()}


def completed_run(root, hostile=False):
    """A minimal saved run record, not a replacement for pipeline integration."""
    root.mkdir(parents=True)
    (root / "final_results").mkdir()
    columns = ["protein_id", "tier", "flags", "is_fusion_model", "rescue_priority",
               "nbd_confidence", "fusion_members", "fusion_evidence_rule"]
    hostile_id = '<script>alert("candidate")</script>' if hostile else "p2"
    rows = [
        ["p1", "TIER_1A_HIGH_CONFIDENCE", "PASS", "False", "0", "HIGH", "", ""],
        [hostile_id, "TIER_2A_HIGH_PRIORITY_RESCUE", "SHORT_PROTEIN;HARD_END", "False", "100", "HIGH", "", ""],
        ["fusion1", "FUSION_RESCUE", "FUSION_SPANNING_ALIGNMENT;REVIEW_REQUIRED", "True", "100", "HIGH", "p1;p2", "member-spanning alignment"],
    ]
    for name, selected in (("nlr_final_report.tsv", rows), ("nlr_strict_candidates.tsv", rows[:1])):
        with (root / "final_results" / name).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(selected)
    final = {str(p.relative_to(root)): digest(p) for p in (root / "final_results").iterdir()}
    stages = {}
    for number in range(8):
        path = root / f"results/stage{number}/result.txt"
        path.parent.mkdir(parents=True)
        path.write_text(f"Stage {number}\n")
        outputs = {str(path.relative_to(root)): digest(path)}
        if number == 7:
            outputs.update(final)
        stages[str(number)] = {"status": "complete", "outputs": outputs,
                               "started": "2026-01-01T00:00:00Z", "finished": "2026-01-01T00:01:00Z"}
    config = {"NLR_PROFILE": "ILLUMINA", "MIN_NBD_ALI_LEN": 180, "REPORT_PLOT_BACKEND": "matplotlib"}
    manifest = {
        "status": "complete", "signature": "recorded-signature", "version": "0.2.0a2",
        "sample": '</title><script>alert("sample")</script>' if hostile else "sample1",
        "species_id": "Example_species", "assembly_version": "assembly1",
        "started": "2026-01-02T12:00:00Z", "finished": "2026-01-02T12:01:00Z",
        "config": config, "final_outputs": final, "threads": 4,
        "runtime": {"python": "3.12.14"},
        "tools": {"hmmscan": {"path": "/recorded/tool", "sha256": "a" * 64,
                               "version_output": "HMMER 3.4"}},
        "inputs": {"PFAM_DB": {"path": "/recorded/Pfam-A.hmm", "sha256": "b" * 64}},
        "warnings": ['<img src="https://invalid.example/x" onerror="alert(1)">'] if hostile else [],
    }
    write_json(root / "manifest.json", manifest)
    write_json(root / "state.json", {"signature": "recorded-signature", "stages": stages})
    write_json(root / "resolved_config.json", config)
    return root


def test_recorded_status_and_verify_are_read_only(tmp_path):
    root = completed_run(tmp_path / "run")
    before = snapshot(root)
    recorded = inspect_run(root)
    assert recorded["status"] == "RECORDED_COMPLETE"
    assert recorded["integrity"]["status"] == "NOT_CHECKED"
    assert recorded["reporter_version"] == __version__
    assert recorded["provenance"]["version"] == "0.2.0a2"
    assert recorded["provenance"]["resolved_settings"]["MIN_NBD_ALI_LEN"] == 180
    checked = inspect_run(root, verify=True)
    assert checked["status"] == "RECORDED_COMPLETE"
    assert checked["integrity"]["status"] == "VERIFIED"
    assert checked["integrity"]["expected_files"] == 10
    assert checked["integrity"]["matched_files"] == 10  # final files are not double-counted
    assert checked["counts"] == {
        "n_report_rows": 3, "n_candidates": 2, "n_unique_candidate_ids": 2,
        "n_unique_report_ids": 3, "n_strict_candidates": 1, "n_fusion_models": 1,
        "n_fusion_review_required": 1, "n_flagged_candidates": 2, "n_review_candidates": 2,
    }
    assert checked["flag_counts"]["REVIEW_REQUIRED"] == 1
    assert checked["tier_counts"]["TIER_1A_HIGH_CONFIDENCE"] == 1
    assert checked["stage_progress"]["complete"] == 8
    assert "process liveness not checked" in format_status(checked)
    assert snapshot(root) == before


def test_changed_and_missing_outputs_counted_once(tmp_path):
    root = completed_run(tmp_path / "run")
    (root / "results/stage2/result.txt").write_text("changed")
    (root / "final_results/nlr_strict_candidates.tsv").unlink()
    checked = inspect_run(root, verify=True)
    assert checked["status"] == "FAILED"
    assert checked["integrity"]["status"] == "FAILED"
    assert checked["integrity"]["changed_count"] == 1
    assert checked["integrity"]["missing_count"] == 1
    assert checked["integrity"]["matched_files"] == 8
    with pytest.raises(ValueError, match="verified outputs"):
        write_html_report(root, tmp_path / "failed.html")
    assert not (tmp_path / "failed.html").exists()


@pytest.mark.parametrize("final_outputs", [None, 7, [], "bad", {}])
def test_malformed_checksum_mapping_returns_failed_inspection(tmp_path, final_outputs):
    root = completed_run(tmp_path / "run")
    manifest = json.loads((root / "manifest.json").read_text())
    manifest["final_outputs"] = final_outputs
    write_json(root / "manifest.json", manifest)
    result = inspect_run(root, verify=True)
    assert result["status"] == "FAILED"
    assert result["integrity"]["invalid_record_count"] >= 1


@pytest.mark.parametrize("relative", ["../outside.txt", "/outside.txt", "results/../../outside.txt", "results\\outside.txt"])
def test_recorded_output_paths_cannot_escape_run(tmp_path, relative):
    root = completed_run(tmp_path / "run")
    state = json.loads((root / "state.json").read_text())
    state["stages"]["2"]["outputs"][relative] = "0" * 64
    write_json(root / "state.json", state)
    result = inspect_run(root, verify=True)
    assert result["integrity"]["status"] == "FAILED"
    assert result["integrity"]["unsafe_count"] == 1


def test_symlinked_outputs_and_metadata_are_not_read(tmp_path):
    root = completed_run(tmp_path / "run")
    output = root / "results/stage2/result.txt"
    outside = tmp_path / "outside.txt"
    output.rename(outside)
    output.symlink_to(outside)
    result = inspect_run(root, verify=True)
    assert result["integrity"]["unsafe_count"] == 1
    (root / "manifest.json").unlink()
    (root / "manifest.json").symlink_to(outside)
    result = inspect_run(root)
    assert result["status"] != "RECORDED_COMPLETE"
    assert any("Symlink" in error for error in result["errors"])


def test_header_only_final_tables_are_valid_empty_results(tmp_path):
    root = completed_run(tmp_path / "run")
    manifest = json.loads((root / "manifest.json").read_text())
    state = json.loads((root / "state.json").read_text())
    for relative in manifest["final_outputs"]:
        path = root / relative
        path.write_text(path.read_text().splitlines()[0] + "\n")
        manifest["final_outputs"][relative] = digest(path)
        state["stages"]["7"]["outputs"][relative] = digest(path)
    write_json(root / "manifest.json", manifest)
    write_json(root / "state.json", state)
    result = inspect_run(root, verify=True)
    assert result["integrity"]["status"] == "VERIFIED"
    assert result["counts"]["n_candidates"] == 0
    assert result["counts"]["n_strict_candidates"] == 0
    assert result["tier_counts"] == {}


def test_unknown_incomplete_failed_and_interrupted_status(tmp_path):
    root = tmp_path / "absent"
    assert inspect_run(root)["status"] == "UNKNOWN"
    assert not root.exists()
    root.mkdir()
    assert inspect_run(root)["status"] == "UNKNOWN"
    write_json(root / "manifest.json", {"status": "running"})
    result = inspect_run(root)
    assert result["status"] == "INCOMPLETE"
    assert "interrupted" in result["status_detail"]
    write_json(root / ".run.lock", {"pid": 999999999, "started": "recorded"})
    result = inspect_run(root)
    assert result["status"] == "RUNNING"
    assert result["lock"]["process_liveness"] == "NOT_CHECKED"
    write_json(root / "manifest.json", {"status": "failed", "error": "exonerate exited 1"})
    assert inspect_run(root)["status"] == "FAILED"
    (root / ".run.lock").unlink()
    (root / "manifest.json").write_text("{bad json")
    assert inspect_run(root)["errors"]


def test_conflicting_stage_and_final_checksums_fail(tmp_path):
    root = completed_run(tmp_path / "run")
    state = json.loads((root / "state.json").read_text())
    state["stages"]["7"]["outputs"]["final_results/nlr_final_report.tsv"] = "c" * 64
    write_json(root / "state.json", state)
    result = inspect_run(root, verify=True)
    assert result["integrity"]["invalid_record_count"] == 1
    assert result["status"] == "FAILED"


def test_html_is_escaped_self_contained_deterministic_and_preserves_run(tmp_path):
    root = completed_run(tmp_path / "run", hostile=True)
    before = snapshot(root)
    output, duplicate = tmp_path / "first.html", tmp_path / "second.html"
    result = write_html_report(root, output)
    write_html_report(root, duplicate)
    text = output.read_text()
    assert result["integrity"]["status"] == "VERIFIED"
    assert result["reporter_version"] == __version__
    assert "&lt;script&gt;" in text and "&lt;img" in text
    assert "<script" not in text and "<img" not in text and "<link" not in text
    assert "Content-Security-Policy" in text
    assert "default-src 'none'" in text
    assert "Original inputs and executable files were not rehashed" in text
    assert "latest invocation started" in text
    assert "Created by FuNLR reporter" in text and "recorded FuNLR 0.2.0a2 run" in text
    assert output.read_bytes() == duplicate.read_bytes()
    assert snapshot(root) == before
    with pytest.raises(FileExistsError):
        write_html_report(root, output)
    assert snapshot(root) == before


def test_html_review_table_preserves_each_source_row(tmp_path):
    """Count assertions alone miss deferred generators repeating the last row."""
    class Tables(HTMLParser):
        def __init__(self):
            super().__init__()
            self.tables = []
            self.table = self.row = self.cell = None

        def handle_starttag(self, tag, attrs):
            if tag == "table":
                self.table = []
            elif tag == "tr":
                self.row = []
            elif tag in {"th", "td"}:
                self.cell = []

        def handle_data(self, data):
            if self.cell is not None:
                self.cell.append(data)

        def handle_endtag(self, tag):
            if tag in {"th", "td"}:
                self.row.append("".join(self.cell))
                self.cell = None
            elif tag == "tr":
                self.table.append(self.row)
                self.row = None
            elif tag == "table":
                self.tables.append(self.table)
                self.table = None

    root = completed_run(tmp_path / "run", hostile=True)
    output = tmp_path / "review.html"
    write_html_report(root, output)
    parser = Tables()
    parser.feed(output.read_text())
    review = next(table for table in parser.tables if table[0][0] == "Protein ID")
    assert review[1:] == [
        ['<script>alert("candidate")</script>', "TIER_2A_HIGH_PRIORITY_RESCUE", "100",
         "SHORT_PROTEIN;HARD_END", "HIGH", ""],
        ["fusion1", "FUSION_RESCUE", "100", "FUSION_SPANNING_ALIGNMENT;REVIEW_REQUIRED",
         "HIGH", "p1;p2"],
    ]


@pytest.mark.parametrize("directory", ["results", "final_results", "work", "logs", "other"])
def test_report_rejects_pipeline_subdirectories(tmp_path, directory):
    root = completed_run(tmp_path / "run")
    (root / directory).mkdir(exist_ok=True)
    with pytest.raises(ValueError):
        write_html_report(root, root / directory / "report.html")
    assert not (root / directory / "report.html").exists()


def test_report_allows_new_html_directly_in_run_root(tmp_path):
    root = completed_run(tmp_path / "run")
    before = snapshot(root)
    result = write_html_report(root, root / "report.html")
    assert Path(result["path"]).is_file()
    assert {key: value for key, value in snapshot(root).items() if key != "report.html"} == before


def test_report_rejects_symlink_and_missing_parent(tmp_path):
    root = completed_run(tmp_path / "run")
    target = tmp_path / "linked.html"
    target.symlink_to(tmp_path / "absent.html")
    with pytest.raises(ValueError, match="symlink"):
        write_html_report(root, target)
    directory = tmp_path / "linked-directory"
    directory.symlink_to(root, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        write_html_report(root, directory / "output.html")
    with pytest.raises(ValueError, match="parent directory"):
        write_html_report(root, tmp_path / "not-created" / "report.html")
    assert not (tmp_path / "not-created").exists()


def make_batch(tmp_path):
    root = tmp_path / "batch"
    completed_run(root / "samples/A")
    manifest = {"schema_version": 1, "identity_sha256": "b" * 64,
                "samples": [{"sample_id": "A", "species_id": "Species_A",
                             "output_dir": "samples/A", "metadata": {"ORIGIN": "sample metadata"}}]}
    state = {"schema_version": 1, "identity_sha256": "b" * 64, "status": "COMPLETE",
             "samples": [{"sample_id": "A", "output_dir": "samples/A", "status": "COMPLETE"}]}
    write_json(root / "batch_manifest.json", manifest)
    write_json(root / "batch_state.json", state)
    return root


def test_batch_inspects_live_children_and_preserves_metadata(tmp_path):
    root = make_batch(tmp_path)
    before = snapshot(root)
    result = inspect_run(root, verify=True)
    assert result["kind"] == "batch"
    assert result["status"] == "RECORDED_COMPLETE"
    assert result["integrity"]["status"] == "VERIFIED"
    assert result["integrity"]["matched_files"] == 10
    assert result["samples"][0]["metadata"] == {"ORIGIN": "sample metadata"}
    assert snapshot(root) == before
    (root / "samples/A/results/stage3/result.txt").write_text("tampered")
    result = inspect_run(root, verify=True)
    assert result["status"] == "FAILED"
    assert result["samples"][0]["status"] == "FAILED"
    assert result["integrity"]["changed_count"] == 1
    assert "A: FAILED" in format_status(result)
    with pytest.raises(ValueError, match="one sample run"):
        write_html_report(root, tmp_path / "batch.html")


@pytest.mark.parametrize("change", ["identity", "schema", "missing_sample", "duplicate", "case_duplicate", "wrong_output"])
def test_batch_state_cannot_claim_completion_when_inconsistent(tmp_path, change):
    root = make_batch(tmp_path)
    state = json.loads((root / "batch_state.json").read_text())
    if change == "identity":
        state["identity_sha256"] = "c" * 64
    elif change == "schema":
        state["schema_version"] = 100
    elif change == "missing_sample":
        state["samples"] = []
    elif change in {"duplicate", "case_duplicate"}:
        state["samples"].append({"sample_id": "A" if change == "duplicate" else "a",
                                 "output_dir": "samples/A", "status": "COMPLETE"})
    elif change == "wrong_output":
        state["samples"][0]["output_dir"] = "samples/other"
    write_json(root / "batch_state.json", state)
    result = inspect_run(root, verify=True)
    assert result["status"] == "INCOMPLETE"
    assert result["errors"]
    assert result["integrity"]["status"] == "FAILED"


def test_batch_manifest_path_escape_is_rejected(tmp_path):
    root = make_batch(tmp_path)
    manifest = json.loads((root / "batch_manifest.json").read_text())
    manifest["samples"][0]["output_dir"] = "../outside"
    write_json(root / "batch_manifest.json", manifest)
    result = inspect_run(root, verify=True)
    assert result["status"] != "RECORDED_COMPLETE"
    assert result["samples"] == []
    assert any("within the run" in item for item in result["errors"])
