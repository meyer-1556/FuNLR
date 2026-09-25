"""Only recorded evidence enters verified HTML; resource labels stay literal."""
import hashlib
import json

import pytest

from funlr.report import inspect_run, write_html_report, _elapsed
from test_report import completed_run, snapshot


def add_record(root, relative, content, stage):
    path = root / relative
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(content,sort_keys=True) + "\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    state = json.loads((root/"state.json").read_text())
    state["stages"][str(stage)]["outputs"][relative] = digest
    (root/"state.json").write_text(json.dumps(state))
    if relative.startswith("final_results/"):
        manifest = json.loads((root/"manifest.json").read_text())
        manifest["final_outputs"][relative] = digest
        (root/"manifest.json").write_text(json.dumps(manifest))
    return path


def test_recorded_summaries_resource_scope_and_html_are_read_only(tmp_path):
    root = completed_run(tmp_path/"run")
    add_record(root,"results/stage0/input_qc.json",{
        "schema_version":1,"genome":{"n_sequences":3,"n50_bp":1000},
        "proteins":{"n_sequences":4,"n_starting_with_m":2},"gff3":{"n_gene_features":5}},0)
    add_record(root,"final_results/evidence/evidence_summary.json",{
        "schema_version":1,"discovery_channels":{"nbd_strict":7,'<script>evil</script>':2},
        "rescue":{"priority":{"miniprot_requested":8,"exonerate_requested":5,
                                 "exonerate_selected":2,"exonerate_not_selected":3}}},7)
    before = snapshot(root)
    inspected = inspect_run(root,verify=True)
    assert inspected["integrity"]["status"] == "VERIFIED"
    assert inspected["input_qc"]["genome"]["n50_bp"] == 1000
    assert inspected["evidence_summary"]["discovery_channels"]["nbd_strict"] == 7
    resources = inspected["resource_observations"]
    assert resources["latest_invocation_elapsed_seconds"] == 60
    assert resources["cpu_seconds"] is resources["peak_memory_bytes"] is None
    assert resources["recorded_output_paths"] == 12
    state = json.loads((root/"state.json").read_text())
    unique_paths = {p for stage in state["stages"].values() for p in stage["outputs"]}
    assert resources["observed_recorded_output_bytes"] == sum((root/p).stat().st_size for p in unique_paths)
    assert all(r["elapsed_seconds"] == 60 for r in inspected["stage_progress"]["stages"])
    output = tmp_path/"report.html"
    write_html_report(root,output)
    html = output.read_text()
    assert "Input observations" in html and "Discovery evidence and rescue accounting" in html
    assert "&lt;script&gt;evil&lt;/script&gt;" in html and "<script>evil</script>" not in html
    assert "Peak memory was not measured" in html
    assert snapshot(root) == before


def test_unrecorded_summary_cannot_appear_as_verified_evidence(tmp_path):
    root = completed_run(tmp_path/"run")
    path = root/"final_results/evidence/evidence_summary.json"
    path.parent.mkdir()
    path.write_text('{"schema_version":1,"discovery_channels":{"invented":999}}')
    result = inspect_run(root,verify=True)
    assert result["integrity"]["status"] == "VERIFIED"
    assert result["evidence_summary"] is None


def test_changed_recorded_summary_prevents_html(tmp_path):
    root = completed_run(tmp_path/"run")
    path = add_record(root,"final_results/evidence/evidence_summary.json",{"schema_version":1},7)
    path.write_text('{"schema_version":1,"discovery_channels":{"changed":999}}')
    result = inspect_run(root,verify=True)
    assert result["integrity"]["status"] == "FAILED"
    assert result["integrity"]["changed_count"] == 1
    with pytest.raises(ValueError,match="verified outputs"):
        write_html_report(root,tmp_path/"invalid.html")


@pytest.mark.parametrize("started,finished", [
    (None,None),("invalid","invalid"),("2026-01-01T00:00:01Z","2026-01-01T00:00:00Z"),
    ("2026-01-01T00:00:00Z","2026-01-01T00:01:00"),
])
def test_unusable_timestamps_are_not_measured_zero(started,finished):
    assert _elapsed(started,finished) is None
