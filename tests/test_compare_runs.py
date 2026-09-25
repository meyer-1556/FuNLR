"""The acceptance comparator must reveal differences rather than normalize them."""
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_runs.py"
spec = importlib.util.spec_from_file_location("funlr_compare_runs", SCRIPT)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)

TABLE = "final_results/nlr_final_report.tsv"
FASTA = "results/stage4/all_tiered_candidates.faa"
BED = "final_results/nlr_candidates.bed"
IDS = "results/stage1/union_candidate_ids.txt"
COUNTS = "results/stage5/summaries/query_counts.tsv"


def write(root, relative, content):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content.encode("utf-8") if isinstance(content, str) else content)
    return path


@pytest.fixture
def runs(tmp_path):
    before, after = tmp_path / "before", tmp_path / "after"
    for root in (before, after):
        root.mkdir()
        write(root, TABLE, "protein_id\ttier\tflags\tstart\nq1\tTIER_1A_HIGH_CONFIDENCE\tPASS\t100\nq2\tTIER_3_ARCHITECTURAL_VARIANT\tSHORT_PROTEIN\t200\n")
        write(root, FASTA, ">q1\nMAAA\n>q2\nMCCC\n")
        write(root, BED, "scaf1\t99\t300\tq1\t0\t+\n")
        write(root, IDS, "q1\nq2\n")
    return before, after


def compare(runs, required=None):
    return comparison.compare_runs(*runs, required_files={TABLE, FASTA, BED, IDS} if required is None else required)


def test_matching_selected_artifacts(runs):
    report = compare(runs)
    assert report["status"] == "match"
    assert report["matching_files"] == report["checked_files"] == 4
    assert report["differences"] == []


@pytest.mark.parametrize("replace", [
    lambda text: text.replace("\ttier\tflags", "\tflags\ttier"),
    lambda text: "\n".join([text.splitlines()[0], text.splitlines()[2], text.splitlines()[1]]) + "\n",
    lambda text: text.replace("\t100\n", "\t100.0\n"),
    lambda text: text.replace("\tPASS\t", "\tNaN\t"),
    lambda text: text.replace("\tPASS\t", "\t\t"),
    lambda text: text.replace("\n", "\r\n"),
    lambda text: text.replace("q1\t", '"q1"\t'),
])
def test_schema_values_order_and_format_differences_fail(runs, replace):
    before, after = runs
    write(after, TABLE, replace((before / TABLE).read_text()))
    report = compare(runs)
    assert report["status"] == "different"
    assert report["differences"][0]["path"] == TABLE
    assert report["differences"][0]["status"] == "different"


def test_same_tsv_values_with_new_line_endings_is_explicit_format_difference(runs):
    before, after = runs
    write(after, TABLE, (before / TABLE).read_bytes().replace(b"\n", b"\r\n"))
    assert compare(runs)["differences"][0]["format_only"] is True


@pytest.mark.parametrize("replacement", [
    ">q1\nMA\nAA\n>q2\nMCCC\n",  # same sequences, different wrapping
    ">q2\nMCCC\n>q1\nMAAA\n",  # different record order
    ">q1 description\nMAAA\n>q2\nMCCC\n",  # header text is preserved
    ">q1\nMAAT\n>q2\nMCCC\n",  # residue change
])
def test_fasta_is_not_rewrapped_reordered_or_header_stripped(runs, replacement):
    write(runs[1], FASTA, replacement)
    report = compare(runs)
    assert report["status"] == "different"
    assert report["differences"][0]["path"] == FASTA


def test_bed_coordinate_change_and_id_order_are_compared(runs):
    write(runs[1], BED, "scaf1\t100\t300\tq1\t0\t+\n")
    write(runs[1], IDS, "q2\nq1\n")
    assert {item["path"] for item in compare(runs)["differences"]} == {BED, IDS}


def test_missing_and_unexpected_outputs_fail(runs):
    (runs[1] / FASTA).unlink()
    extra = "results/stage3/tables/TIER_UNEXPECTED.tsv"
    write(runs[1], extra, "protein_id\ttier\nq1\tTIER_UNEXPECTED\n")
    differences = {item["path"]: item for item in compare(runs)["differences"]}
    assert differences[FASTA]["missing_from"] == ["candidate"]
    assert differences[extra]["missing_from"] == ["reference"]


def test_two_empty_directories_cannot_pass(tmp_path):
    before, after = tmp_path / "a", tmp_path / "b"
    before.mkdir()
    after.mkdir()
    report = comparison.compare_runs(before, after)
    assert report["status"] == "different"
    assert report["matching_files"] == 0
    assert len(report["differences"]) == len(comparison.REQUIRED_FILES)
    assert all(item["missing_from"] == ["reference", "candidate"] for item in report["differences"])


def test_zero_candidate_exports_are_valid_but_empty_tsv_is_not(runs):
    for root in runs:
        write(root, FASTA, "")
        write(root, BED, "")
        write(root, IDS, "")
        write(root, TABLE, "protein_id\ttier\tflags\n")
    assert compare(runs)["status"] == "match"
    for root in runs:
        write(root, TABLE, "")
    report = compare(runs)
    assert report["status"] == "different"
    assert report["differences"][0]["status"] == "invalid"


def test_only_declared_path_metadata_is_excluded(runs):
    before, after = runs
    write(before, COUNTS, 'label\tfile\tseqs\ndedup\t"/old/path\tquoted"\t2\n')
    write(after, COUNTS, "label\tfile\tseqs\ndedup\t/new/path\t2\n")
    for root, timestamp in ((before, "yesterday"), (after, "today")):
        write(root, "results/stage3/metadata/run_info.txt", timestamp)
        write(root, "results/stage5/stage5_manifest.tsv", "label\tpath\nqueries\t" + str(root) + "\n")
    report = compare(runs)
    assert report["status"] == "match"
    assert report["exclusions"]["columns"][COUNTS] == ["file"]
    assert "results/stage5/stage5_manifest.tsv" in report["exclusions"]["files"]
    write(after, COUNTS, "label\tfile\tseqs\ndedup\t/new/path\t3\n")
    assert compare(runs)["differences"][0]["path"] == COUNTS


def test_metadata_column_does_not_hide_format_or_header_changes(runs):
    write(runs[0], COUNTS, "label\tfile\tseqs\ndedup\t/old/path\t2\n")
    write(runs[1], COUNTS, "label\tfile\tseqs\r\ndedup\t/new/path\t2\r\n")
    assert compare(runs)["differences"][0]["format_only"] is True
    write(runs[1], COUNTS, "file\tlabel\tseqs\n/new/path\tdedup\t2\n")
    assert compare(runs)["differences"][0]["header_changed"] is True


def test_invalid_tsv_shape_is_not_accepted_even_if_identical(runs):
    for root in runs:
        write(root, TABLE, "protein_id\ttier\tflags\nq1\tTIER_1A_HIGH_CONFIDENCE\n")
    assert compare(runs)["differences"][0]["status"] == "invalid"


def test_cli_json_and_nonzero_difference_exit(runs, tmp_path):
    target = tmp_path / "comparison.json"
    # The selected unit fixture is intentionally incomplete for full-run acceptance.
    assert comparison.main([str(runs[0]), str(runs[1]), "--json", str(target)]) == 1
    assert json.loads(target.read_text())["status"] == "different"


def test_cli_refuses_to_overwrite_scientific_results(runs):
    path = runs[0] / TABLE
    original = path.read_bytes()
    assert comparison.main([str(runs[0]), str(runs[1]), "--json", str(path)]) == 2
    assert path.read_bytes() == original


def test_cli_preserves_run_provenance(runs):
    path = write(runs[0], "manifest.json", '{"status":"complete"}\n')
    original = path.read_bytes()
    assert comparison.main([str(runs[0]), str(runs[1]), "--json", str(path)]) == 2
    assert path.read_bytes() == original


def test_same_run_and_nonexistent_input_are_invalid(runs, tmp_path):
    with pytest.raises(ValueError, match="different run directories"):
        comparison.compare_runs(runs[0], runs[0])
    with pytest.raises(ValueError, match="existing FuNLR run"):
        comparison.compare_runs(runs[0], tmp_path / "missing")


def test_v3_and_legacy_architecture_changes_are_reported_without_conversion(runs):
    name = "results/stage2/architecture_summary.tsv"
    write(runs[0], name, "protein_id\thas_nbd\thas_repeat\torder_ok\nq1\tTrue\tTrue\tTrue\n")
    write(runs[1], name, "protein_id\thas_nbd\thas_sensor\torder_invalid\nq1\tTrue\tTrue\tFalse\n")
    difference = next(item for item in compare(runs)["differences"] if item["path"] == name)
    assert difference["status"] == "different"
    assert difference["header_changed"] is True
    assert difference["format_only"] is False
    assert "results/stage4/all_tiered_candidates.faa" in comparison.required_for(runs[0])
    assert "results/stage4/tiered_candidates.faa" in comparison.required_for(runs[1])


def test_v3_miniprot_query_schema_and_legacy_protein_id_schema_are_valid(runs):
    name = "results/stage5/tables/miniprot_summary.tsv"
    write(runs[0], name, "protein_id\tscaffold\tstart\tend\nq1\tctg\t1\t20\n")
    write(runs[1], name, "query\tscaffold\tstart\tend\nq1\tctg\t1\t20\n")
    difference = next(item for item in compare(runs)["differences"] if item["path"] == name)
    assert difference["status"] == "different"
    assert difference["header_changed"] is True


def test_v3_tier_calls_require_dynamic_table_and_exports_even_when_missing_twice(runs):
    tier = "TIER_1C_CANONICAL_NON_SSFR_SENSOR"
    for root in runs:
        write(root, "results/stage3/tables/tiered_candidates.tsv", f"protein_id\ttier\tflags\trescue_priority\nq1\t{tier}\tPASS\t50\n")
    report = comparison.compare_runs(*runs)
    missing = {item["path"]: item for item in report["differences"] if item["status"] == "missing"}
    for name in (f"results/stage3/tables/{tier}.tsv", f"results/stage4/{tier}.ids", f"results/stage4/{tier}.faa"):
        assert missing[name]["missing_from"] == ["reference", "candidate"]
    assert "results/stage5/comprehensive/tables/miniprot_summary.tsv" in missing
    assert "results/stage6/comprehensive/tables/fusion_evidence.tsv" in missing


def test_v3_path_columns_are_masked_without_hiding_counts_or_status(runs):
    files = {
        "results/stage1/priority_sources_summary.tsv": "source\tcount\tfile\nnbd_strict\t2\t{path}\n",
        "results/stage5/comprehensive/summaries/query_counts.tsv": "label\tfile\tseqs\ndedup\t{path}\t2\n",
        "final_results/integration/export_status.tsv": "export\tstatus\tpath\ngff3\tWRITTEN_EXISTING_MODELS\t{path}\n",
    }
    for root in runs:
        for name, content in files.items():
            write(root, name, content.format(path=root))
    assert compare(runs)["status"] == "match"
    name = "final_results/integration/export_status.tsv"
    write(runs[1], name, files[name].format(path=runs[1]).replace("WRITTEN_EXISTING_MODELS", "SKIPPED"))
    assert compare(runs)["differences"][0]["path"] == name


def test_integration_tags_and_gene_id_lists_are_scientific_outputs(runs):
    for root in runs:
        write(root, "final_results/integration/input.withNLR.gff3", "ctg\tx\tmRNA\t1\t20\t.\t+\t.\tID=q1;NLR_tier=TIER_1A\n")
        write(root, "final_results/integration/nlr_genes_T1.txt", "q1\n")
        write(root, "final_results/domain_plots/R_sessionInfo.txt", str(root))
    assert compare(runs)["status"] == "match"
    write(runs[1], "final_results/integration/nlr_genes_T1.txt", "")
    write(runs[1], "final_results/integration/input.withNLR.gff3", "ctg\tx\tmRNA\t1\t21\t.\t+\t.\tID=q1;NLR_tier=TIER_1A\n")
    assert {item["path"] for item in compare(runs)["differences"]} == {
        "final_results/integration/nlr_genes_T1.txt", "final_results/integration/input.withNLR.gff3"}
