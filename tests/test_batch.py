"""TSV boundary, preflight, and sequential execution contracts.

The mocked runner checks cohort orchestration; existing pipeline tests and the
separate native cohort validation exercise biological stage execution.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import pytest
import yaml

from funlr import batch


def arguments(sheet, output, config=None, *flags):
    parser = argparse.ArgumentParser()
    batch.add_batch_arguments(parser)
    argv = ["--samples", str(sheet), "--outdir", str(output)]
    if config:
        argv += ["--config", str(config)]
    return parser.parse_args([*argv, *flags])


def inputs(folder):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "genome.fa").write_text(">ctg\n" + "ATG" * 20 + "\n")
    (folder / "proteins.fa").write_text(">Q\n" + "M" * 20 + "\n")
    (folder / "genes.gff3").write_text(
        "##gff-version 3\n"
        "ctg\ttest\tgene\t1\t60\t.\t+\t.\tID=G\n"
        "ctg\ttest\tmRNA\t1\t60\t.\t+\t.\tID=Q;Parent=G\n"
        "ctg\ttest\tCDS\t1\t60\t.\t+\t0\tParent=Q\n"
    )
    for name in ("nbd.hmm", "pfam.hmm"):
        (folder / name).write_text("HMMER3/f [TEST INPUT; no native search is performed]\n//\n")


def fixture(tmp_path, rows="a\tfungus one\nb\tfungus two\n"):
    source = tmp_path / "source"
    inputs(source)
    config = source / "shared.yaml"
    config.write_text(yaml.safe_dump({
        "inputs": {"genome": "genome.fa", "proteins": "proteins.fa", "gff3": "genes.gff3"},
        "databases": {"pfam": "pfam.hmm", "nbd_hmms": "nbd.hmm"},
        "asm": {"asm_enable": 0},
        "execution": {"cpu_threads": 2, "ram_gb": 12},
    }))
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("SAMPLE_ID\tSPECIES_ID\n" + rows)
    return sheet, tmp_path / "cohort", config


def fake_runner(calls, fail=None, interrupt=None):
    def run(args, file_inputs, scientific, warnings, counts):
        calls.append(args)
        if args.sample == fail:
            raise RuntimeError("explicit native-tool failure")
        if args.sample == interrupt:
            raise KeyboardInterrupt()
        root = Path(args.outdir)
        final = root / "final_results"
        final.mkdir(parents=True, exist_ok=True)
        report = "protein_id\ttier\tis_fusion_model\nQ\tTIER_1A_HIGH_CONFIDENCE\tFalse\nF\tFUSION_RESCUE\tTrue\n"
        (final / "nlr_final_report.tsv").write_text(report)
        (final / "nlr_strict_candidates.tsv").write_text(report)
        manifest = {"inputs": {key: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for key, path in file_inputs.items()}}
        (root / "manifest.json").write_text(json.dumps(manifest))
    return run


def test_template_and_option_aliases():
    template = batch.sample_template()
    assert template.endswith("\n")
    assert next(csv.reader([template], delimiter="\t")) == list(batch.COLUMNS)
    parser = argparse.ArgumentParser()
    batch.add_batch_arguments(parser)
    args = parser.parse_args(["--input_tsv", "samples.tsv", "--output-dir", "out", "--cpu-threads", "3", "--ram_gb", "16"])
    assert (args.samples, args.outdir, args.threads, args.ram_gb) == ("samples.tsv", "out", 3, 16)
    assert args.keep_going is False
    for alias in ("--stop-on-error", "--stop_on_error"):
        assert parser.parse_args(["--samples", "samples.tsv", "-o", "out", alias]).keep_going is False
        with pytest.raises(SystemExit):
            parser.parse_args(["--samples", "samples.tsv", "-o", "out", alias, "--keep-going"])


@pytest.mark.parametrize("sheet_text,match", [
    ("SAMPLE_ID\tSAMPLE_ID\tSPECIES_ID\na\ta\tsp\n", "Duplicate"),
    (" sample_id \tSAMPLE_ID\tSPECIES_ID\na\tb\tsp\n", "Duplicate"),
    ("SAMPLE_ID\tSPECIES_ID\tFUNLR_PROFILE\na\tsp\tfungi_general\n", "Unsupported"),
    ("SAMPLE_ID\tSPECIES_ID\tREPEATS_PATH\na\tsp\tpath\n", "Unsupported"),
    ("SAMPLE_ID\tSPECIES_ID\tEVAL_NBD\na\tsp\t1e-8\n", "Unsupported"),
    ("SAMPLE_ID\tSPECIES_ID\tPROTEOME_PATH\na\tsp\tpath\n", "Unsupported"),
    ('SAMPLE_ID\t"SPECIES_ID\t"\na\tsp\n', "Control"),
    ("SAMPLE_ID\na\n", "Missing sample-sheet"),
    ("SAMPLE_ID\tSPECIES_ID\na\tsp\textra\n", "cells"),
    ("SAMPLE_ID\tSPECIES_ID\na\tsp\nA\tsp\n", "case-colliding"),
    ("SAMPLE_ID\tSPECIES_ID\n../escape\tsp\n", "Unsafe"),
    ("SAMPLE_ID\tSPECIES_ID\n.\tsp\n", "Unsafe"),
    ("SAMPLE_ID\tSPECIES_ID\n..\tsp\n", "Unsafe"),
    ("SAMPLE_ID\tSPECIES_ID\n/a\tsp\n", "Unsafe"),
    ("SAMPLE_ID\tSPECIES_ID\n" + "a" * 129 + "\tsp\n", "Unsafe"),
    ("SAMPLE_ID\tSPECIES_ID\na\tNA\n", "SPECIES_ID is required"),
    ('SAMPLE_ID\tSPECIES_ID\na\t"two\nlines"\n', "Multiline"),
    ('SAMPLE_ID\tSPECIES_ID\na\t"two\twords"\n', "Control"),
    ("SAMPLE_ID\tSPECIES_ID\na\tbad\x00cell\n", "Control"),
    ("SAMPLE_ID\tSPECIES_ID\tORGANISM_KINGDOM\na\tsp\tPlantae\n", "ORGANISM_KINGDOM"),
    ("SAMPLE_ID\tSPECIES_ID\tORGANISM_KARYOTE\na\tsp\tprokaryote\n", "ORGANISM_KARYOTE"),
    ("SAMPLE_ID\tSPECIES_ID\tNLR_PROFILE\na\tsp\tunknown\n", "NLR_PROFILE"),
    ("SAMPLE_ID\tSPECIES_ID\tDISCOVERY_MODE\na\tsp\tunknown\n", "DISCOVERY_MODE"),
    ("SAMPLE_ID\tSPECIES_ID\tASSEMBLY_PROFILE\tNLR_PROFILE\na\tsp\tillumina\tHIFI\n", "conflicting ASSEMBLY_PROFILE"),
])
def test_sheet_rejection_is_before_writes_or_execution(tmp_path, monkeypatch, sheet_text, match):
    sheet, output, config = fixture(tmp_path)
    sheet.write_text(sheet_text)
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("runner must not be called"))
    with pytest.raises(ValueError, match=match):
        batch.run_batch(arguments(sheet, output, config))
    assert not output.exists()


@pytest.mark.parametrize("flag,status", [("--dry-run", "PLANNED"), ("--validate-only", "VALIDATED")])
def test_no_write_modes_validate_inputs_without_tools(tmp_path, monkeypatch, flag, status):
    sheet, output, config = fixture(tmp_path)
    before = {path: path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("runner must not be called"))
    import funlr.runner
    monkeypatch.setattr(funlr.runner, "tool_inventory", lambda *args: pytest.fail("tools must not run"))
    result = batch.run_batch(arguments(sheet, output, config, flag))
    assert result["status"] == status
    assert result["sample_count"] == 2
    assert result["samples"][0]["input_counts"] == {"contigs": 1, "proteins": 1, "transcripts": 1}
    assert not output.exists()
    assert before == {path: path.read_bytes() for path in tmp_path.rglob("*") if path.is_file()}


def test_later_invalid_sample_prevents_first_sample_from_running(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    sheet.write_text("SAMPLE_ID\tSPECIES_ID\tGENOME_PATH\na\tsp\t\nb\tsp\tmissing.fa\n")
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("no sample may run before whole-sheet validation"))
    with pytest.raises(ValueError, match="Sample b: Missing or empty"):
        batch.run_batch(arguments(sheet, output, config))
    assert not output.exists()


def test_profile_and_path_precedence_preserve_explicit_yaml_values(tmp_path):
    sheet, output, config = fixture(tmp_path)
    settings = yaml.safe_load(config.read_text())
    settings["thresholds"] = {"eval_pfam": "2e-6"}
    settings["execution"].update(resume=True, dry_run=True)
    config.write_text(yaml.safe_dump(settings))
    other = tmp_path / "row-data"
    inputs(other)
    sheet.write_text(
        "SAMPLE_ID\tSPECIES_ID\tPROTEINS_PATH\tGENOME_PATH\tASSEMBLY_PROFILE\tPLOIDY\tORGANISM_KINGDOM\tORGANISM_KARYOTE\n"
        "a\tfungus one\trow-data/proteins.fa\tNone\tHIFI\tunknown\tFunga\tEUKARYOTE\n"
        "b\tfungus two\tNA\tnull\tILLUMINA\t1\t\t\n"
    )
    result = batch.run_batch(arguments(sheet, output, config, "--validate-only", "--threads", "4"))
    a, b = result["samples"]
    assert a["resolved_config"]["inputs"]["genome"] == str(config.parent / "genome.fa")
    assert a["resolved_config"]["inputs"]["proteins"] == str(other / "proteins.fa")
    assert b["resolved_config"]["inputs"]["proteins"] == str(config.parent / "proteins.fa")
    assert a["resolved_config"]["thresholds"]["eval_pfam"] == "2e-6"
    assert a["resolved_config"]["architecture"]["nbd_conf_cov_high"] == 0.5
    assert b["resolved_config"]["architecture"]["nbd_conf_cov_high"] == 0.45
    assert a["resolved_config"]["execution"] == {"output_dir": str(output / "samples/a"), "cpu_threads": 4, "ram_gb": 12, "resume": False, "dry_run": False}
    assert a["metadata"]["PLOIDY"] == "unknown" and b["metadata"]["PLOIDY"] == "1"


def test_normalized_headers_profiles_and_discovery_presets_keep_explicit_thresholds(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    settings = yaml.safe_load(config.read_text())
    settings["thresholds"] = {"eval_pfam": "2e-6"}
    settings["discovery"] = {"discovery_mode": "BROAD"}
    config.write_text(yaml.safe_dump(settings))
    sheet.write_text(
        " sample_id \tspecies_id\t nlr_profile \tassembly_profile\tdiscovery_mode\n"
        "a\tfungus one\t hifi \tHIFI\tstrict\n"
        "b\tfungus two\tillumina\tNA\tbalanced\n"
        "c\tfungus three\tnull\tNone\t\n"
    )
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("no runner in preflight"))
    result = batch.run_batch(arguments(sheet, output, config, "--dry-run"))
    a, b, c = [sample["resolved_config"] for sample in result["samples"]]
    assert a["context"]["profile"] == "HIFI"
    assert a["architecture"]["nbd_conf_cov_high"] == 0.5
    assert b["architecture"]["nbd_conf_cov_high"] == 0.45
    assert [sample["discovery"]["discovery_mode"] for sample in (a, b, c)] == ["STRICT", "BALANCED", "BROAD"]
    assert [sample["discovery"]["pfam_enrich_min_log2fc"] for sample in (a, b, c)] == [1.5, 1.0, 0.5]
    assert all(sample["thresholds"]["eval_pfam"] == "2e-6" for sample in (a, b, c))
    assert not output.exists()
    settings["discovery"]["pfam_enrich_min_log2fc"] = 2.25
    config.write_text(yaml.safe_dump(settings))
    explicit = batch.run_batch(arguments(sheet, output, config, "--validate-only"))
    assert all(sample["resolved_config"]["discovery"]["pfam_enrich_min_log2fc"] == 2.25 for sample in explicit["samples"])


def test_masked_genome_fallback_uses_tsv_relative_paths_and_warns_when_primary_wins(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    row_data = tmp_path / "row-data"
    inputs(row_data)
    sheet.write_text(
        "SAMPLE_ID\tSPECIES_ID\tGENOME_PATH\tMASKED_GENOME_PATH\n"
        "a\tsp\tNA\trow-data/genome.fa\n"
        "b\tsp\trow-data/genome.fa\tunused-missing.fa\n"
        "c\tsp\tnull\tNone\n"
    )
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("no runner in preflight"))
    result = batch.run_batch(arguments(sheet, output, config, "--validate-only"))
    a, b, c = result["samples"]
    assert a["resolved_config"]["inputs"]["genome"] == b["resolved_config"]["inputs"]["genome"] == str(row_data / "genome.fa")
    assert c["resolved_config"]["inputs"]["genome"] == str(config.parent / "genome.fa")
    assert not any("MASKED_GENOME_PATH" in warning for warning in a["warnings"])
    assert "Both GENOME_PATH and MASKED_GENOME_PATH are set; using GENOME_PATH." in b["warnings"]
    assert not output.exists()
    sheet.write_text(sheet.read_text().replace("a\tsp\tNA\trow-data/genome.fa", "a\tsp\tNA\tmissing.fa"))
    with pytest.raises(ValueError, match="Sample a: Missing or empty"):
        batch.run_batch(arguments(sheet, output, config, "--dry-run"))
    assert not output.exists()


def test_new_row_settings_are_frozen_and_rechecked_on_resume(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    sheet.write_text(" sample_id \tspecies_id\tnlr_profile\tdiscovery_mode\tmasked_genome_path\na\tsp\tHIFI\tSTRICT\tsource/genome.fa\n")
    calls = []
    monkeypatch.setattr(batch, "run_sample", fake_runner(calls))
    batch.run_batch(arguments(sheet, output, config))
    assert (calls[0].scientific["NLR_PROFILE"], calls[0].scientific["DISCOVERY_MODE"]) == ("HIFI", "STRICT")
    manifest = (output / "batch_manifest.json").read_bytes()
    batch.run_batch(arguments(sheet, output, config, "--resume"))
    assert calls[-1].resume is True
    assert (output / "batch_manifest.json").read_bytes() == manifest
    before = {path: path.read_bytes() for path in output.rglob("*") if path.is_file()}
    calls.clear()
    sheet.write_text(sheet.read_text().replace("STRICT", "BROAD"))
    with pytest.raises(ValueError, match="Cannot resume"):
        batch.run_batch(arguments(sheet, output, config, "--resume"))
    assert not calls
    assert before == {path: path.read_bytes() for path in output.rglob("*") if path.is_file()}


def test_input_output_ancestor_collision_before_writes(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("runner must not run"))
    # --resume would otherwise require a batch manifest, but input ownership
    # collision is diagnosed before any manifest read or output mutation.
    with pytest.raises(ValueError, match="must not contain each other"):
        batch.run_batch(arguments(sheet, tmp_path, config, "--resume"))
    assert not (tmp_path / "batch_manifest.json").exists()


def test_sequential_success_copies_sheet_and_preserves_frozen_manifest_on_resume(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    calls = []
    monkeypatch.setattr(batch, "run_sample", fake_runner(calls))
    result = batch.run_batch(arguments(sheet, output, config))
    assert [args.sample for args in calls] == ["a", "b"]
    assert result["status"] == "COMPLETE"
    assert result["samples"][0]["n_candidates"] == 1
    assert result["samples"][0]["n_report_rows"] == 2
    assert result["samples"][0]["n_strict_candidates"] == 2
    assert result["samples"][0]["n_fusion_models"] == 1
    assert (output / "samples.tsv").read_bytes() == sheet.read_bytes()
    manifest_bytes = (output / "batch_manifest.json").read_bytes()
    config_bytes = (output / "configs/a.yaml").read_bytes()
    calls.clear()
    resumed = batch.run_batch(arguments(sheet, output, config, "--resume"))
    assert resumed["status"] == "COMPLETE"
    assert len(calls) == 2 and all(args.resume for args in calls)
    assert (output / "batch_manifest.json").read_bytes() == manifest_bytes
    assert (output / "configs/a.yaml").read_bytes() == config_bytes
    assert not (output / ".batch.lock").exists()


@pytest.mark.parametrize("keep_going,statuses,calls_expected", [(False, ["FAILED", "NOT_RUN"], ["a"]), (True, ["FAILED", "COMPLETE"], ["a", "b"])])
def test_failures_record_unrun_samples_and_return_failure(tmp_path, monkeypatch, keep_going, statuses, calls_expected):
    sheet, output, config = fixture(tmp_path)
    calls = []
    monkeypatch.setattr(batch, "run_sample", fake_runner(calls, fail="a"))
    args = arguments(sheet, output, config, *(["--keep-going"] if keep_going else []))
    result = batch.run_batch(args)
    assert result["status"] == "FAILED"
    assert [row["status"] for row in result["samples"]] == statuses
    assert [args.sample for args in calls] == calls_expected
    assert result["samples"][0]["n_candidates"] is None
    assert "explicit native-tool failure" in result["samples"][0]["error"]
    assert not (output / ".batch.lock").exists()


def test_interruption_releases_lock_and_records_not_run_samples(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    monkeypatch.setattr(batch, "run_sample", fake_runner([], interrupt="a"))
    with pytest.raises(KeyboardInterrupt):
        batch.run_batch(arguments(sheet, output, config))
    state = json.loads((output / "batch_state.json").read_text())
    assert state["status"] == "INTERRUPTED"
    assert [row["status"] for row in state["samples"]] == ["INTERRUPTED", "NOT_RUN"]
    assert not (output / ".batch.lock").exists()


@pytest.mark.parametrize("change", ["sheet", "config", "input", "saved_config", "symlink", "lock"])
def test_resume_rejects_changed_provenance_before_state_mutation(tmp_path, monkeypatch, change):
    sheet, output, config = fixture(tmp_path)
    calls = []
    monkeypatch.setattr(batch, "run_sample", fake_runner(calls))
    batch.run_batch(arguments(sheet, output, config))
    previous = (output / "batch_state.json").read_bytes()
    calls.clear()
    if change == "sheet":
        sheet.write_text(sheet.read_text().replace("fungus one", "other fungus"))
    elif change == "config":
        config.write_text(config.read_text().replace("cpu_threads: 2", "cpu_threads: 3"))
    elif change == "input":
        path = config.parent / "proteins.fa"
        path.write_text(path.read_text().replace("MMMM", "AAAA", 1))
    elif change == "saved_config":
        path = output / "configs/a.yaml"
        path.write_text(path.read_text().replace("cpu_threads: 2", "cpu_threads: 3"))
    elif change == "symlink":
        (output / "samples/b/escape").symlink_to(config.parent, target_is_directory=True)
    else:
        (output / ".batch.lock").write_text("another process")
    with pytest.raises(ValueError):
        batch.run_batch(arguments(sheet, output, config, "--resume"))
    assert calls == []
    assert (output / "batch_state.json").read_bytes() == previous


def test_unowned_output_and_symlink_root_are_rejected(tmp_path, monkeypatch):
    sheet, output, config = fixture(tmp_path)
    output.mkdir()
    (output / "unrelated.txt").write_text("preserve")
    monkeypatch.setattr(batch, "run_sample", lambda *args: pytest.fail("runner must not run"))
    with pytest.raises(ValueError, match="not empty"):
        batch.run_batch(arguments(sheet, output, config))
    with pytest.raises(ValueError, match="provenance"):
        batch.run_batch(arguments(sheet, output, config, "--resume"))
    link = tmp_path / "output-link"
    link.symlink_to(output, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        batch.run_batch(arguments(sheet, link, config))
    assert (output / "unrelated.txt").read_text() == "preserve"
