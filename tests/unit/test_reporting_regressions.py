"""Regression checks for summary counts and annotation exports."""
import csv
import logging

import pandas as pd

from funlr.config import FunlrConfig
from funlr.core.context import RunContext, RunPaths
from funlr.core.runner import CommandRunner
from funlr.stages.annotation_export import export_annotations
from funlr.stages.rescue_v3 import run_miniprot


def context(tmp_path, inputs=None):
    cfg = FunlrConfig({
        "inputs": inputs or {},
        "execution": {"output_dir": str(tmp_path / "run")},
        "reporting": {"plots": False, "integration": False},
    })
    return RunContext(cfg, RunPaths(tmp_path / "run"),
                      CommandRunner({}, tmp_path / "commands.jsonl"), None,
                      logging.getLogger("reporting-regression-test"))


def test_cds_only_miniprot_output_is_missing_from_hit_audit_and_still_refined(tmp_path, monkeypatch):
    ctx = context(tmp_path)
    fasta = ">Q\nAAAA\n>R\nAAAA\n"
    (ctx.paths.stage_dir(0) / "proteins_clean.faa").write_text(fasta)
    (ctx.paths.stage_dir(4) / "rescue_priority.faa").write_text(fasta)
    tables = ctx.paths.stage_dir(3) / "tables"
    tables.mkdir()
    pd.DataFrame([
        {"protein_id": pid, "tier": "TIER_1A_HIGH_CONFIDENCE",
         "flags": "", "nbd_confidence": "HIGH"}
        for pid in ("Q", "R")
    ]).to_csv(tables / "tiered_candidates.tsv", sep="\t", index=False)

    def emit_partial_gff(_argv, *, stdout_path, stderr_path):
        stdout_path.write_text(
            "##gff-version 3\n"
            "ctg\tminiprot\tCDS\t1\t12\t99\t+\t0\tParent=MP1;Target=Q 1 4\n"
            "ctg\tminiprot\tmRNA\t101\t112\t99\t+\t.\tID=MP2;Target=R 1 4\n"
            "ctg\tminiprot\tCDS\t101\t112\t99\t+\t0\tParent=MP2;Target=R 1 4\n"
        )
        stderr_path.write_text("")

    monkeypatch.setattr(ctx.runner, "resolve", lambda name: name)
    monkeypatch.setattr(ctx.runner, "run", emit_partial_gff)
    run_miniprot(ctx)
    for root in (ctx.paths.stage_dir(5), ctx.paths.stage_dir(5) / "comprehensive"):
        assert (root / "summaries/miniprot_missing_hits.txt").read_text() == "Q\n"
        assert (root / "summaries/miniprot_extra_hits.txt").read_text() == ""
        assert (root / "queries/exonerate_refine_ids.txt").read_text() == "Q\n"
        reasons = pd.read_csv(root / "qc/exonerate_refine_reasons.tsv", sep="\t")
        assert reasons.to_dict("records") == [{"query": "Q", "reason": "no_mrna"}]
        # Retain all feature evidence in the ordinary summary for inspection.
        summary = pd.read_csv(root / "tables/miniprot_summary.tsv", sep="\t")
        assert set(summary["query"]) == {"Q", "R"}


def test_joined_annotation_tsv_uses_lf_without_altering_input_or_values(tmp_path):
    source = tmp_path / "annotations.txt"
    original = b"GeneID\tDescription\r\nG\tknown protein\r\nU\tunrelated protein\r\n"
    source.write_bytes(original)
    ctx = context(tmp_path, {"annotations": str(source)})
    report = pd.DataFrame([{
        "protein_id": "Q", "gene_id": "G", "tier": "TIER_1A_HIGH_CONFIDENCE",
        "nbd_confidence": "HIGH", "rescue_priority": 0, "flags": "",
        "domains_grouped": "NACHT;ANK", "NLR_architecture": "NACHT;ANK",
    }])
    target = ctx.paths.final_dir / "integration"
    export_annotations(ctx, report, target)
    result = target / "annotations.withNLR.tsv"
    assert source.read_bytes() == original
    assert b"\r" not in result.read_bytes()
    assert result.read_bytes().endswith(b"\n")
    with result.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["Description"] for row in rows] == ["known protein", "unrelated protein"]
    assert rows[0]["NLR_tier"] == "TIER_1A_HIGH_CONFIDENCE"
    assert rows[0]["NLR_architecture"] == "NACHT;ANK"
    assert rows[1]["NLR_tier"] == ""
