"""Regression tests for the observed reference rescue defects and updated tracks."""
from pathlib import Path
import logging

import pandas as pd
import pytest

from funlr.config import FunlrConfig
from funlr.core.context import RunContext, RunPaths
from funlr.core.runner import CommandRunner
from funlr.core.errors import StageError
from funlr.stages.rescue_v3 import (
    FEATURE_COLUMNS, fusion_candidates, fusion_evidence, miniprot_features,
    refinement_reasons, run_miniprot, run_exonerate,
)
from funlr.stages.reporting_v3 import REPORT_COLUMNS, integrate_fusions, strict_sets, run_final, domain_plots
from funlr.stages.annotation_export import tagged_gff, rewrite_protein_headers


def member(pid, start, end, nbd, sensor, strand="+"):
    return {"protein_id": pid, "gene_id": pid, "start": start, "end": end, "scaffold": "ctg", "strand": strand, "has_nbd": nbd, "has_sensor": sensor, "nbd_confidence": "HIGH" if nbd else "LOW", "flags": "PASS", "tier": "TIER_1A_HIGH_CONFIDENCE", "protein_length": 100}


def test_reference_fragment_only_alignment_is_rejected():
    raw = Path(__file__).parents[1] / "fixtures/rescue/fragment_only.vulgar.txt"
    hypothesis = {"fusion_id": "PRIORITY__FUSION_000003", "members": "psicub_aus-FUN_003816-T1;psicub_aus-FUN_003851-T1", "member_lengths": "675;196", "member_starts": "39649;132795", "member_ends": "42284;133449", "scaffold": "scaffold_49", "strand": "+"}
    result = fusion_evidence(hypothesis, raw)
    assert not result["accepted"]
    assert result["reason"] == "FRAGMENT_ONLY_OR_WRONG_LOCUS"
    assert result["aligned_query_end"] == 669
    assert result["query_span_fraction"] == pytest.approx(669 / 871)


def test_native_minus_strand_alignment_spans_both_supplied_members():
    raw = Path(__file__).parents[1] / "fixtures/rescue/minus_member_spanning.vulgar.txt"
    hypothesis = {"fusion_id": "PRIORITY__FUSION_000003", "members": "psicub_aus-FUN_004049-T1;psicub_aus-FUN_004048-T1", "member_lengths": "798;700", "member_starts": "40471;38182", "member_ends": "43381;40338", "scaffold": "scaffold_54", "strand": "-"}
    result = fusion_evidence(hypothesis, raw)
    assert result["accepted"]
    assert result["query_span_fraction"] == 1.0
    assert (result["start"], result["end"], result["strand"]) == (38185, 43381, "-")


def test_fusion_members_follow_transcript_direction_on_both_strands():
    seqs = {"N": "N" * 100, "S": "S" * 100}
    for strand, rows in (("+", [member("N", 1, 300, True, False), member("S", 1001, 1300, False, True)]), ("-", [member("S", 1, 300, False, True, "-"), member("N", 1001, 1300, True, False, "-")])):
        hypotheses, queries = fusion_candidates(pd.DataFrame(rows), seqs, {}, "priority")
        assert len(hypotheses) == 1
        assert hypotheses.iloc[0]["members"] == "N;S"
        assert hypotheses.iloc[0]["strand"] == strand
        assert next(iter(queries.values())) == seqs["N"] + seqs["S"]


def test_missing_fusion_member_is_fatal():
    with pytest.raises(StageError, match="Fusion member"):
        fusion_candidates(pd.DataFrame([member("N", 1, 300, True, False), member("S", 1001, 1300, False, True)]), {"N": "N" * 100}, {}, "priority")


def hypothesis():
    return {"fusion_id": "PRIORITY__FUSION_000001", "members": "N;S", "member_lengths": "100;100", "member_starts": "1;1001", "member_ends": "300;1300", "scaffold": "ctg", "strand": "+"}


def test_single_alignment_must_match_all_members_and_source_loci(tmp_path):
    raw = tmp_path / "alignment.txt"
    raw.write_text("vulgar: PRIORITY__FUSION_000001 0 200 . ctg 0 1300 + 500 M 100 300 I 0 700 M 100 300\n")
    assert fusion_evidence(hypothesis(), raw)["accepted"]
    raw.write_text("vulgar: PRIORITY__FUSION_000001 0 200 . wrong 0 1300 + 500 M 100 300 I 0 700 M 100 300\n")
    assert not fusion_evidence(hypothesis(), raw)["accepted"]
    # Bounds span both members, but the second protein and its locus are gaps.
    raw.write_text("vulgar: PRIORITY__FUSION_000001 0 200 . ctg 0 1300 + 500 M 100 300 G 100 0 I 0 1000\n")
    assert not fusion_evidence(hypothesis(), raw)["accepted"]


def test_track_duplicates_collapse_and_unrelated_strict_candidate_survives(tmp_path):
    rows = [member("N", 1, 300, True, False), member("S", 1001, 1300, False, True), member("U", 2001, 2300, True, True)]
    original = pd.DataFrame(rows).reindex(columns=REPORT_COLUMNS)
    original["mp_scaffold"] = "ctg"
    original["comprehensive_mp_scaffold"] = "ctg"
    raw = tmp_path / "alignment"
    raw.write_text("vulgar: PRIORITY__FUSION_000001 0 200 . ctg 0 1300 + 500 M 100 300 I 0 700 M 100 300\n")
    pri = fusion_evidence(hypothesis(), raw)
    comp = {**pri, "fusion_id": "COMPREHENSIVE__FUSION_000001", "track": "comprehensive"}
    out = integrate_fusions(original, pd.DataFrame([{**pri, "track": "priority"}, comp]))
    assert out["is_fusion_model"].sum() == 1
    for frame in strict_sets(out).values():
        assert len(frame) == 2
        assert "U" in set(frame["protein_id"])
    rejected = integrate_fusions(original, pd.DataFrame([{**pri, "track": "priority", "accepted": False}]))
    assert len(strict_sets(rejected)["combined"]) == 3
    assert rejected["fused_into"].eq("").all()


def test_miniprot_real_attributes_and_empty_refinement(tmp_path):
    path = tmp_path / "mp.gff"
    path.write_text("ctg\tminiprot\tmRNA\t1\t100\t99\t+\t.\tID=MP1;Rank=1;Identity=0.9800;Positive=0.9900;Target=Q 1 33\n")
    result = miniprot_features(path)
    assert result.iloc[0][["identity", "positive", "rank"]].tolist() == ["0.9800", "0.9900", "1"]
    assert refinement_reasons(["Q", "R"], pd.DataFrame(columns=FEATURE_COLUMNS), pd.DataFrame()) == {"Q": "no_miniprot_hit", "R": "no_miniprot_hit"}


def context(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">ctg\nAAA\n")
    cfg = FunlrConfig({"inputs": {"genome": str(genome)}, "execution": {"output_dir": str(tmp_path / "run")}, "reporting": {"plots": False, "integration": False}})
    paths = RunPaths(tmp_path / "run")
    return RunContext(cfg, paths, CommandRunner({}, tmp_path / "commands.jsonl"), None, logging.getLogger("test"))


def test_comprehensive_only_and_all_no_hit_stable_outputs(tmp_path):
    ctx = context(tmp_path)
    stage0, stage3, stage4 = (ctx.paths.stage_dir(i) for i in (0, 3, 4))
    (stage0 / "proteins_clean.faa").write_text(">Q\nAAAA\n")
    (stage3 / "tables").mkdir()
    pd.DataFrame([{**member("Q", 1, 12, True, True), "tier": "TIER_1A_HIGH_CONFIDENCE"}]).to_csv(stage3 / "tables/tiered_candidates.tsv", sep="\t", index=False)
    (stage4 / "rescue_priority.faa").write_text("")
    binary = tmp_path / "miniprot"
    binary.write_text("#!/bin/sh\nprintf '##gff-version 3\\n'\n")
    binary.chmod(0o755)
    ctx.runner.tools["miniprot"] = str(binary)
    run_miniprot(ctx)
    assert (ctx.paths.stage_dir(5) / "queries/exonerate_refine_ids.txt").read_text() == ""
    assert (ctx.paths.stage_dir(5) / "comprehensive/queries/exonerate_refine_ids.txt").read_text() == "Q\n"
    # Skip the comprehensive query to test both empty Stage6 tracks coherently.
    (ctx.paths.stage_dir(5) / "comprehensive/queries/exonerate_refine_ids.txt").write_text("")
    run_exonerate(ctx)
    run_final(ctx)
    assert len(pd.read_csv(ctx.paths.final_dir / "nlr_strict_candidates.tsv", sep="\t")) == 1
    assert pd.read_csv(ctx.paths.final_dir / "nlr_strict_candidates.priority.tsv", sep="\t").empty
    assert pd.read_csv(ctx.paths.final_dir / "nlr_strict_candidates.comprehensive.tsv", sep="\t").empty


def test_annotation_tags_and_headers_preserve_original_models(tmp_path):
    source, tagged, raw, final = (tmp_path / name for name in ("in.gff3", "out.gff3", "raw.fa", "final.fa"))
    source.write_text("##gff-version 3\nctg\tx\tmRNA\t1\t300\t.\t+\t.\tID=T;Parent=G;product=protein%20one\nctg\tx\tCDS\t1\t300\t.\t+\t0\tParent=T\n")
    tagged_gff(source, tagged, {"G": {"NLR_tier": "TIER_1A_HIGH_CONFIDENCE", "NLR_architecture": "NACHT;ANK"}})
    assert source.read_text().splitlines()[-1] == tagged.read_text().splitlines()[-1]
    raw.write_text(">T\nAAAA\n")
    rewrite_protein_headers(tagged, raw, final)
    assert final.read_text() == ">T [product=protein one] [NLR_tier=TIER_1A_HIGH_CONFIDENCE] [NLR_architecture=NACHT;ANK]\nAAAA\n"


def test_zero_exit_from_r_without_pngs_is_a_failure(tmp_path):
    ctx = context(tmp_path)
    ctx.config["reporting"]["plot_backend"] = "r"
    hits = pd.DataFrame([{"protein_id": "Q", "domain": "NACHT", "start": 1, "end": 100}])
    hits.to_csv(ctx.paths.stage_dir(2) / "all_domain_hits.tsv", sep="\t", index=False)
    stage7 = ctx.paths.stage_dir(7)
    for name in ("logs", "metadata"):
        (stage7 / name).mkdir()
    binary = tmp_path / "Rscript"
    binary.write_text("#!/bin/sh\necho 'Warning: failed to load graphics DLL' >&2\nexit 0\n")
    binary.chmod(0o755)
    ctx.runner.tools["Rscript"] = str(binary)
    tiered = pd.DataFrame([member("Q", 1, 300, True, True)])
    work = stage7 / "plot_work"
    (work / "individual").mkdir(parents=True)
    for name in ("nlr_domains_by_tier.png", "sensor_distribution.png", "nbd_distribution.png", "asm_presence.png", "length_vs_priority.png", "order_sanity.png", "individual/Q.png"):
        (work / name).write_bytes(b"\x89PNG\r\n\x1a\nold plot")
    with pytest.raises(StageError, match="expected PNG"):
        domain_plots(ctx, tiered, stage7, ctx.paths.final_dir)
