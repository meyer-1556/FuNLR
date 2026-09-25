"""Unit tests for stages 5 (miniprot), 6 (exonerate), 7 (finalize).

Pure logic is tested without a RunContext; empty-input stable-output paths are
tested end-to-end with shim executables standing in for the external tools.
"""

from __future__ import annotations

import logging
import json
import shutil
from collections import OrderedDict
from pathlib import Path

import pandas as pd
import pytest

from funlr.config import FunlrConfig
from funlr.core.context import RunContext, RunPaths
from funlr.core.errors import StageError
from funlr.core.runner import CommandRunner
from funlr.parsers.fasta import dedup_keep_first, fasta_ids
from funlr.parsers.gff import (
    miniprot_locus_summary,
    read_exonerate_features,
    read_miniprot_features,
)
from funlr.stages import stage5_miniprot as s5
from funlr.stages import stage6_exonerate as s6
from funlr.stages import stage7_finalize as s7

FIX = Path(__file__).resolve().parents[1] / "fixtures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_shim(bin_dir: Path, name: str, body: str) -> Path:
    bin_dir.mkdir(parents=True, exist_ok=True)
    p = bin_dir / name
    p.write_text("#!/bin/sh\n" + body + "\n")
    p.chmod(0o755)
    return p


def _version_shim(bin_dir: Path, name: str) -> Path:
    return _make_shim(bin_dir, name, f'echo "{name}-shim 0.1"')


def _make_ctx(
    tmp_path: Path,
    tools: dict[str, str] | None = None,
    include_tier3: bool = False,
    dry_run: bool = False,
    exon_maxn: int = 200,
) -> RunContext:
    out = tmp_path / "out"
    genome = tmp_path / "genome.fa"
    genome.write_text(">ctg1\nACGTACGT\n")
    cfg = FunlrConfig(
        {
            "sample": {"sample_id": "testsample", "species_id": "Test_sp"},
            "inputs": {"genome": str(genome)},
            "rescue": {"include_tier3": include_tier3, "exon_maxn": exon_maxn},
            "execution": {"output_dir": str(out), "cpu_threads": 4},
            "tools": tools or {},
            "reporting": {"plots": False, "integration": False},
        }
    )
    logger = logging.getLogger("funlr.test")
    paths = RunPaths(output_dir=out)
    runner = CommandRunner(tools or {}, out / "logs" / "commands.jsonl", logger, dry_run=dry_run)
    return RunContext(
        config=cfg, paths=paths, runner=runner, state=None, logger=logger, dry_run=dry_run
    )


def _tier_fastas(dir_: Path) -> "OrderedDict[str, Path]":
    return OrderedDict(
        [
            ("tier2A", dir_ / "tier2A_high_priority_rescue.faa"),
            ("tier2B", dir_ / "tier2B_rescue_candidate.faa"),
            ("tier1B", dir_ / "tier1B_needs_review.faa"),
            ("tier3", dir_ / "tier3_architectural_variant.faa"),
        ]
    )


def _copy_tier_fixtures(dir_: Path) -> "OrderedDict[str, Path]":
    dir_.mkdir(parents=True, exist_ok=True)
    mapping = {
        "tier2A": "s05_tier2A.faa",
        "tier2B": "s05_tier2B.faa",
        "tier1B": "s05_tier1B.faa",
        "tier3": "s05_tier3.faa",
    }
    tiers = _tier_fastas(dir_)
    for label, src in mapping.items():
        shutil.copy(FIX / src, tiers[label])
    return tiers


# ---------------------------------------------------------------------------
# Stage 5: query concatenation
# ---------------------------------------------------------------------------

def test_concat_query_fastas_order_and_tier3_toggle(tmp_path):
    tiers = _copy_tier_fixtures(tmp_path)

    out = tmp_path / "raw_no_t3.faa"
    s5.concat_query_fastas(tiers, out, include_tier3=False)
    # Legacy order: tier2A, tier2B, tier1B (tier3 excluded by default).
    assert fasta_ids(out) == ["p1", "p2", "p2", "p3", "p4"]

    out3 = tmp_path / "raw_t3.faa"
    s5.concat_query_fastas(tiers, out3, include_tier3=True)
    assert fasta_ids(out3) == ["p1", "p2", "p2", "p3", "p4", "p5"]


def test_concat_query_fastas_skips_missing_and_empty(tmp_path):
    tiers = _copy_tier_fixtures(tmp_path)
    tiers["tier2B"].unlink()  # missing
    tiers["tier1B"].write_text("")  # empty
    out = tmp_path / "raw.faa"
    s5.concat_query_fastas(tiers, out, include_tier3=True)
    assert fasta_ids(out) == ["p1", "p2", "p5"]


# ---------------------------------------------------------------------------
# Stage 5: dedup + counts + refine ids
# ---------------------------------------------------------------------------

def test_dedup_keep_first_counts_and_stability(tmp_path):
    tiers = _copy_tier_fixtures(tmp_path)
    raw = tmp_path / "raw.faa"
    s5.concat_query_fastas(tiers, raw, include_tier3=False)
    dedup = tmp_path / "dedup.faa"
    total, kept = dedup_keep_first(raw, dedup)
    assert (total, kept) == (5, 4)
    assert fasta_ids(dedup) == ["p1", "p2", "p3", "p4"]
    # First occurrence wins: p2 sequence comes from tier2A (MCCC), not tier2B (MGGG).
    text = dedup.read_text()
    assert ">p2\nMCCC\n" in text
    assert "MGGG" not in text
    # Header description of first record is preserved.
    assert ">p1 some description\n" in text


def test_query_counts_table(tmp_path):
    tiers = _copy_tier_fixtures(tmp_path)
    raw = tmp_path / "raw.faa"
    s5.concat_query_fastas(tiers, raw, include_tier3=False)
    dedup = tmp_path / "dedup.faa"
    dedup_keep_first(raw, dedup)
    files = OrderedDict(
        [
            ("tier2A", tiers["tier2A"]),
            ("tier2B", tiers["tier2B"]),
            ("tier1B", tiers["tier1B"]),
            ("tier3", tiers["tier3"]),
            ("raw_concat", raw),
            ("dedup", dedup),
        ]
    )
    out = tmp_path / "query_counts.tsv"
    s5.write_query_counts(files, out)
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["label", "file", "seqs"]
    assert list(df["label"]) == ["tier2A", "tier2B", "tier1B", "tier3", "raw_concat", "dedup"]
    assert list(df["seqs"]) == [2, 2, 1, 1, 5, 4]


def test_miniprot_features_to_locus_summary():
    feats = read_miniprot_features(FIX / "s05_miniprot.gff3")
    assert list(feats["protein_id"]) == ["p2", "p2", "p2", "p9"]
    loc = miniprot_locus_summary(feats).set_index("protein_id")
    p2 = loc.loc["p2"]
    assert (p2["scaffold"], p2["start"], p2["end"], p2["strand"], p2["n_features"]) == (
        "ctg1", 100, 900, "+", 3,
    )
    p9 = loc.loc["p9"]
    assert (p9["scaffold"], p9["start"], p9["end"], p9["strand"], p9["n_features"]) == (
        "ctg2", 50, 200, "-", 1,
    )


def test_refine_ids_sorted_difference():
    assert s5.compute_refine_ids(["c", "a", "b"], {"b"}) == ["a", "c"]
    assert s5.compute_refine_ids([], {"x"}) == []
    assert s5.compute_refine_ids(["a"], {"a"}) == []


def test_hit_ids_from_summary(tmp_path):
    summary = tmp_path / "miniprot_summary.tsv"
    summary.write_text("protein_id\tscaffold\tstart\tend\tstrand\np2\tctg1\t1\t9\t+\n\np9\tctg2\t5\t8\t-\n")
    assert s5.hit_ids_from_summary(summary) == {"p2", "p9"}
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    assert s5.hit_ids_from_summary(empty) == set()
    assert s5.hit_ids_from_summary(tmp_path / "missing.tsv") == set()


# ---------------------------------------------------------------------------
# Stage 6: exonerate parsing / merging / caps
# ---------------------------------------------------------------------------

def test_exonerate_pid_regex_recovery():
    feats = read_exonerate_features(FIX / "s06_exonerate_per_query.gff")
    # Comment lines filtered; 'sequence p1 ;' and 'Query= p2' forms recovered.
    assert list(feats["protein_id"]) == ["p1", "p1", "p2"]
    assert list(feats["feature"]) == ["gene", "exon", "gene"]


def test_exonerate_pid_fallback_forms():
    feats = read_exonerate_features(FIX / "s06_exonerate_fallback.gff")
    # First identifier-like token fallback; row with non-integer coords skipped.
    assert list(feats["protein_id"]) == ["XYZ123"]
    assert feats.iloc[0]["scaffold"] == "ctg3"


def test_merge_gffs_filters_comment_lines(tmp_path):
    merged = tmp_path / "merged.gff3"
    n = s6.merge_gffs(
        [("p1", FIX / "s06_exonerate_per_query.gff"), ("p2", FIX / "s06_exonerate_per_query.gff")], merged
    )
    assert n == 6
    lines = merged.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert len(lines) == 7  # header + 6 feature lines
    assert sum(1 for ln in lines if ln.startswith("#")) == 1
    assert read_exonerate_features(merged)["protein_id"].tolist() == ["p1"] * 3 + ["p2"] * 3


def test_apply_exon_maxn():
    ids = ["a", "b", "c"]
    assert s6.apply_exon_maxn(ids, 0) == ["a", "b", "c"]
    assert s6.apply_exon_maxn(ids, 2) == ["a", "b"]
    assert s6.apply_exon_maxn(ids, 200) == ["a", "b", "c"]


# ---------------------------------------------------------------------------
# Stage 7: loaders, merge, strict set, QC, BEDs
# ---------------------------------------------------------------------------

def test_load_miniprot_summary_rename_and_empty(tmp_path):
    mp = s7.load_miniprot_summary(FIX / "s07_miniprot_summary.tsv")
    assert list(mp.columns) == ["protein_id", "mp_scaffold", "mp_start", "mp_end", "mp_strand"]
    row = mp.set_index("protein_id").loc["p2"]
    assert (row["mp_scaffold"], row["mp_start"], row["mp_end"], row["mp_strand"]) == (
        "ctg1", 1900, 2700, "-",
    )
    empty = s7.load_miniprot_summary(tmp_path / "missing.tsv")
    assert empty.empty
    assert list(empty.columns) == ["protein_id", "mp_scaffold", "mp_start", "mp_end", "mp_strand"]


def test_load_exonerate_summary_direct():
    ex = s7.load_exonerate_summary(FIX / "s07_exonerate_summary.tsv", FIX / "s07_exonerate_merged.gff3")
    # 'junk' column dropped by the legacy keep-list.
    assert list(ex.columns) == [
        "protein_id", "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features",
    ]
    row = ex.set_index("protein_id").loc["p3"]
    assert (row["exo_scaffold"], row["exo_start"], row["exo_end"], row["exo_strand"]) == (
        "ctg2", 40, 410, "+",
    )


def test_load_exonerate_summary_gff_fallback(tmp_path):
    # Summary missing/empty -> parse merged GFF via Target=/Parent= regexes.
    ex = s7.load_exonerate_summary(
        tmp_path / "absent.tsv", FIX / "s07_exonerate_merged.gff3"
    ).set_index("protein_id")
    assert set(ex.index) == {"p3", "p7"}  # 'no_id_here' line skipped
    p3 = ex.loc["p3"]
    assert (p3["exo_scaffold"], p3["exo_start"], p3["exo_end"], p3["exo_strand"]) == (
        "ctg2", 40, 410, "+",
    )
    p7 = ex.loc["p7"]
    assert (p7["exo_scaffold"], p7["exo_start"], p7["exo_end"], p7["exo_strand"]) == (
        "ctg5", 10, 90, "-",
    )
    assert pd.isna(p3["n_features"])

    # Neither summary nor gff -> empty frame with legacy columns.
    empty = s7.load_exonerate_summary(tmp_path / "absent.tsv", tmp_path / "absent.gff3")
    assert empty.empty
    assert list(empty.columns) == [
        "protein_id", "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features",
    ]


def _merged_report() -> pd.DataFrame:
    tiered = pd.read_csv(FIX / "s07_tiered.tsv", sep="\t")
    mp = s7.load_miniprot_summary(FIX / "s07_miniprot_summary.tsv")
    ex = s7.load_exonerate_summary(FIX / "s07_exonerate_summary.tsv", FIX / "s07_exonerate_merged.gff3")
    return s7.merge_final_report(tiered, mp, ex)


def test_merge_final_report_column_selection():
    out = _merged_report()
    # Only WANT columns, in WANT order; 'custom_note' dropped.
    assert list(out.columns) == [
        "protein_id", "gene_id", "scaffold", "start", "end", "strand",
        "tier", "flags", "domains_grouped",
        "mp_scaffold", "mp_start", "mp_end", "mp_strand",
        "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features",
    ]
    by_pid = out.set_index("protein_id")
    assert by_pid.loc["p2", "mp_start"] == 1900
    assert pd.isna(by_pid.loc["p1", "mp_start"])
    assert by_pid.loc["p3", "exo_start"] == 40
    assert pd.isna(by_pid.loc["p1", "exo_scaffold"])


def test_merge_final_report_requires_protein_id():
    bad = pd.DataFrame({"gene_id": ["g1"]})
    with pytest.raises(StageError):
        s7.merge_final_report(bad, pd.DataFrame(columns=["protein_id"]), pd.DataFrame(columns=["protein_id"]))


@pytest.mark.parametrize(
    "value,expected",
    [
        (float("nan"), False),
        ("PASS", False),
        ("", False),
        ("CONTIG_END", False),
        ("NON_NLR_ANNOT", True),
        ("CONTIG_END|NON_NLR_ANNOT", True),
        ("A;NON_NLR_ANNOT,B", True),
        # Legacy substring fallback quirk (line 2294): infix matches count.
        ("XNON_NLR_ANNOTX", True),
        ("NON_NLR_ANNOT_LIKE", True),
    ],
)
def test_has_non_nlr_edge_cases(value, expected):
    assert s7.has_non_nlr(value) is expected


def test_strict_candidates_filter():
    out = _merged_report()
    strict = s7.strict_candidates(out)
    # TIER_4A excluded by tier; p3 excluded via NON_NLR_ANNOT flag; p5 (empty flags) kept.
    assert set(strict["protein_id"]) == {"p1", "p2", "p5"}


def test_qc_summary_metric_rows():
    tiered = pd.read_csv(FIX / "s07_tiered.tsv", sep="\t")
    strict = s7.strict_candidates(_merged_report())
    qc = s7.build_qc_summary(tiered, strict)
    metrics = dict(zip(qc["metric"], qc["value"]))
    assert metrics["n_total_candidates"] == 5
    for tier in [
        "TIER_1A_HIGH_CONFIDENCE",
        "TIER_2A_HIGH_PRIORITY_RESCUE",
        "TIER_2B_RESCUE_CANDIDATE",
        "TIER_4A_HOUSEKEEPING",
        "TIER_3_ARCHITECTURAL_VARIANT",
    ]:
        assert metrics[f"tier_count::{tier}"] == 1
    assert metrics["has_nbd_true"] == 3        # TRUE strings: p1, p2, p5
    assert metrics["has_repeat_true"] == 2     # inferred from domains: WD40 (p1), TPR (p5)
    assert metrics["order_ok_true"] == 2       # TRUE strings: p1, p3 ('NA' -> False)
    assert metrics["n_strict_candidates"] == 3
    assert metrics["flag::CONTIG_END"] == 2    # p2, p3 ('PASS' dropped)
    assert metrics["flag::NON_NLR_ANNOT"] == 2  # p3, p4
    # Row order: totals first, flags last.
    assert qc["metric"].iloc[0] == "n_total_candidates"
    # Tied counts keep first-occurrence order (legacy sorted() is stable).
    assert qc["metric"].iloc[-2:].tolist() == ["flag::CONTIG_END", "flag::NON_NLR_ANNOT"]


def test_qc_summary_na_paths():
    tiered = pd.DataFrame(
        {"protein_id": ["x", "y"], "tier": ["TIER_1A_HIGH_CONFIDENCE"] * 2}
    )
    qc = s7.build_qc_summary(tiered, tiered)
    metrics = dict(zip(qc["metric"], qc["value"]))
    assert metrics["n_total_candidates"] == 2
    assert metrics["has_nbd_true"] == "NA"     # no column and no domains to infer from
    assert metrics["has_repeat_true"] == "NA"
    assert metrics["order_ok_true"] == 0       # all-NA series -> counts only True
    assert not any(m.startswith("flag::") for m in metrics)  # no flags column


def test_bed_line_formats():
    out = _merged_report()
    orig = s7.bed_lines_original(out)
    assert len(orig) == 5
    # 0-based start, name = pid|tier, score 0, strand from source.
    assert orig[0] == "ctg1\t99\t900\tp1|TIER_1A_HIGH_CONFIDENCE\t0\t+\n"
    mp = s7.bed_lines_miniprot(out)
    assert mp == ["ctg1\t1899\t2700\tp2|MINIPROT\t0\t-\n"]
    exo = s7.bed_lines_exonerate(out)
    assert exo == ["ctg2\t39\t410\tp3|EXONERATE\t0\t+\n"]


def test_tier_summary_frame():
    summ = s7.tier_summary(_merged_report())
    assert list(summ.columns) == ["tier", "count"]
    assert dict(zip(summ["tier"], summ["count"]))["TIER_1A_HIGH_CONFIDENCE"] == 1
    empty = s7.tier_summary(pd.DataFrame({"protein_id": ["x"]}))
    assert list(empty.columns) == ["tier", "count"]
    assert empty.empty


# ---------------------------------------------------------------------------
# Stage 5: end-to-end empty-input stable-output path (shim miniprot)
# ---------------------------------------------------------------------------

def test_stage5_empty_queries_stable_outputs(tmp_path):
    miniprot = _version_shim(tmp_path / "bin", "miniprot")
    ctx = _make_ctx(tmp_path, tools={"miniprot": str(miniprot)})
    # No stage 4 FASTAs at all -> empty raw concat -> stable empty outputs.
    outs = s5.run(ctx)
    stage5 = ctx.paths.stage_dir(5)
    gff = stage5 / "gff" / "miniprot.gff3"
    summary = stage5 / "tables" / "miniprot_summary.tsv"
    refine = stage5 / "queries" / "exonerate_refine_ids.txt"
    manifest = stage5 / "stage5_manifest.tsv"
    run_info = stage5 / "metadata" / "run_info.json"
    assert gff.read_text() == "##gff-version 3\n"
    assert pd.read_csv(summary, sep="\t").empty
    assert refine.read_text() == ""
    assert (stage5 / "queries" / "rescue_queries.faa").read_text() == ""
    assert all(p.is_file() for p in outs)
    assert manifest.is_file()
    labels = [ln.split("\t")[0] for ln in manifest.read_text().splitlines()]
    assert "tables/miniprot_summary.tsv" in labels
    assert "queries/exonerate_refine_ids.txt" in labels
    info = json.loads(run_info.read_text())
    assert info["stage"] == 5 and info["track"] == "priority"
    assert info["query_count"] == 0
    assert refine in outs


def test_stage5_full_run_with_shim_miniprot(tmp_path):
    bin_dir = tmp_path / "bin"
    miniprot = _make_shim(
        bin_dir,
        "miniprot",
        'if [ "$1" = "--version" ]; then echo "miniprot-shim 0.1"; exit 0; fi\n'
        f'cat "{FIX}/s05_miniprot.gff3"',
    )
    ctx = _make_ctx(tmp_path, tools={"miniprot": str(miniprot)})
    _copy_tier_fixtures(ctx.paths.stage_dir(4))

    s5.run(ctx)
    stage5 = ctx.paths.stage_dir(5)

    qc = pd.read_csv(stage5 / "summaries" / "query_counts.tsv", sep="\t")
    assert list(qc["seqs"]) == [5, 4]  # raw versus deduplicated; Tier3 excluded

    feats = pd.read_csv(stage5 / "tables" / "miniprot_features.tsv", sep="\t")
    assert len(feats) == 4
    summary = pd.read_csv(stage5 / "tables" / "miniprot_summary.tsv", sep="\t")
    assert set(summary["query"]) == {"p2", "p9"}
    counts = pd.read_csv(stage5 / "tables" / "miniprot_hit_counts.tsv", sep="\t")
    assert dict(zip(counts["feature"], counts["count"])) == {"mRNA": 2, "CDS": 2}

    # Refine list = dedup query ids minus miniprot-hit ids, sorted.
    refine = (stage5 / "queries" / "exonerate_refine_ids.txt").read_text().splitlines()
    assert refine == ["p1", "p3", "p4"]
    # miniprot stderr captured to logs/miniprot.log.
    assert (stage5 / "logs" / "miniprot.log").is_file()


def test_stage5_dry_run_plans_without_outputs(tmp_path):
    miniprot = _version_shim(tmp_path / "bin", "miniprot")
    ctx = _make_ctx(tmp_path, tools={"miniprot": str(miniprot)}, dry_run=True)
    _copy_tier_fixtures(ctx.paths.stage_dir(4))
    outs = s5.run(ctx)
    assert outs == []
    stage5 = ctx.paths.stage_dir(5)
    # The command runner's dry-run path writes neither output nor a command log.
    assert not (stage5 / "gff" / "miniprot.gff3").exists()
    assert not (tmp_path / "out" / "logs" / "commands.jsonl").exists()


# ---------------------------------------------------------------------------
# Stage 6: end-to-end empty-input stable-output paths (shim seqkit/exonerate)
# ---------------------------------------------------------------------------

def _stage6_ctx(tmp_path: Path, dry_run: bool = False) -> RunContext:
    bin_dir = tmp_path / "bin"
    tools = {
        "seqkit": str(_version_shim(bin_dir, "seqkit")),
        "exonerate": str(_version_shim(bin_dir, "exonerate")),
    }
    ctx = _make_ctx(tmp_path, tools=tools, dry_run=dry_run)
    (ctx.paths.stage_dir(0) / "proteins_clean.faa").write_text(">p1\nMAAA\n")
    (ctx.paths.stage_dir(5) / "queries").mkdir(parents=True, exist_ok=True)
    (ctx.paths.stage_dir(5) / "queries" / "rescue_queries.faa").write_text(">p1\nMAAA\n")
    return ctx


def test_stage6_missing_refine_ids_stable_outputs(tmp_path):
    ctx = _stage6_ctx(tmp_path)
    outs = s6.run(ctx)
    stage6 = ctx.paths.stage_dir(6)
    for rel in ["gff/exonerate_merged.gff3", "tables/exonerate_summary.tsv",
                "summaries/exonerate_counts.tsv"]:
        p = stage6 / rel
        assert p.is_file()
        if p.suffix == ".tsv":
            assert pd.read_csv(p, sep="\t").empty
        else:
            assert p.read_text() == "##gff-version 3\n"
    assert (stage6 / "stage6_manifest.tsv").is_file()
    info = json.loads((stage6 / "metadata" / "run_info.json").read_text())
    assert info["stage"] == 6 and info["selected_original"] == 0
    assert all(path.is_file() for path in outs)
    assert (stage6 / "comprehensive/tables/exonerate_summary.tsv").is_file()


def test_stage6_empty_refine_ids_stable_outputs(tmp_path):
    ctx = _stage6_ctx(tmp_path)
    (ctx.paths.stage_dir(5) / "queries" / "exonerate_refine_ids.txt").write_text("")
    s6.run(ctx)
    stage6 = ctx.paths.stage_dir(6)
    assert (stage6 / "gff" / "exonerate_merged.gff3").read_text() == "##gff-version 3\n"
    assert (stage6 / "stage6_manifest.tsv").is_file()


def test_stage6_missing_rescue_fasta_raises(tmp_path):
    ctx = _stage6_ctx(tmp_path)
    (ctx.paths.stage_dir(5) / "queries" / "rescue_queries.faa").write_text("")
    (ctx.paths.stage_dir(5) / "queries" / "exonerate_refine_ids.txt").write_text("p1\n")
    with pytest.raises(StageError):
        s6.run(ctx)


# ---------------------------------------------------------------------------
# Stage 7: end-to-end run (pure pandas)
# ---------------------------------------------------------------------------

def _stage7_ctx(tmp_path: Path, dry_run: bool = False) -> RunContext:
    ctx = _make_ctx(tmp_path, dry_run=dry_run)
    s3_tables = ctx.paths.stage_dir(3) / "tables"
    s3_tables.mkdir(parents=True, exist_ok=True)
    runtime_tiered = pd.read_csv(FIX / "s07_tiered.tsv", sep="\t")
    runtime_tiered["tier"] = runtime_tiered["tier"].replace({"TIER_3_ARCHITECTURAL_VARIANT": "TIER_3B_ARCHITECTURAL_VARIANT"})
    runtime_tiered.to_csv(s3_tables / "tiered_candidates.tsv", sep="\t", index=False)
    s5_tables = ctx.paths.stage_dir(5) / "tables"
    s5_tables.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIX / "s07_miniprot_summary.tsv", s5_tables / "miniprot_summary.tsv")
    s6_gff = ctx.paths.stage_dir(6) / "gff"
    s6_gff.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIX / "s07_exonerate_merged.gff3", s6_gff / "exonerate_merged.gff3")
    # No exonerate_summary.tsv on purpose: exercises the GFF fallback path.
    return ctx


def test_stage7_full_run(tmp_path):
    ctx = _stage7_ctx(tmp_path)
    outs = s7.run(ctx)
    stage7 = ctx.paths.stage_dir(7)
    final = ctx.paths.final_dir

    report7 = stage7 / "tables" / "nlr_final_report.tsv"
    report_f = final / "nlr_final_report.tsv"
    assert report7.read_text() == report_f.read_text()
    report = pd.read_csv(report7, sep="\t")
    assert len(report) == 5
    assert "custom_note" not in report.columns
    # Exonerate loci came from the GFF fallback.
    assert report.set_index("protein_id").loc["p3", "exo_start"] == 40

    strict = pd.read_csv(stage7 / "tables" / "nlr_strict_candidates.tsv", sep="\t")
    assert set(strict["protein_id"]) == {"p1", "p2", "p3", "p5"}  # updated runtime retains NON_NLR_ANNOT
    assert (final / "nlr_strict_candidates.tsv").read_text() == (
        stage7 / "tables" / "nlr_strict_candidates.tsv"
    ).read_text()

    for name in ["tier_summary.tsv", "qc_summary.tsv"]:
        assert (stage7 / "summaries" / name).read_text() == (final / name).read_text()
    qc = pd.read_csv(final / "qc_summary.tsv", sep="\t")
    metrics = dict(zip(qc["metric"], qc["value"]))
    assert metrics["n_total_candidates"] == 5
    assert metrics["n_strict_candidates_combined"] == 4

    for name in ["nlr_candidates.bed", "nlr_miniprot_loci.bed", "nlr_exonerate_loci.bed"]:
        assert (stage7 / "beds" / name).read_text() == (final / name).read_text()
    assert (final / "nlr_candidates.bed").read_text().splitlines()[0] == (
        "ctg1\t99\t900\tp1|TIER_1A_HIGH_CONFIDENCE\t0\t+"
    )

    manifest = (stage7 / "stage7_manifest.tsv").read_text().splitlines()
    assert manifest[0] == "label\tpath"
    labels = [ln.split("\t")[0] for ln in manifest[1:]]
    assert "tables/nlr_final_report.tsv" in labels
    assert "tables/nlr_strict_candidates.comprehensive.tsv" in labels
    assert (stage7 / "metadata" / "run_info.json").is_file()
    assert report7 in outs and report_f in outs


def test_stage7_missing_tiered_raises(tmp_path):
    ctx = _stage7_ctx(tmp_path)
    (ctx.paths.stage_dir(3) / "tables" / "tiered_candidates.tsv").unlink()
    with pytest.raises(StageError):
        s7.run(ctx)


def test_stage7_dry_run_missing_upstream_returns_empty(tmp_path):
    ctx = _stage7_ctx(tmp_path, dry_run=True)
    (ctx.paths.stage_dir(3) / "tables" / "tiered_candidates.tsv").unlink()
    assert s7.run(ctx) == []

# Portable execution corrections shared with the first prototype.

def test_stage5_zero_byte_miniprot_sends_all_queries_to_exonerate(tmp_path):
    miniprot = _make_shim(
        tmp_path / "bin", "miniprot",
        'if [ "$1" = "--version" ]; then echo miniprot-shim; fi\nexit 0',
    )
    ctx = _make_ctx(tmp_path, tools={"miniprot": str(miniprot)})
    _copy_tier_fixtures(ctx.paths.stage_dir(4))
    outputs = s5.run(ctx)
    stage5 = ctx.paths.stage_dir(5)
    assert (stage5 / "queries/exonerate_refine_ids.txt").read_text().splitlines() == ["p1", "p2", "p3", "p4"]
    for name in ("miniprot_features.tsv", "miniprot_summary.tsv", "miniprot_hit_counts.tsv"):
        assert pd.read_csv(stage5 / "tables" / name, sep="\t").empty
    assert all(path.is_file() for path in outputs)


def test_no_rescue_stage5_to_stage6_finishes_with_schemas(tmp_path):
    tools = {name: str(_version_shim(tmp_path / "bin", name))
             for name in ("miniprot", "seqkit", "exonerate")}
    ctx = _make_ctx(tmp_path, tools=tools)
    s5.run(ctx)
    # Empty Stage5 dedup FASTA and absent Stage0 proteins are fine with no work.
    outputs = s6.run(ctx)
    assert all(path.is_file() for path in outputs)
    table = pd.read_csv(ctx.paths.stage_dir(6) / "tables/exonerate_summary.tsv", sep="\t")
    assert table.empty and table.columns.tolist() == s6.SUMMARY_COLUMNS


def _stage6_runtime_ctx(tmp_path, failing_tool=None):
    """A literal-ID SeqKit stand-in, unlike the original regex-only shim."""
    import sys
    ctx = _stage6_ctx(tmp_path)
    queries = ctx.paths.stage_dir(5) / "queries"
    (queries / "exonerate_refine_ids.txt").write_text("p1\n")
    seqkit = tmp_path / "bin/seqkit"
    seqkit.write_text(f"#!{sys.executable}\n" + '''
import sys
from pathlib import Path
args = sys.argv[1:]
if "--version" in args:
    print("seqkit literal-ID shim")
    raise SystemExit(0)
''' + ("raise SystemExit(9)\n" if failing_tool == "seqkit" else "") + '''
if "-f" in args:
    ids = set(Path(args[args.index("-f") + 1]).read_text().splitlines())
else:
    ids = {args[args.index("-p") + 1]}
for entry in Path(args[-1]).read_text().split(">")[1:]:
    if entry.split()[0] in ids:
        print(">" + entry, end="")
''')
    seqkit.chmod(0o755)
    exonerate = tmp_path / "bin/exonerate"
    exonerate.write_text(f"#!{sys.executable}\n" + '''
import sys
if "--version" in sys.argv:
    print("exonerate shim")
    raise SystemExit(0)
''' + ("raise SystemExit(9)\n" if failing_tool == "exonerate" else "") + '''
print("C4 Alignment:")
print("ctg1\\texonerate\\tgene\\t100\\t900\\t580\\t+\\t.\\tgene_id 0 ; sequence p1 ;")
print("ctg1\\texonerate\\texon\\t100\\t300\\t.\\t+\\t.\\tgene_id 0 ;")
print("ctg1\\texonerate\\texon\\tinvalid\\t300\\t.\\t+\\t.\\tgene_id 0 ;")
print("vulgar: alignment prose")
''')
    exonerate.chmod(0o755)
    return ctx


def test_stage6_exact_extraction_known_identity_feature_rows_only_and_no_stale_merge(tmp_path):
    ctx = _stage6_runtime_ctx(tmp_path)
    stale_dir = ctx.paths.stage_dir(6) / "gff/per_query"
    stale_dir.mkdir(parents=True)
    (stale_dir / "stale.exonerate.gff").write_text(
        "ctg9\texonerate\tgene\t1\t9\t.\t+\t.\tsequence stale ;\n")
    outputs = s6.run(ctx)
    stage6 = ctx.paths.stage_dir(6)
    merged = stage6 / "gff/exonerate_merged.gff3"
    assert len(merged.read_text().splitlines()) == 3
    assert "FUNLR_query=p1;LegacyAttrs=" in merged.read_text()
    assert "C4 Alignment:" not in merged.read_text()
    assert "stale" not in merged.read_text()
    raw = stage6 / "gff/per_query/p1.exonerate.gff"
    assert "C4 Alignment:" in raw.read_text()
    assert raw in outputs
    summary = pd.read_csv(stage6 / "tables/exonerate_summary.tsv", sep="\t")
    assert summary["protein_id"].tolist() == ["p1"]
    assert summary["n_features"].tolist() == [2]
    fallback = s7.load_exonerate_summary(tmp_path / "missing.tsv", merged)
    assert fallback["protein_id"].tolist() == ["p1"]


@pytest.mark.parametrize("tool", ["exonerate"])
def test_stage6_external_tool_failures_are_fatal(tmp_path, tool):
    import subprocess
    ctx = _stage6_runtime_ctx(tmp_path, failing_tool=tool)
    with pytest.raises(subprocess.CalledProcessError):
        s6.run(ctx)


def test_stage6_missing_extracted_ids_are_fatal(tmp_path):
    ctx = _stage6_runtime_ctx(tmp_path)
    (ctx.paths.stage_dir(5) / "queries/exonerate_refine_ids.txt").write_text("p1\nabsent\n")
    with pytest.raises(StageError, match="absent"):
        s6.run(ctx)


def test_merge_gffs_roundtrips_known_query_id(tmp_path):
    raw = tmp_path / "raw.gff"
    raw.write_text("ctg1\texonerate\texon\t1\t9\t.\t+\t.\tgene_id 0 ;\n")
    merged = tmp_path / "merged.gff3"
    pid = "query:with|punctuation%"
    s6.merge_gffs([(pid, raw)], merged)
    assert read_exonerate_features(merged)["protein_id"].tolist() == [pid]
    assert s7.load_exonerate_summary(tmp_path / "missing.tsv", merged)["protein_id"].tolist() == [pid]


def test_stage7_empty_report_preserves_strict_and_final_schemas(tmp_path):
    ctx = _stage7_ctx(tmp_path)
    tiered = ctx.paths.stage_dir(3) / "tables/tiered_candidates.tsv"
    frame = pd.read_csv(tiered, sep="\t").iloc[:0]
    frame.to_csv(tiered, sep="\t", index=False)
    s7.run(ctx)
    final = pd.read_csv(ctx.paths.final_dir / "nlr_final_report.tsv", sep="\t")
    strict = pd.read_csv(ctx.paths.final_dir / "nlr_strict_candidates.tsv", sep="\t")
    assert strict.empty and final.empty
    assert strict.columns.tolist() == final.columns.tolist()
    assert "protein_id" in strict.columns
    assert (ctx.paths.final_dir / "nlr_candidates.bed").read_text() == ""
