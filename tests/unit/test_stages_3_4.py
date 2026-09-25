"""Unit tests for stages 3 (tier) and 4 (export).

Pure logic (flags, tiers, rescue score, repeat proximity, BED writing, ID
extraction) is tested without a RunContext; stage run() paths are tested with
fixture inputs and a seqkit shim executable.
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from funlr.config import FunlrConfig
from funlr.core.context import RunContext, RunPaths
from funlr.core.errors import StageError
from funlr.core.runner import CommandRunner
from funlr.stages import stage3_tier as s3
from funlr.stages import stage4_export as s4

FIX = Path(__file__).resolve().parents[1] / "fixtures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ctx(tmp_path: Path, tools: dict[str, str] | None = None,
              dry_run: bool = False) -> RunContext:
    out = tmp_path / "out"
    cfg = FunlrConfig(
        {
            "sample": {"sample_id": "testsample", "species_id": "Test_sp"},
            "execution": {"output_dir": str(out), "cpu_threads": 2},
            "tools": tools or {},
        }
    )
    logger = logging.getLogger("funlr.test")
    paths = RunPaths(output_dir=out)
    runner = CommandRunner(tools or {}, out / "logs" / "commands.jsonl", logger,
                           dry_run=dry_run)
    return RunContext(config=cfg, paths=paths, runner=runner, state=None,
                      logger=logger, dry_run=dry_run)


def _seqkit_shim(bin_dir: Path) -> Path:
    """Shim seqkit: handles --version and `grep -f ids prot` (single-line FASTA)."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    p = bin_dir / "seqkit"
    p.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "--version" ] || [ "$1" = "version" ]; then\n'
        '  echo "seqkit-shim 0.1"; exit 0\n'
        "fi\n"
        '# argv: grep -f <ids> <prot>\n'
        'awk \'NR==FNR{ids[$1]=1;next} /^>/{keep=(substr($1,2) in ids)} keep\' "$3" "$4"\n'
    )
    p.chmod(0o755)
    return p


def _stage3_inputs(ctx: RunContext) -> None:
    """Copy the s03 fixtures into the stage0/stage2 layout of a ctx."""
    stage0 = ctx.paths.results_dir / "stage0"
    stage2 = ctx.paths.results_dir / "stage2"
    stage0.mkdir(parents=True, exist_ok=True)
    stage2.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIX / "s03_master.tsv", stage0 / "master_table.tsv")
    shutil.copy(FIX / "s03_contig_lengths.tsv", stage0 / "contig_lengths.tsv")
    arch = pd.read_csv(FIX / "s03_arch.tsv", sep="\t")
    arch["has_sensor"] = arch["has_repeat"]
    arch["has_ssfr_sensor"] = arch["has_repeat"]
    arch["has_non_ssfr_sensor"] = False
    arch["nbd_confidence"] = arch["has_nbd"].map({True: "HIGH", False: "LOW"})
    arch.to_csv(stage2 / "architecture_summary.tsv", sep="\t", index=False)
    for directory, names in ((stage0, ["input_manifest.tsv", "software_manifest.tsv"]),
                             (stage2, ["parse_filter_report.tsv", "run_info.txt"])):
        for name in names:
            (directory / name).write_text("fixture_metadata\n")


def _row(**kw) -> dict:
    return kw


# ---------------------------------------------------------------------------
# Stage 3: contig-end / rescue-flank / length flags
# ---------------------------------------------------------------------------


def test_contig_end_flag_both_ends():
    r = _row(contig_len=100000, start=5000, end=8000)
    assert s3.contig_end_flag(r, 10000) is True  # start within END
    r = _row(contig_len=100000, start=40000, end=95000)
    assert s3.contig_end_flag(r, 10000) is True  # contig_len - end within END
    r = _row(contig_len=100000, start=40000, end=45000)
    assert s3.contig_end_flag(r, 10000) is False
    # boundary: start exactly at END counts; END+1 does not
    assert s3.contig_end_flag(_row(contig_len=100000, start=10000, end=20000), 10000)
    assert not s3.contig_end_flag(_row(contig_len=100000, start=10001, end=20000), 10000)
    # missing coords -> False
    assert s3.contig_end_flag(_row(contig_len=None, start=1, end=2), 10000) is False


def test_rescue_possible_vs_no_flank_boundaries():
    # left_flank = start - 1; >= MIN_RESCUE_FLANK_BP means rescue possible.
    assert s3.rescue_possible_flag(_row(contig_len=100000, start=5001, end=95000), 5000)
    assert not s3.rescue_possible_flag(_row(contig_len=6000, start=5000, end=5500), 5000)
    # right flank exactly at the threshold also counts
    assert s3.rescue_possible_flag(_row(contig_len=10000, start=1000, end=5000), 5000)
    # any NA coord -> False
    assert s3.rescue_possible_flag(_row(contig_len=10000, start=None, end=5), 5000) is False


def test_length_flags():
    assert s3.length_flags(_row(protein_length=499)) == ["SHORT_PROTEIN"]
    assert s3.length_flags(_row(protein_length=500)) == []
    assert s3.length_flags(_row(protein_length=3000)) == []
    assert s3.length_flags(_row(protein_length=3001)) == ["VERY_LONG_PROTEIN"]
    assert s3.length_flags(_row(protein_length=None)) == []
    assert s3.length_flags(_row(protein_length="junk")) == []


# ---------------------------------------------------------------------------
# Stage 3: NBD_ONLY strictness + the empty-domain NaN quirk
# ---------------------------------------------------------------------------


def test_nbd_only_strictness():
    base = dict(has_nbd=True, has_repeat=False, order_ok=True,
                contig_len=100000, start=40000, end=46000, protein_length=900)
    # Only NBD tokens -> NBD_ONLY
    r = dict(base, domains_grouped="NACHT", domains_raw="NACHT")
    assert "NBD_ONLY" in s3.illumina_arch_flags(r, 10000, 5000)
    # Other domain besides NBD -> NOT NBD_ONLY
    r = dict(base, domains_grouped="NACHT;Fbox", domains_raw="NACHT;Fbox")
    assert "NBD_ONLY" not in s3.illumina_arch_flags(r, 10000, 5000)
    # Repeat present -> NOT NBD_ONLY
    r = dict(base, has_repeat=True, domains_grouped="NACHT;WD40")
    assert "NBD_ONLY" not in s3.illumina_arch_flags(r, 10000, 5000)


def test_nbd_only_nan_quirk():
    # Legacy quirk: empty domains_* read back from TSV are NaN; str(NaN) ->
    # token "nan" counts as an "other domain", so zero-domain-hit NBD
    # candidates are NOT NBD_ONLY.
    r = _row(domains_grouped=float("nan"), domains_raw=float("nan"))
    assert s3.has_other_domains_besides_nbd(r) is True
    # Truly empty strings (in-memory) -> no other domains
    assert s3.has_other_domains_besides_nbd(_row(domains_grouped="", domains_raw="")) is False
    # falls back to domains_raw when domains_grouped is empty
    assert s3.has_other_domains_besides_nbd(_row(domains_grouped="", domains_raw="NACHT")) is False


# ---------------------------------------------------------------------------
# Stage 3: tier truth table
# ---------------------------------------------------------------------------

TIER_CASES = [
    # (name, row, expected tier)
    ("1A", _row(has_nbd=True, has_repeat=True, flags="PASS"),
     "TIER_1A_HIGH_CONFIDENCE"),
    ("1B via CONTIG_END", _row(has_nbd=True, has_repeat=True, flags="CONTIG_END;RESCUE_POSSIBLE"),
     "TIER_1B_NEEDS_REVIEW"),
    ("1B via SHORT_PROTEIN", _row(has_nbd=True, has_repeat=True, flags="SHORT_PROTEIN"),
     "TIER_1B_NEEDS_REVIEW"),
    ("2A", _row(has_nbd=True, has_repeat=False, flags="NBD_ONLY;CONTIG_END;RESCUE_POSSIBLE"),
     "TIER_2A_HIGH_PRIORITY_RESCUE"),
    ("2B via RESCUE_NO_FLANK",
     _row(has_nbd=True, has_repeat=False, flags="NBD_ONLY;CONTIG_END;RESCUE_NO_FLANK"),
     "TIER_2B_RESCUE_CANDIDATE"),
    ("2B via REPEAT_NEARBY", _row(has_nbd=True, has_repeat=False, flags="NBD_ONLY;REPEAT_NEARBY"),
     "TIER_2B_RESCUE_CANDIDATE"),
    ("2B via NBD_ONLY", _row(has_nbd=True, has_repeat=False, flags="NBD_ONLY"),
     "TIER_2B_RESCUE_CANDIDATE"),
    ("3", _row(has_nbd=True, has_repeat=False, flags="PASS"),
     "TIER_3_ARCHITECTURAL_VARIANT"),
    ("FP domain + repeat -> 3", _row(has_nbd=True, has_repeat=True, flags="FP_DOMAIN_PRESENT"),
     "TIER_3_ARCHITECTURAL_VARIANT"),
    ("4A", _row(has_nbd=False, has_repeat=True, flags="PASS"),
     "TIER_4A_REPEAT_ONLY_NO_NBD"),
    ("4B", _row(has_nbd=True, has_repeat=False, flags="NON_NLR_ANNOT"),
     "TIER_4B_LIKELY_STAND_HOUSEKEEPING"),
    ("4B precedence over NBD+repeat",
     _row(has_nbd=True, has_repeat=True, flags="NON_NLR_ANNOT"),
     "TIER_4B_LIKELY_STAND_HOUSEKEEPING"),
    ("4C", _row(has_nbd=False, has_repeat=False, flags="PASS"),
     "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL"),
]


@pytest.mark.parametrize("name,row,expected", TIER_CASES,
                         ids=[c[0] for c in TIER_CASES])
def test_assign_tier_truth_table(name, row, expected):
    assert s3.assign_tier(row) == expected


# ---------------------------------------------------------------------------
# Stage 3: rescue_score terms and clamps
# ---------------------------------------------------------------------------

SCORE_CASES = [
    # base scores
    (_row(tier="TIER_2A_HIGH_PRIORITY_RESCUE", flags="PASS"), 90),
    (_row(tier="TIER_2B_RESCUE_CANDIDATE", flags="PASS"), 70),
    (_row(tier="TIER_1B_NEEDS_REVIEW", flags="PASS"), 60),
    (_row(tier="TIER_1A_HIGH_CONFIDENCE", flags="PASS"), 0),
    # individual terms
    (_row(tier="TIER_1A_HIGH_CONFIDENCE", flags="CONTIG_END"), 15),
    (_row(tier="TIER_1A_HIGH_CONFIDENCE", flags="RESCUE_POSSIBLE"), 10),
    (_row(tier="TIER_2B_RESCUE_CANDIDATE", flags="RESCUE_NO_FLANK"), 40),
    (_row(tier="TIER_1A_HIGH_CONFIDENCE", flags="REPEAT_NEARBY"), 10),
    (_row(tier="TIER_1A_HIGH_CONFIDENCE", flags="SHORT_PROTEIN"), 10),
    (_row(tier="TIER_2B_RESCUE_CANDIDATE", flags="FP_DOMAIN_PRESENT"), 45),
    (_row(tier="TIER_1B_NEEDS_REVIEW", flags="NON_NLR_ANNOT"), 25),
    (_row(tier="TIER_1B_NEEDS_REVIEW", flags="VERY_LONG_PROTEIN"), 50),
    # combined
    (_row(tier="TIER_1B_NEEDS_REVIEW", flags="CONTIG_END;RESCUE_POSSIBLE"), 85),
    # upper clamp: 90 + 15 + 10 = 115 -> 100
    (_row(tier="TIER_2A_HIGH_PRIORITY_RESCUE", flags="CONTIG_END;RESCUE_POSSIBLE"), 100),
    # lower clamp: 70 + 10 + 10 - 30 - 25 - 35 - 10 = -10 -> 0
    (_row(tier="TIER_2B_RESCUE_CANDIDATE",
          flags="REPEAT_NEARBY;SHORT_PROTEIN;RESCUE_NO_FLANK;FP_DOMAIN_PRESENT;"
                "NON_NLR_ANNOT;VERY_LONG_PROTEIN"), 0),
]


@pytest.mark.parametrize("row,expected", SCORE_CASES)
def test_rescue_score_terms_and_clamps(row, expected):
    assert s3.rescue_score(row) == expected


# ---------------------------------------------------------------------------
# Stage 3: repeat proximity heuristic
# ---------------------------------------------------------------------------


def _prox_df(rows):
    return pd.DataFrame(rows, columns=["scaffold", "has_nbd", "has_repeat",
                                       "_start_i", "_end_i"])


def test_repeat_nearby_window_edges():
    # candidate at 20000-21000; window 10000-31000 with repeat_prox_bp=10000
    rows = [
        ("scf1", True, False, 20000, 21000),   # candidate
        ("scf1", False, True, 9000, 10000),    # repeat ends exactly at win_s -> nearby
    ]
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [True, False]

    rows[1] = ("scf1", False, True, 9000, 9999)  # one bp outside -> not nearby
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [False, False]

    rows[1] = ("scf1", False, True, 31000, 32000)  # starts exactly at win_e -> nearby
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [True, False]

    rows[1] = ("scf1", False, True, 31001, 32000)  # one bp outside
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [False, False]


def test_repeat_nearby_neighbor_count_edge():
    # neighbors=1: repeat two positions away is out of range even within window
    rows = [
        ("scf1", True, False, 20000, 21000),   # candidate (j=1 after sorting)
        ("scf1", False, False, 15000, 16000),  # filler non-repeat
        ("scf1", False, True, 25000, 26000),   # repeat within bp window
    ]
    df = _prox_df(rows).sort_values("_start_i")  # repeat, filler... reorder by start
    res = s3.compute_repeat_nearby(df, 10000, 1)
    # sorted: filler(15000), candidate(20000), repeat(25000): repeat is j+1 -> nearby
    assert bool(res[df.index[1]]) is True

    rows = [
        ("scf1", True, False, 20000, 21000),
        ("scf1", False, False, 15000, 16000),
        ("scf1", False, False, 22000, 23000),  # second filler pushes repeat to j+2
        ("scf1", False, True, 25000, 26000),
    ]
    df = _prox_df(rows).sort_values("_start_i")
    res = s3.compute_repeat_nearby(df, 10000, 1)
    cand_idx = df.index[df["_start_i"] == 20000][0]
    assert bool(res[cand_idx]) is False


def test_repeat_nearby_cross_scaffold_and_candidate_filters():
    rows = [
        ("scf1", True, False, 20000, 21000),   # candidate on scf1
        ("scf2", False, True, 20050, 21050),   # repeat on scf2, overlapping coords
    ]
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [False, False]

    # candidates that carry repeats themselves are never scanned
    rows = [
        ("scf1", True, True, 20000, 21000),    # has_repeat -> not a candidate
        ("scf1", False, True, 25000, 26000),
    ]
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [False, False]

    # no-NBD proteins are never candidates
    rows = [
        ("scf1", False, False, 20000, 21000),
        ("scf1", False, True, 25000, 26000),
    ]
    res = s3.compute_repeat_nearby(_prox_df(rows), 10000, 10)
    assert res.tolist() == [False, False]


# ---------------------------------------------------------------------------
# Stage 3: non-NLR annotation regex tolerance
# ---------------------------------------------------------------------------


def test_bad_regex_tolerance():
    warnings: list[str] = []
    assert s3._compile_or_none("([", warnings.append) is None
    assert warnings and "bad regex pattern ignored" in warnings[0]
    assert s3._compile_or_none("") is None
    assert s3._compile_or_none("helicase") is not None


def test_non_nlr_annot_flag_sources():
    re_kw = s3._compile_or_none("replication")
    re_dom = s3._compile_or_none("Pkinase")
    re_go = s3._compile_or_none("chromosome")
    assert s3.non_nlr_annot_flag(
        _row(Description="DNA replication factor", Preferred_name="x", GOs="", domains_raw=""),
        re_kw, re_go, re_dom) is True
    assert s3.non_nlr_annot_flag(
        _row(Description="", Preferred_name="", GOs="", domains_raw="NACHT;Pkinase"),
        re_kw, re_go, re_dom) is True
    assert s3.non_nlr_annot_flag(
        _row(Description="", Preferred_name="", GOs="GO:chromosome", domains_raw=""),
        re_kw, re_go, re_dom) is True
    assert s3.non_nlr_annot_flag(
        _row(Description="hypothetical", Preferred_name="x", GOs="", domains_raw="NACHT"),
        re_kw, re_go, re_dom) is False
    # None regexes never hit
    assert s3.non_nlr_annot_flag(_row(Description="replication"), None, None, None) is False


# ---------------------------------------------------------------------------
# Stage 3: BED writing
# ---------------------------------------------------------------------------


def test_write_bed_format(tmp_path):
    frame = pd.DataFrame(
        [
            {"scaffold": "scfB", "start": 5000, "end": 8000, "strand": "+",
             "protein_id": "t1b_end", "tier": "TIER_1B_NEEDS_REVIEW",
             "rescue_priority": 85},
            {"scaffold": "scfA", "start": 1, "end": 900, "strand": "-",
             "protein_id": "p2", "tier": "TIER_2A_HIGH_PRIORITY_RESCUE",
             "rescue_priority": 100},
        ]
    )
    out = tmp_path / "x.bed"
    s3.write_bed(out, frame)
    lines = out.read_text().splitlines()
    # 0-based start, name=protein_id|tier, score=int(rescue_priority*10)
    assert lines[0] == "scfB\t4999\t8000\tt1b_end|TIER_1B_NEEDS_REVIEW\t850\t+"
    assert lines[1] == "scfA\t0\t900\tp2|TIER_2A_HIGH_PRIORITY_RESCUE\t1000\t-"


# ---------------------------------------------------------------------------
# Stage 3: full table build on the fixture set (every tier reached)
# ---------------------------------------------------------------------------

EXPECTED_TIERS = {
    "t1a": "TIER_1A_HIGH_CONFIDENCE",
    "t1b_end": "TIER_1B_NEEDS_REVIEW",
    "t1b_short": "TIER_1B_NEEDS_REVIEW",
    "t2a": "TIER_2A_HIGH_PRIORITY_RESCUE",
    "t2b_noflank": "TIER_2B_RESCUE_CANDIDATE",
    "t2b_nbdonly": "TIER_2B_RESCUE_CANDIDATE",
    "t2b_repeatnear": "TIER_2B_RESCUE_CANDIDATE",
    "rep1": "TIER_4A_REPEAT_ONLY_NO_NBD",
    "t3": "TIER_3_ARCHITECTURAL_VARIANT",
    "t3_fp": "TIER_3_ARCHITECTURAL_VARIANT",
    "t3_nan": "TIER_3_ARCHITECTURAL_VARIANT",
    "t4a": "TIER_4A_REPEAT_ONLY_NO_NBD",
    "t4b": "TIER_4B_LIKELY_STAND_HOUSEKEEPING",
    "t4c": "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL",
    "t_orphan": "TIER_1A_HIGH_CONFIDENCE",
}

DEFAULT_FP = FunlrConfig()["false_positive"]


def _build_fixture_table():
    master = pd.read_csv(FIX / "s03_master.tsv", sep="\t")
    arch = pd.read_csv(FIX / "s03_arch.tsv", sep="\t")
    clen = pd.read_csv(FIX / "s03_contig_lengths.tsv", sep="\t", header=None,
                       names=["scaffold", "contig_len"])
    return s3.build_tiered_table(
        master, arch, clen,
        contig_end_bp=10000, min_rescue_flank_bp=5000,
        repeat_prox_bp=10000, repeat_prox_neighbors=10,
        non_nlr_keywords_regex=DEFAULT_FP["non_nlr_keywords_regex"],
        non_nlr_domain_regex=DEFAULT_FP["non_nlr_domain_regex"],
        non_nlr_go_regex=DEFAULT_FP["non_nlr_go_regex"],
    )


def test_build_tiered_table_truth_table():
    df = _build_fixture_table().set_index("protein_id")
    for pid, tier in EXPECTED_TIERS.items():
        assert df.loc[pid, "tier"] == tier, pid
    # spot-check flags and scores
    assert df.loc["t1b_end", "flags"] == "CONTIG_END;RESCUE_POSSIBLE"
    assert df.loc["t1b_end", "rescue_priority"] == 85
    assert df.loc["t2a", "rescue_priority"] == 100  # clamped from 115
    assert df.loc["t2b_noflank", "flags"] == "NBD_ONLY;CONTIG_END;RESCUE_NO_FLANK"
    assert df.loc["t2b_repeatnear", "flags"] == "NBD_ONLY;REPEAT_NEARBY"
    assert df.loc["t1a", "flags"] == "PASS"
    # NaN-domain quirk through the real read path: t3_nan is not NBD_ONLY
    assert "NBD_ONLY" not in df.loc["t3_nan", "flags"]
    # every tier name is exercised
    assert set(EXPECTED_TIERS.values()) == set(df["tier"].unique())


def test_build_tiered_table_missing_columns():
    master = pd.read_csv(FIX / "s03_master.tsv", sep="\t")
    arch = pd.read_csv(FIX / "s03_arch.tsv", sep="\t")
    clen = pd.read_csv(FIX / "s03_contig_lengths.tsv", sep="\t", header=None,
                       names=["scaffold", "contig_len"])
    kw = dict(contig_end_bp=10000, min_rescue_flank_bp=5000,
              repeat_prox_bp=10000, repeat_prox_neighbors=10)
    with pytest.raises(StageError, match="architecture_summary.tsv missing columns"):
        s3.build_tiered_table(master, arch.drop(columns=["has_nbd"]), clen, **kw)
    with pytest.raises(StageError, match="master_table.tsv missing columns"):
        s3.build_tiered_table(master.drop(columns=["scaffold"]), arch, clen, **kw)


def test_build_tiered_table_bad_regex_runs():
    master = pd.read_csv(FIX / "s03_master.tsv", sep="\t")
    arch = pd.read_csv(FIX / "s03_arch.tsv", sep="\t")
    clen = pd.read_csv(FIX / "s03_contig_lengths.tsv", sep="\t", header=None,
                       names=["scaffold", "contig_len"])
    df = s3.build_tiered_table(
        master, arch, clen, contig_end_bp=10000, min_rescue_flank_bp=5000,
        repeat_prox_bp=10000, repeat_prox_neighbors=10,
        non_nlr_keywords_regex="([",  # bad regex: warn-and-ignore
        non_nlr_domain_regex="", non_nlr_go_regex="",
    )
    # with the keywords regex ignored, t4b loses NON_NLR_ANNOT and falls to 1A
    assert df.set_index("protein_id").loc["t4b", "tier"] == "TIER_1A_HIGH_CONFIDENCE"


def test_build_tiered_table_one_to_one_merge_error():
    master = pd.read_csv(FIX / "s03_master.tsv", sep="\t")
    arch = pd.read_csv(FIX / "s03_arch.tsv", sep="\t")
    clen = pd.read_csv(FIX / "s03_contig_lengths.tsv", sep="\t", header=None,
                       names=["scaffold", "contig_len"])
    dup_master = pd.concat([master, master.iloc[[0]]], ignore_index=True)
    with pytest.raises(StageError, match="one_to_one"):
        s3.build_tiered_table(dup_master, arch, clen, contig_end_bp=10000,
                              min_rescue_flank_bp=5000, repeat_prox_bp=10000,
                              repeat_prox_neighbors=10)


# ---------------------------------------------------------------------------
# Stage 3: run() — outputs, dry-run, missing inputs
# ---------------------------------------------------------------------------


def test_stage3_run_full(tmp_path):
    ctx = _make_ctx(tmp_path)
    _stage3_inputs(ctx)
    outputs = s3.run(ctx)
    stage3 = ctx.paths.stage_dir(3)
    tiered = stage3 / "tables" / "tiered_candidates.tsv"
    assert tiered in outputs
    df = pd.read_csv(tiered, sep="\t")
    assert not {"flags_list", "_start_i", "_end_i"}.intersection(df.columns)
    expected = {
        "t1a": "TIER_1A_HIGH_CONFIDENCE", "t1b_end": "TIER_1B_NEEDS_REVIEW",
        "t1b_short": "TIER_4F_VERY_SHORT_FRAGMENT", "t2a": "TIER_2A_HIGH_PRIORITY_RESCUE",
        "t2b_noflank": "TIER_2B_RESCUE_CANDIDATE", "t2b_nbdonly": "TIER_2B_RESCUE_CANDIDATE",
        "t2b_repeatnear": "TIER_2B_RESCUE_CANDIDATE", "rep1": "TIER_4A_REPEAT_ONLY_NO_NBD",
        "t3": "TIER_3B_ARCHITECTURAL_VARIANT", "t3_fp": "TIER_1A_HIGH_CONFIDENCE",
        "t3_nan": "TIER_2B_RESCUE_CANDIDATE", "t4a": "TIER_4A_REPEAT_ONLY_NO_NBD",
        "t4b": "TIER_1A_HIGH_CONFIDENCE", "t4c": "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL",
        "t_orphan": "TIER_1A_HIGH_CONFIDENCE",
    }
    assert df.set_index("protein_id")["tier"].to_dict() == expected
    assert list(df.columns) == ["protein_id"] + sorted(set(df.columns) - {"protein_id"})
    for tier in set(expected.values()):
        assert (stage3 / "tables" / f"{tier}.tsv").is_file()
    rescue = pd.read_csv(stage3 / "tables/rescue_priority.tsv", sep="\t")
    assert set(rescue.protein_id) == {"t1b_end", "t2a", "t2b_noflank", "t2b_nbdonly", "t2b_repeatnear", "t3_nan"}
    assert rescue.rescue_priority.tolist() == sorted(rescue.rescue_priority, reverse=True)
    bed = (stage3 / "beds/nlr_candidates.bed").read_text().splitlines()
    assert len(bed) == 14 and all("t_orphan" not in line for line in bed)
    assert (stage3 / "beds/by_tier/TIER_2A_HIGH_PRIORITY_RESCUE.bed").read_text().count("\n") == 1
    assert (stage3 / "qc/duplicate_protein_ids.summary.tsv").read_text() == "protein_id\tduplicate_rows\tunique_signatures\n"
    info = (stage3 / "metadata/run_info.txt").read_text()
    assert "stage: 3\n" in info and "genome: Test_sp\n" in info
    assert "nlr_profile: ILLUMINA\n" in info and "repeat_prox_neighbors: 10\n" in info


def test_stage3_run_missing_input_stageerror(tmp_path):
    ctx = _make_ctx(tmp_path)
    with pytest.raises(StageError, match="Required input missing/empty"):
        s3.run(ctx)


def test_stage3_dry_run_missing_inputs(tmp_path, caplog):
    ctx = _make_ctx(tmp_path, dry_run=True)
    with caplog.at_level(logging.INFO):
        assert s3.run(ctx) == []
    assert any("PLANNED" in r.message for r in caplog.records)


def test_stage3_dry_run_with_inputs_preserves_outputs(tmp_path):
    ctx = _make_ctx(tmp_path, dry_run=True)
    _stage3_inputs(ctx)
    outputs = s3.run(ctx)
    assert outputs == []
    assert not (ctx.paths.stage_dir(3) / "tables" / "tiered_candidates.tsv").exists()


# ---------------------------------------------------------------------------
# Stage 4: extract_ids / export_faa
# ---------------------------------------------------------------------------


def test_extract_ids_semantics(tmp_path):
    tsv = tmp_path / "in.tsv"
    tsv.write_text("protein_id\tother\np1\tx\np2\ty\n\t\np3\tz\n")
    out = tmp_path / "out.ids"
    s4.extract_ids(tsv, out)
    # tail -n +2 | cut -f1 | awk 'NF>0'
    assert out.read_text().splitlines() == ["p1", "p2", "p3"]
    # missing/empty input -> empty output
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    s4.extract_ids(empty, out)
    assert out.read_text() == ""
    s4.extract_ids(tmp_path / "nonexistent.tsv", out)
    assert out.read_text() == ""


def test_export_faa_with_seqkit_shim(tmp_path):
    shim = _seqkit_shim(tmp_path / "bin")
    ctx = _make_ctx(tmp_path, tools={"seqkit": str(shim)})
    ids = tmp_path / "x.ids"
    ids.write_text("t1a\nt4c\n")
    out = tmp_path / "x.faa"
    s4.export_faa(ctx, ids, FIX / "s04_proteins.faa", out)
    assert out.read_text() == ">t1a\n" + "M" + "A" * 29 + "\n>t4c\n" + "M" + "A" * 29 + "\n"


def test_export_faa_empty_ids_no_call(tmp_path):
    shim = _seqkit_shim(tmp_path / "bin")
    ctx = _make_ctx(tmp_path, tools={"seqkit": str(shim)})
    ids = tmp_path / "empty.ids"
    ids.write_text("")
    out = tmp_path / "x.faa"
    s4.export_faa(ctx, ids, FIX / "s04_proteins.faa", out)
    assert out.read_text() == ""  # pre-created empty
    # no seqkit invocation recorded
    log_path = ctx.paths.logs_dir / "commands.jsonl"
    assert not log_path.exists() or "seqkit" not in log_path.read_text()


def test_export_faa_failure_propagates(tmp_path):
    shim = tmp_path / "broken_seqkit"
    shim.write_text("#!/bin/sh\nexit 17\n")
    shim.chmod(0o755)
    ctx = _make_ctx(tmp_path, tools={"seqkit": str(shim)})
    ids = tmp_path / "queries.ids"
    ids.write_text("t1a\n")
    with pytest.raises(subprocess.CalledProcessError):
        s4.export_faa(ctx, ids, FIX / "s04_proteins.faa", tmp_path / "out.faa")


def test_empty_architecture_tiers_and_exports_finish(tmp_path):
    """A true zero-candidate discovery remains a valid, readable final handoff."""
    ctx = _stage4_ctx(tmp_path)
    _stage3_inputs(ctx)
    arch_path = ctx.paths.results_dir / "stage2" / "architecture_summary.tsv"
    pd.read_csv(arch_path, sep="\t").head(0).to_csv(arch_path, sep="\t", index=False)
    stage0 = ctx.paths.results_dir / "stage0"
    shutil.copy(FIX / "s04_proteins.faa", stage0 / "proteins_clean.faa")

    assert s3.run(ctx)
    tables = ctx.paths.results_dir / "stage3" / "tables"
    table = pd.read_csv(tables / "tiered_candidates.tsv", sep="\t")
    assert table.empty
    assert {"protein_id", "tier", "flags", "rescue_priority"}.issubset(table.columns)
    assert pd.read_csv(tables / "rescue_priority.tsv", sep="\t").empty
    assert (ctx.paths.results_dir / "stage3" / "beds" / "nlr_candidates.bed").read_text() == ""

    assert s4.run(ctx)
    stage4 = ctx.paths.stage_dir(4)
    assert len(list(stage4.glob("*.faa"))) == 2
    assert all(path.stat().st_size == 0 for path in stage4.glob("*.faa"))
    summary = pd.read_csv(stage4 / "fasta_export_summary.tsv", sep="\t")
    assert len(summary) == 2
    assert summary["seqs"].sum() == 0


# ---------------------------------------------------------------------------
# Stage 4: run() — manifest/summary, missing inputs, dry-run
# ---------------------------------------------------------------------------


def _stage4_ctx(tmp_path, dry_run=False):
    shim = _seqkit_shim(tmp_path / "bin")
    ctx = _make_ctx(tmp_path, tools={"seqkit": str(shim)}, dry_run=dry_run)
    stage0 = ctx.paths.stage_dir(0)
    for name in ("input_manifest.tsv", "software_manifest.tsv"):
        (stage0 / name).write_text("fixture_metadata\n")
    return ctx


def _make_tier_tsv(path: Path, ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("protein_id\ttier\n" + "".join(f"{i}\tT\n" for i in ids))


def test_stage4_run_exports(tmp_path):
    ctx = _stage4_ctx(tmp_path)
    stage0 = ctx.paths.stage_dir(0)
    shutil.copy(FIX / "s04_proteins.faa", stage0 / "proteins_clean.faa")
    tables = ctx.paths.results_dir / "stage3/tables"
    _make_tier_tsv(tables / "TIER_1A_HIGH_CONFIDENCE.tsv", ["t1a", "t_orphan"])
    _make_tier_tsv(tables / "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL.tsv", ["t4c"])
    _make_tier_tsv(tables / "tiered_candidates.tsv", ["t1a", "t4c", "t_orphan"])
    _make_tier_tsv(tables / "rescue_priority.tsv", [])
    outputs = s4.run(ctx)
    stage4 = ctx.paths.stage_dir(4)
    fasta = stage4 / "TIER_1A_HIGH_CONFIDENCE.faa"
    assert fasta in outputs and fasta.read_text().count(">") == 2
    assert (stage4 / "TIER_1A_HIGH_CONFIDENCE.ids").read_text().splitlines() == ["t1a", "t_orphan"]
    assert (stage4 / "tiered_candidates.faa").read_text().count(">") == 3
    assert (stage4 / "rescue_priority.faa").read_text() == ""
    assert not (stage4 / "TIER_2A_HIGH_PRIORITY_RESCUE.faa").exists()
    manifest = pd.read_csv(stage4 / "fasta_export_manifest.tsv", sep="\t")
    assert manifest.label.tolist() == ["TIER_1A_HIGH_CONFIDENCE", "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL", "tiered_candidates", "rescue_priority"]
    summary = pd.read_csv(stage4 / "fasta_export_summary.tsv", sep="\t").set_index("file")
    assert summary.loc["TIER_1A_HIGH_CONFIDENCE.faa", "seqs"] == 2
    assert summary.loc["tiered_candidates.faa", "seqs"] == 3
    assert summary.loc["rescue_priority.faa", "seqs"] == 0


def test_stage4_run_missing_prot_stageerror(tmp_path):
    ctx = _stage4_ctx(tmp_path)
    tables = ctx.paths.results_dir / "stage3" / "tables"
    _make_tier_tsv(tables / "tiered_candidates.tsv", ["t1a"])
    with pytest.raises(StageError, match="Missing/empty proteins FASTA"):
        s4.run(ctx)


def test_stage4_run_missing_tables_dir_stageerror(tmp_path):
    ctx = _stage4_ctx(tmp_path)
    stage0 = ctx.paths.results_dir / "stage0"
    stage0.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIX / "s04_proteins.faa", stage0 / "proteins_clean.faa")
    with pytest.raises(StageError, match="Stage 3 tables directory not found"):
        s4.run(ctx)


def test_stage4_dry_run_preserves_existing_outputs(tmp_path):
    ctx = _stage4_ctx(tmp_path, dry_run=True)
    tables = ctx.paths.results_dir / "stage3" / "tables"
    _make_tier_tsv(tables / "TIER_2A_HIGH_PRIORITY_RESCUE.tsv", ["t2a"])
    stage4 = ctx.paths.stage_dir(4)
    existing = stage4 / "tier2A_high_priority_rescue.faa"
    existing.write_text(">old\nAAAA\n")
    assert s4.run(ctx) == []
    assert existing.read_text() == ">old\nAAAA\n"
    assert not (stage4 / "fasta_export_manifest.tsv").exists()
    assert not (ctx.paths.logs_dir / "commands.jsonl").exists()


def test_stage4_dry_run_no_tables_returns_empty(tmp_path, caplog):
    ctx = _stage4_ctx(tmp_path, dry_run=True)
    with caplog.at_level(logging.INFO):
        assert s4.run(ctx) == []
    assert any("PLANNED" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# Integration: stage 3 -> stage 4 on the fixture set
# ---------------------------------------------------------------------------


def test_integration_stage3_then_stage4(tmp_path):
    ctx = _stage4_ctx(tmp_path)
    _stage3_inputs(ctx)
    stage0 = ctx.paths.results_dir / "stage0"
    shutil.copy(FIX / "s04_proteins.faa", stage0 / "proteins_clean.faa")

    out3 = s3.run(ctx)
    assert out3
    out4 = s4.run(ctx)
    stage4 = ctx.paths.stage_dir(4)

    def ids_in(faa: Path) -> set[str]:
        return {line[1:] for line in faa.read_text().splitlines()
                if line.startswith(">")}

    assert ids_in(stage4 / "TIER_2A_HIGH_PRIORITY_RESCUE.faa") == {"t2a"}
    assert ids_in(stage4 / "TIER_2B_RESCUE_CANDIDATE.faa") == {
        "t2b_noflank", "t2b_nbdonly", "t2b_repeatnear", "t3_nan"}
    assert ids_in(stage4 / "TIER_1B_NEEDS_REVIEW.faa") == {"t1b_end"}
    assert ids_in(stage4 / "TIER_1A_HIGH_CONFIDENCE.faa") == {"t1a", "t_orphan", "t3_fp", "t4b"}
    assert ids_in(stage4 / "TIER_3B_ARCHITECTURAL_VARIANT.faa") == {"t3"}
    assert ids_in(stage4 / "TIER_4A_REPEAT_ONLY_NO_NBD.faa") == {"rep1", "t4a"}
    assert ids_in(stage4 / "TIER_4F_VERY_SHORT_FRAGMENT.faa") == {"t1b_short"}
    assert ids_in(stage4 / "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL.faa") == {"t4c"}
    assert ids_in(stage4 / "tiered_candidates.faa") == set(EXPECTED_TIERS)
    assert ids_in(stage4 / "rescue_priority.faa") == {"t1b_end", "t2a", "t2b_noflank", "t2b_nbdonly", "t2b_repeatnear", "t3_nan"}
    manifest = pd.read_csv(stage4 / "fasta_export_manifest.tsv", sep="\t")
    assert manifest.label.tolist() == [
        "TIER_1A_HIGH_CONFIDENCE", "TIER_1B_NEEDS_REVIEW", "TIER_2A_HIGH_PRIORITY_RESCUE",
        "TIER_2B_RESCUE_CANDIDATE", "TIER_3B_ARCHITECTURAL_VARIANT", "TIER_4A_REPEAT_ONLY_NO_NBD",
        "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL", "TIER_4F_VERY_SHORT_FRAGMENT", "tiered_candidates", "rescue_priority"]
    assert (stage4 / "run_info.txt").is_file()
    assert stage4 / "run_info.txt" in out4
