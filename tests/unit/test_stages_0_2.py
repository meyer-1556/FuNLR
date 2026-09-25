"""Unit tests for stages 0-2 pure logic (legacy-faithful behavior)."""

from __future__ import annotations

from pathlib import Path
import logging
from types import SimpleNamespace

import pandas as pd
import pytest

from funlr.parsers.hmmer import parse_domtblout
from funlr.stages import stage0_validate as s0
from funlr.stages import stage1_discover as s1
from funlr.stages import stage2_architecture as s2
from funlr.config import FunlrConfig
from funlr.core.context import RunPaths

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"

EVAL_NBD = 1e-4
MIN_NBD_ALI = 150
EVAL_USE = 1e-3
MIN_DOM_ALI = 20


# ---------------------------------------------------------------------------
# Stage 0 helpers
# ---------------------------------------------------------------------------


def test_build_contig_lengths_two_columns_no_header(tmp_path):
    out = tmp_path / "contig_lengths.tsv"
    n = s0.build_contig_lengths(FIXTURES / "s00_genome.fai", out)
    assert n == 2
    assert out.read_text() == "scf1\t5000\nscf2\t3000\n"


def test_build_protein_files_clean_and_lengths(tmp_path):
    proteins = tmp_path / "proteins.fa"
    proteins.write_text(">p1 some description\nAAAA\nBBBB\n>p2\n" + "C" * 61 + "\n")
    faa = tmp_path / "proteins_clean.faa"
    lens = tmp_path / "protein_lengths.tsv"
    n = s0.build_protein_files(proteins, faa, lens)
    assert n == 2
    # IDs only in header, 60-column wrap
    assert faa.read_text() == f">p1\nAAAABBBB\n>p2\n{'C' * 60}\nC\n"
    assert lens.read_text() == "protein_id\tprotein_length\np1\t8\np2\t61\n"


def test_build_master_table_merges_coords_lens_eggnog(tmp_path):
    coords = tmp_path / "id_coords.tsv"
    coords.write_text(
        "protein_id\ttranscript_id\tgene_id\tscaffold\tstart\tend\tstrand\n"
        "p1\tp1\tg1\tscf1\t100\t900\t+\n"
        "p2\tp2\tg2\tscf1\t1000\t2000\t-\n"
    )
    lens = tmp_path / "protein_lengths.tsv"
    lens.write_text("protein_id\tprotein_length\np1\t300\np2\t250\n")
    eggnog = tmp_path / "eggnog.tsv"
    eggnog.write_text(
        "p1\tx\tx\tx\tx\tx\tS\tNACHT protein\tnachtA\tGO:1\tx\tK00001\n"
        "p2\tx\tx\tx\tx\tx\tA\tABC transporter\tabcA\tGO:2\tx\tK00002\n"
    )
    light = tmp_path / "eggnog_light.tsv"
    master_out = tmp_path / "master_table.tsv"
    master, eg = s0.build_master_table(coords, lens, eggnog, light, master_out)

    assert master.shape[0] == 2
    assert list(master.columns) == [
        "protein_id", "transcript_id", "gene_id", "scaffold", "start", "end",
        "strand", "protein_length", "COG_category", "Description",
        "Preferred_name", "GOs", "KEGG_ko",
    ]
    row1 = master.set_index("protein_id").loc["p1"]
    assert row1["protein_length"] == 300
    assert row1["Description"] == "NACHT protein"  # left merge keeps p1's eggnog row
    # eggnog light table written with header even when parsed
    assert light.read_text().splitlines()[0] == (
        "protein_id\tCOG_category\tDescription\tPreferred_name\tGOs\tKEGG_ko"
    )


def test_build_master_table_eggnog_missing_failsoft(tmp_path):
    coords = tmp_path / "id_coords.tsv"
    coords.write_text(
        "protein_id\ttranscript_id\tgene_id\tscaffold\tstart\tend\tstrand\n"
        "p1\tp1\tg1\tscf1\t100\t900\t+\n"
    )
    lens = tmp_path / "protein_lengths.tsv"
    lens.write_text("protein_id\tprotein_length\np1\t300\n")
    light = tmp_path / "eggnog_light.tsv"
    master_out = tmp_path / "master_table.tsv"
    master, eg = s0.build_master_table(coords, lens, None, light, master_out)
    assert eg.empty
    assert master.shape[0] == 1
    assert pd.isna(master.loc[0, "Description"])


def test_empty_architecture_has_complete_schema():
    hits = pd.DataFrame(columns=[
        "protein_id", "domain", "acc", "dom_i_evalue", "start", "end", "ali_len", "source",
    ])
    summary = s2.summarize_architecture([], hits, set())
    assert summary.empty
    assert list(summary.columns) == [
        "protein_id", "domains_raw", "domains_grouped", "has_nbd_pfam", "has_repeat",
        "order_ok", "flag_fp_domains", "flag_rare_lrr", "has_nbd_stage1", "has_nbd",
    ]


def test_stage2_empty_candidate_fasta_skips_scans(tmp_path, monkeypatch):
    from funlr.stages import _architecture_run
    paths = RunPaths(tmp_path / "out")
    stage1 = paths.stage_dir(1)
    (stage1 / "union_candidates.faa").write_text("")
    (stage1 / "nbd_candidate_ids_strict.txt").write_text("")
    (stage1 / "nbd_hits.tsv").write_text("protein_id\tnbd_hmm\tdom_i_evalue\tali_start\tali_end\tali_len\n")
    hmm = tmp_path / "profiles.hmm"
    hmm.write_text("already validated and prepared by the orchestrator\n")
    cfg = FunlrConfig({"databases": {"pfam": str(hmm), "custom_hmms": str(hmm)}, "asm": {"asm_enable": 0}})
    log = logging.getLogger("funlr.test.empty")
    ctx = SimpleNamespace(
        config=cfg, paths=paths, runner=None, cpu_threads=1, dry_run=False,
        stage_logger=lambda *args: log,
    )
    monkeypatch.setattr(_architecture_run, "require_tools", lambda *args: {})

    # runner=None makes any attempted external scan fail the test.
    outputs = s2.run(ctx)
    stage2 = paths.stage_dir(2)
    assert stage2 / "architecture_summary.tsv" in outputs
    assert pd.read_csv(stage2 / "architecture_summary.tsv", sep="\t").empty
    assert pd.read_csv(stage2 / "all_domain_hits.tsv", sep="\t").empty
    assert (stage2 / "pfam.domtblout").read_text() == ""
    assert (stage2 / "custom.domtblout").read_text() == ""


# ---------------------------------------------------------------------------
# Stage 1: domtblout filtering, counters, boundaries
# ---------------------------------------------------------------------------


def test_parse_domtblout_counters_and_boundaries(tmp_path):
    hits, ids, report = (tmp_path / "h.tsv", tmp_path / "ids.txt", tmp_path / "rep.txt")
    df, rep = s1.write_nbd_parse_outputs(
        FIXTURES / "s01_nbd.domtblout", EVAL_NBD, MIN_NBD_ALI, hits, ids, report
    )
    # Counters: 5 non-comment lines (incl. 1 malformed), 2 kept, 1 high i-eval,
    # 1 short alignment; malformed line is counted in total only.
    assert rep == {
        "source": "nbd",
        "total_domtbl_lines": 5,
        "kept_after_filters": 2,
        "filtered_high_i_eval": 1,
        "filtered_short_align": 1,
    }
    # Boundary: i-Evalue exactly at threshold (1e-4) kept; ali_len exactly at
    # min (150) kept -> prot1 (ali 100-249) survives.
    assert set(df["protein_id"]) == {"prot1", "prot3"}
    row = df[df["protein_id"] == "prot1"].iloc[0]
    assert row["dom_i_evalue"] == pytest.approx(1e-4)
    assert row["ali_len"] == 150

    # nbd_hits.tsv has legacy stage-1 column names/order
    assert hits.read_text().splitlines()[0] == (
        "protein_id\tnbd_hmm\tdom_i_evalue\tali_start\tali_end\tali_len"
    )
    assert ids.read_text() == "prot1\nprot3\n"
    # nbd_parse_report.txt: exact legacy key order
    assert report.read_text().splitlines() == [
        "total_domtbl_lines\t5",
        "kept_after_filters\t2",
        "filtered_high_i_eval\t1",
        "filtered_short_align\t1",
        "unique_candidate_proteins\t2",
    ]


def test_parse_domtblout_missing_file():
    df, rep = parse_domtblout("/nonexistent/domtblout", EVAL_NBD, MIN_NBD_ALI, "nbd")
    assert df.empty
    assert rep["total_domtbl_lines"] == 0
    assert rep["kept_after_filters"] == 0


# ---------------------------------------------------------------------------
# Stage 1: eggNOG priority list
# ---------------------------------------------------------------------------


def test_eggnog_priority_keyword_matching(tmp_path):
    out = tmp_path / "eggnog_priority_ids.txt"
    prio = s1.write_eggnog_priority_ids(FIXTURES / "s01_master.tsv", out)
    got = set(prio["protein_id"])
    # NACHT (Description), WD40 (Preferred_name), gasdermin (case-insensitive)
    assert got == {"prot1", "prot3", "prot4"}
    # file has no header, one id per line
    assert out.read_text().splitlines() == list(prio["protein_id"])


def test_eggnog_priority_missing_master_writes_empty(tmp_path):
    out = tmp_path / "prio.txt"
    prio = s1.write_eggnog_priority_ids(tmp_path / "no_master.tsv", out)
    assert prio.empty
    assert out.read_text() == ""


def test_eggnog_priority_master_without_protein_id_writes_empty(tmp_path):
    master = tmp_path / "master.tsv"
    master.write_text("Description\nNACHT protein\n")
    out = tmp_path / "prio.txt"
    prio = s1.write_eggnog_priority_ids(master, out)
    assert prio.empty
    assert out.read_text() == ""


def test_eggnog_priority_master_without_text_columns_writes_empty(tmp_path):
    master = tmp_path / "master.tsv"
    master.write_text("protein_id\tscaffold\np1\tscf1\n")
    out = tmp_path / "prio.txt"
    prio = s1.write_eggnog_priority_ids(master, out)
    assert prio.empty
    assert out.read_text() == ""


# ---------------------------------------------------------------------------
# Stage 1: union / dedup of candidate ids
# ---------------------------------------------------------------------------


def test_write_union_ids_sorted_unique_skips_empty(tmp_path):
    f1 = tmp_path / "nbd.ids"
    f1.write_text("protB\nprotA\nprotB\n\n")
    f2 = tmp_path / "prio.ids"
    f2.write_text("protC\nprotA\n")
    out = tmp_path / "union.txt"
    ids = s1.write_union_ids([f1, f2], out)
    assert ids == ["protA", "protB", "protC"]
    assert out.read_text() == "protA\nprotB\nprotC\n"


def test_write_union_ids_missing_input_tolerated(tmp_path):
    f1 = tmp_path / "a.ids"
    f1.write_text("protA\n")
    out = tmp_path / "union.txt"
    ids = s1.write_union_ids([f1, tmp_path / "missing.ids"], out)
    assert ids == ["protA"]


# ---------------------------------------------------------------------------
# Stage 2: group_domain mapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("NACHT", "NACHT"),
        ("nacht", "NACHT"),
        ("NB-ARC", "NB-ARC"),
        ("Ank_2", "Ank"),
        ("Ankyrin_repeat", "Ank"),
        ("TPR", "TPR"),
        ("TPR-like", "TPR"),
        # legacy \bTPR\b: '_' is a word char, so Pfam-style TPR_1 does NOT match
        ("TPR_1", "TPR_1"),
        ("WD40", "WD40"),
        ("WD_40", "WD40"),
        ("WD-repeat_1", "WD40"),
        ("HEAT", "HEAT"),
        ("Kelch_1", "Kelch"),
        ("LRR_8", "LRR"),
        ("LRRNT_2", "LRR"),
        ("Pkinase", "Pkinase"),  # unmapped domains pass through unchanged
        ("AAA_10", "AAA_10"),
    ],
)
def test_group_domain_mapping(raw, expected):
    assert s2.group_domain(raw) == expected


# ---------------------------------------------------------------------------
# Stage 2: parse + architecture summary
# ---------------------------------------------------------------------------


def _union_prots() -> list[str]:
    return [line[1:].split()[0] for line in (FIXTURES / "s02_union.faa").read_text().splitlines() if line.startswith(">")]


def test_parse_domain_hits_concat_restrict_report(tmp_path):
    all_hits = tmp_path / "all_domain_hits.tsv"
    report = tmp_path / "parse_filter_report.tsv"
    df, rep = s2.parse_domain_hits(
        {"pfam": FIXTURES / "s02_pfam.domtblout", "custom": FIXTURES / "s02_custom.domtblout"},
        EVAL_USE,
        MIN_DOM_ALI,
        set(_union_prots()),
        all_hits,
        report,
    )
    # 8 pfam kept + 1 custom kept = 9; protX is not in the union FASTA -> dropped
    assert df.shape[0] == 8
    assert "protX" not in set(df["protein_id"])
    assert set(df["source"]) == {"pfam", "custom"}

    # two-row report, legacy column order
    assert list(rep.columns) == [
        "source", "total_domtbl_lines", "kept_after_filters",
        "filtered_high_i_eval", "filtered_short_align",
    ]
    pf = rep.set_index("source").loc["pfam"]
    assert (pf["total_domtbl_lines"], pf["kept_after_filters"]) == (10, 8)
    assert (pf["filtered_high_i_eval"], pf["filtered_short_align"]) == (1, 1)
    cu = rep.set_index("source").loc["custom"]
    assert (cu["total_domtbl_lines"], cu["kept_after_filters"]) == (1, 1)


def test_parse_domain_hits_empty_custom_domtblout(tmp_path):
    empty = tmp_path / "custom.domtblout"
    empty.write_text("")
    df, rep = s2.parse_domain_hits(
        {"pfam": FIXTURES / "s02_pfam.domtblout", "custom": empty},
        EVAL_USE, MIN_DOM_ALI, set(_union_prots()),
        tmp_path / "hits.tsv", tmp_path / "rep.tsv",
    )
    cu = rep.set_index("source").loc["custom"]
    assert cu["total_domtbl_lines"] == 0
    assert "custom" not in set(df["source"])


@pytest.fixture(scope="module")
def arch_hits(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("arch")
    df, _ = s2.parse_domain_hits(
        {"pfam": FIXTURES / "s02_pfam.domtblout", "custom": FIXTURES / "s02_custom.domtblout"},
        EVAL_USE, MIN_DOM_ALI, set(_union_prots()),
        tmp / "hits.tsv", tmp / "rep.tsv",
    )
    return df


def test_architecture_summary_columns_and_rows(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, {"protZ"})
    assert list(out.columns) == [
        "protein_id", "domains_raw", "domains_grouped", "has_nbd_pfam",
        "has_repeat", "order_ok", "flag_fp_domains", "flag_rare_lrr",
        "has_nbd_stage1", "has_nbd",
    ]
    # all union proteins retained, in FASTA order
    assert list(out["protein_id"]) == ["protA", "protB", "protC", "protD", "protZ"]


def test_architecture_summary_nbd_repeat_correct_order(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, set()).set_index("protein_id")
    a = out.loc["protA"]
    # (a) NACHT then WD40 -> NBD + repeat in correct order
    assert a["has_nbd_pfam"] and a["has_repeat"] and a["order_ok"]
    assert a["has_nbd"]
    assert a["domains_grouped"] == "NACHT;WD40"  # sorted by start


def test_architecture_summary_reversed_order_fails(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, set()).set_index("protein_id")
    b = out.loc["protB"]
    # (b) WD40 (start 5) before NACHT (start 200) -> order_ok False
    assert b["has_nbd_pfam"] and b["has_repeat"] and not b["order_ok"]


def test_architecture_summary_zero_hit_protein_retained(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, {"protZ"}).set_index("protein_id")
    z = out.loc["protZ"]
    # (c) zero-domain protein kept; has_nbd comes from stage-1 ids
    assert z["domains_raw"] == "" and z["domains_grouped"] == ""
    assert not z["has_nbd_pfam"] and not z["has_repeat"]
    assert z["order_ok"] and not z["flag_fp_domains"] and not z["flag_rare_lrr"]
    assert z["has_nbd_stage1"] and z["has_nbd"]


def test_architecture_summary_fp_pattern_flagging(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, set()).set_index("protein_id")
    # (d) ABC_tran and Pkinase are FP domains
    assert out.loc["protC", "flag_fp_domains"]
    assert out.loc["protD", "flag_fp_domains"]
    # LRR-only-ish protein flagged rare LRR; no FP on clean NLR-like protA
    assert out.loc["protC", "flag_rare_lrr"]
    assert not out.loc["protA", "flag_fp_domains"]
    # protD: Pkinase is neither NBD nor repeat; unmapped custom domain passes through
    assert not out.loc["protD", "has_nbd"]
    assert not out.loc["protD", "has_repeat"]
    assert out.loc["protD", "domains_grouped"] == "Pkinase;CustomFungal1"


def test_architecture_summary_has_nbd_union_of_stage1_and_pfam(arch_hits):
    out = s2.summarize_architecture(_union_prots(), arch_hits, {"protC"}).set_index("protein_id")
    # stage1-only NBD call
    assert out.loc["protC", "has_nbd"]
    # pfam-only NBD call
    assert out.loc["protA", "has_nbd"] and not out.loc["protA", "has_nbd_stage1"]


def test_architecture_summary_empty_hits():
    cols = ["protein_id", "domain", "acc", "dom_i_evalue", "start", "end", "ali_len", "source"]
    empty = pd.DataFrame(columns=cols)
    out = s2.summarize_architecture(["p1", "p2"], empty, {"p2"})
    assert out.shape[0] == 2
    assert out["has_nbd"].tolist() == [False, True]
    assert out["order_ok"].tolist() == [True, True]


def test_write_sorted_ids(tmp_path):
    out = tmp_path / "stage1_nbd.ids"
    s2.write_sorted_ids({"protB", "protA", "protA"}, out)
    assert out.read_text() == "protA\nprotB\n"


def test_read_id_set_missing_or_empty(tmp_path):
    missing = tmp_path / "nope.ids"
    assert s2.read_id_set(missing) == set()
    empty = tmp_path / "empty.ids"
    empty.write_text("")
    assert s2.read_id_set(empty) == set()
