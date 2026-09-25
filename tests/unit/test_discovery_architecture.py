"""Updated reference contracts and boundary cases that the real reference omits."""
from pathlib import Path

import pandas as pd
import pytest

from funlr.core.errors import StageError
from funlr.stages._inputs import reconcile_coordinates, parse_header_eggnog
from funlr.stages._discovery import parse_calls, merge_calls, enrich_pfams, description_priority
from funlr.stages._architecture import build_architecture
from funlr.stages._discovery_run import build_union, write_ids
from funlr.stages._settings import DEFAULTS


def domtbl(domain, protein, start, end, *, evalue=1e-8, bits=50, hmm_len=300, accession="-"):
    fields = ["0"] * 22
    for pos, value in {0: domain, 1: accession, 2: hmm_len, 3: protein, 4: "-", 5: 1000,
                       12: evalue, 13: bits, 17: start, 18: end}.items():
        fields[pos] = str(value)
    return " ".join(fields) + "\n"


def stages(tmp_path, proteins=("p1",)):
    s0, s1, s2 = (tmp_path / f"stage{n}" for n in range(3))
    for folder in (s0, s1, s2):
        folder.mkdir()
    fasta = "".join(f">{protein}\n{'M' * 600}\n" for protein in proteins)
    (s0 / "proteins_clean.faa").write_text(fasta)
    (s1 / "union_candidates.faa").write_text(fasta)
    for name in ("pfam.domtblout", "custom.domtblout", "asm.domtblout", "pfam_bg.domtblout", "stage1_nbd.ids"):
        (s2 / name).write_text("")
    return s0, s1, s2


def test_reconcile_protein_ids_prefers_cds_then_transcript_then_gene(tmp_path):
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(">cds1\nM\n>tx1\nM\n>tx2\nM\n>gene3\nM\n")
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "ctg\tx\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "ctg\tx\tCDS\t1\t90\t.\t+\t0\tParent=tx1;protein_id=prefix|cds1\n"
        "ctg\tx\tmRNA\t100\t190\t.\t+\t.\tID=tx2;Parent=gene2\n"
        "ctg\tx\tCDS\t100\t190\t.\t+\t0\tParent=tx2;orig_protein_id=absent\n"
        "ctg\tx\tmRNA\t200\t290\t.\t-\t.\tID=tx3;Parent=gene3\n"
        "ctg\tx\tmRNA\t300\t390\t.\t+\t.\tID=tx4;Parent=gene4\n"
        "ctg\tx\tCDS\t300\t390\t.\t+\t0\tParent=tx4;protein_id=missing4\n")
    coords, counts = reconcile_coordinates(gff, proteins)
    assert coords["protein_id"].tolist() == ["cds1", "tx2", "gene3", "missing4"]
    assert coords["mapping_source"].tolist() == ["CDS_PROTEIN_ID", "TRANSCRIPT_ID_MATCH", "GENE_ID_MATCH", "CDS_PROTEIN_ID_NOT_IN_FASTA"]
    assert counts["rows_with_protein_id_not_in_fasta"] == 1


def test_eggnog_resolves_named_columns_and_requires_header(tmp_path):
    path = tmp_path / "annotations.tsv"
    path.write_text("## provenance\n#query\tPFAMs\tDescription\np1\tNACHT,WD40\treceptor\n")
    row = parse_header_eggnog(path).iloc[0]
    assert row["PFAMs"] == "NACHT,WD40"
    assert row["Description"] == "receptor"
    assert row["GOs"] == ""
    path.write_text("p1\tNACHT\treceptor\n")
    with pytest.raises(StageError, match="#query"):
        parse_header_eggnog(path)


def test_strict_relaxed_boundaries_and_pfam_discovery(tmp_path):
    s0, s1, _ = stages(tmp_path)
    (s1 / "nbd_whole.domtblout").write_text(
        domtbl("NACHT", "strict", 10, 189, evalue=1e-5) +
        domtbl("NACHT", "relaxed", 10, 129, evalue=1e-3) +
        domtbl("NACHT", "too_short", 10, 128, evalue=1e-8) +
        domtbl("NACHT", "too_weak", 10, 300, evalue=0.0011))
    (s1 / "pfam_nbd.domtblout").write_text(domtbl("NB-ARC", "pfam", 1, 120, evalue=1e-3))
    parse_calls(s0, s1, DEFAULTS)
    assert (s1 / "nbd_candidate_ids_strict.txt").read_text() == "strict\n"
    assert (s1 / "nbd_candidate_ids_relaxed.txt").read_text() == "strict\nrelaxed\n"
    assert (s1 / "pfam_nbd_candidate_ids.txt").read_text() == "pfam\n"
    assert len(pd.read_csv(s1 / "nbd_hits_detailed.tsv", sep="\t")) == 4


def test_no_nbd_hits_keep_headed_outputs_and_enrichment(tmp_path):
    s0, s1, _ = stages(tmp_path)
    (s0 / "master_table.tsv").write_text("protein_id\tPFAMs\np1\tWD40\n")
    (s1 / "nbd_whole.domtblout").write_text("# zero hits\n")
    (s1 / "pfam_nbd.domtblout").write_text("")
    parse_calls(s0, s1, DEFAULTS)
    merge_calls(s0, s1, DEFAULTS)
    enrich_pfams(s0, s1, DEFAULTS)
    for name in ("nbd_hits.tsv", "nbd_hits_detailed.tsv", "pfam_nbd_hits.tsv", "nbd_hits_with_master_strict.tsv", "pfam_enriched_in_nbd.tsv"):
        assert pd.read_csv(s1 / name, sep="\t").empty
    assert (s1 / "nbd_candidate_ids_strict.txt").read_text() == ""
    assert (s1 / "pfam_priority_ids.txt").read_text() == ""


def test_optional_eggnog_empty_light_table_keeps_empty_description_ids(tmp_path):
    s0, s1, _ = stages(tmp_path)
    (s0 / "eggnog_light.tsv").write_text("protein_id\tCOG_category\tDescription\tPreferred_name\tGOs\tKEGG_ko\tPFAMs\n")
    description_priority(s0, s1, DEFAULTS)
    assert (s1 / "desc_priority_ids.txt").read_text() == ""


def test_balanced_union_and_all_structural_sources_respect_empty(tmp_path):
    s0, s1, _ = stages(tmp_path, ("nbd", "enriched", "description", "pfam"))
    for name, ids in {"nbd_candidate_ids_relaxed.txt": ["nbd", "absent"],
                      "nbd_candidate_ids_strict.txt": ["nbd"], "pfam_nbd_candidate_ids.txt": ["pfam"],
                      "pfam_priority_ids.txt": ["enriched"], "desc_priority_ids.txt": ["description", "nbd"]}.items():
        write_ids(s1 / name, ids)
    assert build_union(s0, s1, DEFAULTS) == ["enriched", "nbd", "pfam"]
    assert (s1 / "union_candidate_ids.missing_from_fasta.txt").read_text() == "absent\n"
    values = {**DEFAULTS, "INCLUDE_DESC_KEYWORDS_IN_UNION": 1, "DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT": 0}
    assert "description" in build_union(s0, s1, values)
    (s1 / "pfam_nbd_candidate_ids.txt").write_text("")
    values.update(DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT=1, DESC_KEYWORD_STRUCTURAL_HINT_MODE="ALL",
                  DESC_KEYWORD_STRUCTURAL_HINT_SOURCES="RELAXED_NBD,PFAM_NBD")
    build_union(s0, s1, values)
    assert (s1 / "desc_priority_ids.structural_hint.txt").read_text() == ""


def test_architecture_effector_state_does_not_leak_to_empty_candidates(tmp_path):
    s0, s1, s2 = stages(tmp_path, ("empty_first", "effector", "empty_after", "orientation_bad"))
    (s2 / "pfam.domtblout").write_text(
        domtbl("NACHT", "effector", 100, 279) + domtbl("TPR_1", "effector", 289, 320) +
        domtbl("NACHT", "orientation_bad", 100, 279) + domtbl("TPR_1", "orientation_bad", 288, 320))
    (s2 / "custom.domtblout").write_text(domtbl("Goodbye__class__test", "effector", 10, 90))
    values = {**DEFAULTS, "EVAL_USE": "1e-3"}
    out = build_architecture(s0, s1, s2, values).set_index("protein_id")
    assert bool(out.loc["effector", "has_effector_nterm"])
    assert bool(out.loc["effector", "has_ssfr_sensor"])
    assert not bool(out.loc["effector", "order_invalid"])
    assert bool(out.loc["orientation_bad", "order_invalid"])
    for pid in ("empty_first", "empty_after"):
        assert not bool(out.loc[pid, "has_effector"])
        assert not bool(out.loc[pid, "has_effector_nterm"])
        assert not bool(out.loc[pid, "has_nbd"])
    assert out.loc["effector", "nbd_confidence"] == "HIGH"


def test_asm_thresholds_and_pfam_medium_confidence(tmp_path):
    s0, s1, s2 = stages(tmp_path, ("accepted", "low_score", "short"))
    (s2 / "pfam.domtblout").write_text(domtbl("NACHT", "accepted", 50, 149, evalue=1e-4))
    (s2 / "asm.domtblout").write_text(
        domtbl("ASM_TEST", "accepted", 1, 12, evalue=1e-3, bits=20) +
        domtbl("ASM_TEST", "low_score", 1, 12, evalue=1e-3, bits=19.9) +
        domtbl("ASM_TEST", "short", 1, 11, evalue=1e-3, bits=20))
    out = build_architecture(s0, s1, s2, {**DEFAULTS, "EVAL_USE": "1e-3"}).set_index("protein_id")
    assert out.loc["accepted", "nbd_confidence"] == "MEDIUM"
    assert out["asm_present"].to_dict() == {"accepted": True, "low_score": False, "short": False}
    assert not out["has_sensor"].any()  # ASM is supporting evidence only.


def test_empty_architecture_retains_full_schema_and_asm_summary(tmp_path):
    s0, s1, s2 = stages(tmp_path, ())
    out = build_architecture(s0, s1, s2, {**DEFAULTS, "EVAL_USE": "1e-3"})
    assert out.empty and len(out.columns) == 30
    assert pd.read_csv(s2 / "asm_hits_summary.tsv", sep="\t").empty
    assert pd.read_csv(s2 / "architecture_summary.tsv", sep="\t").empty
    assert pd.read_csv(s2 / "pfam_enriched_in_nbd_stage2.tsv", sep="\t").empty
