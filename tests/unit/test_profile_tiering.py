"""Updated reference decision rules, confidence boundaries, and export contracts."""
import os
from pathlib import Path

import pandas as pd
import pytest

from funlr.config import FunlrConfig
from funlr.core.errors import StageError
from funlr.stages.tiering import TierPolicy, build_profile_table, write_profile_outputs
from funlr.stages import stage4_export


def row(**overrides):
    result = dict(protein_id="p", has_nbd=True, has_nbd_any=True, nbd_ok=True,
                  nbd_confidence="HIGH", has_sensor=False, has_ssfr_sensor=False,
                  has_non_ssfr_sensor=False, flags="PASS", protein_length=900)
    result.update(overrides)
    return result


@pytest.mark.parametrize("changes,tier", [
    ({"has_nbd_any": False, "has_sensor": True, "has_lrr": True}, "4A_REPEAT_ONLY_NO_NBD"),
    ({"has_nbd_any": False}, "4C_NO_NBD_NO_REPEAT_LOW_SIGNAL"),
    ({"has_lrr": True}, "4D_EXCLUDED_LRR"),
    ({"order_invalid": True}, "4E_ORIENTATION_INVALID"),
    ({"non_nlr_annot": True, "nbd_confidence": "MEDIUM"}, "4B_LIKELY_STAND_HOUSEKEEPING"),
    ({"flags": "VERY_SHORT_PROTEIN"}, "4F_VERY_SHORT_FRAGMENT"),
    ({"has_sensor": True, "has_ssfr_sensor": True}, "1A_HIGH_CONFIDENCE"),
    ({"has_sensor": True, "has_ssfr_sensor": True, "flags": "HARD_END"}, "1B_NEEDS_REVIEW"),
    ({"has_sensor": True, "has_non_ssfr_sensor": True}, "1C_CANONICAL_NON_SSFR_SENSOR"),
    ({"has_any_repeat_region": True}, "2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE"),
    ({"has_any_repeat_region": True, "flags": "HARD_END"}, "2A_HIGH_PRIORITY_RESCUE"),
    ({"flags": "HARD_END;RESCUE_NO_FLANK"}, "2B_RESCUE_CANDIDATE"),
    ({"nbd_ok": False, "nbd_confidence": "LOW", "flags": "HARD_END;RESCUE_NO_FLANK"}, "2C_LOW_PRIORITY_FRAGMENT"),
    ({"nbd_ok": False, "nbd_confidence": "LOW", "flags": "HARD_END"}, "2B_RESCUE_CANDIDATE"),
    ({"nbd_ok": False, "nbd_confidence": "LOW", "has_sensor": True, "has_ssfr_sensor": True}, "3B_ARCHITECTURAL_VARIANT"),
    ({"flags": "NBD_ONLY_UNKNOWN_SENSOR"}, "2B_RESCUE_CANDIDATE"),
    ({"flags": "SHORT_PROTEIN"}, "2C_LOW_PRIORITY_FRAGMENT"),
    ({"integrated_domain_cterm": True}, "3A_INTEGRATED_DECOY"),
    ({}, "3B_ARCHITECTURAL_VARIANT"),
])
def test_updated_tier_decision_order(changes, tier):
    assert TierPolicy().assign_tier(row(**changes)) == "TIER_" + tier


def test_hifi_context_and_unknown_sensor():
    policy = TierPolicy(profile="HIFI", hard_end_bp=1000, subtelo_bp=50000)
    assert policy.assign_tier(row(flags="NBD_ONLY_UNKNOWN_SENSOR")) == "TIER_4G_SUSPECTED_PSEUDOGENE"
    assert policy.assign_tier(row(flags="HARD_END")) == "TIER_2B_RESCUE_CANDIDATE"
    assert policy.assign_tier(row(flags="HARD_END", has_any_repeat_region=True)) == "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE"
    # Subtelomeric location alone is biological context, not a tier demotion.
    assert policy.assign_tier(row(flags="SUBTELO", has_sensor=True, has_ssfr_sensor=True)) == "TIER_1A_HIGH_CONFIDENCE"


def test_repeat_aware_orientation_exception_keeps_zero_score():
    policy = TierPolicy()
    candidate = row(order_invalid=True, has_any_repeat_region=True, flags="ORIENTATION_INVALID")
    candidate["tier"] = policy.assign_tier(candidate)
    assert candidate["tier"] == "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE"
    assert policy.rescue_score(candidate) == 0


def test_flags_use_end_distance_length_boundaries_and_unknown_orientation():
    policy = TierPolicy(hard_end_bp=1000, subtelo_bp=5000)
    candidate = row(contig_len=100000, start=1001, end=4000, domains_grouped="NACHT", protein_length=332)
    assert policy.build_flags(candidate) == ["NBD_ONLY_UNKNOWN_SENSOR", "HARD_END", "RESCUE_POSSIBLE", "SUBTELO", "VERY_SHORT_PROTEIN", "NBD_ONLY_TRUNCATED"]
    assert policy.length_flags(row(protein_length=333)) == ["SHORT_PROTEIN"]
    assert policy.length_flags(row(protein_length=500)) == []
    candidate.update(start=1002, protein_length=900, has_sensor=True, sensor_is_cterm_of_nbd=pd.NA)
    assert policy.build_flags(candidate) == ["SUBTELO"]
    candidate["sensor_is_cterm_of_nbd"] = False
    assert policy.build_flags(candidate) == ["DOMAIN_ORDER_WEIRD", "SUBTELO"]


@pytest.mark.parametrize("flag,increment", [("HAS_EFFECTOR_NTERM",20), ("HAS_EFFECTOR_ANY",0),
    ("HARD_END",15), ("RESCUE_POSSIBLE",10), ("RESCUE_NO_FLANK",-30),
    ("REPEAT_NEARBY",10), ("SUBTELO",5), ("ASM_PRESENT",25),
    ("INTEGRATED_DOMAIN_CTERM",-10), ("SHORT_PROTEIN",10),
    ("VERY_SHORT_PROTEIN",-30), ("VERY_LONG_PROTEIN",-10),
    ("FP_DOMAIN_PRESENT",-25), ("FP_DOMAIN_OVERLAPS_NBD",-50), ("NON_NLR_ANNOT",-15)])
def test_score_terms_preserve_substring_semantics(flag, increment):
    candidate = row(tier="TIER_1C_CANONICAL_NON_SSFR_SENSOR", flags=flag)
    assert TierPolicy().rescue_score(candidate) == max(0, min(100, 60 + increment))


def inputs(arch_rows=None):
    records = arch_rows or [dict(protein_id="p", has_nbd=True, has_sensor=True,
        has_ssfr_sensor=True, has_non_ssfr_sensor=False, domains_raw="NACHT;WD40",
        domains_grouped="NACHT;WD40", nbd_confidence="LOW", has_nbd_stage1=True)]
    arch = pd.DataFrame(records)
    master = pd.DataFrame([dict(protein_id=pid, scaffold="s", start=40000+i*15000,
        end=42000+i*15000, strand="+", protein_length=900) for i,pid in enumerate(arch.protein_id)])
    contigs = pd.DataFrame([{"scaffold": "s", "contig_len": 200000}])
    return master, arch, contigs


def test_explicit_all_low_batch_is_never_promoted():
    args = inputs()
    frame = build_profile_table(*args, TierPolicy())
    assert frame.nbd_confidence.tolist() == ["LOW"]
    assert frame.nbd_ok.tolist() == [False]
    assert frame.tier.tolist() == ["TIER_3B_ARCHITECTURAL_VARIANT"]
    assert args[1].columns.tolist()[-2:] == ["nbd_confidence", "has_nbd_stage1"]


def test_absent_confidence_uses_legacy_evidence_only():
    master, arch, contigs = inputs()
    frame = build_profile_table(master, arch.drop(columns="nbd_confidence"), contigs, TierPolicy())
    assert frame.nbd_confidence.tolist() == ["HIGH"]
    assert frame.tier.tolist() == ["TIER_1A_HIGH_CONFIDENCE"]


def test_low_confidence_does_not_acquire_repeat_rescue_context():
    master, arch, contigs = inputs()
    low = arch.iloc[0].to_dict()
    low.update(protein_id="q", has_sensor=False, has_ssfr_sensor=False, domains_raw="NACHT", domains_grouped="NACHT")
    arch = pd.concat([arch, pd.DataFrame([low])], ignore_index=True)
    master = pd.concat([master, pd.DataFrame([dict(protein_id="q",scaffold="s",start=45000,end=47000,strand="+",protein_length=900)])],ignore_index=True)
    frame = build_profile_table(master,arch,contigs,TierPolicy()).set_index("protein_id")
    assert not frame.loc["q", "repeat_nearby"]
    arch.loc[arch.protein_id.eq("q"), "nbd_confidence"] = "MEDIUM"
    frame = build_profile_table(master,arch,contigs,TierPolicy()).set_index("protein_id")
    assert frame.loc["q", "repeat_nearby"]
    assert frame.loc["q", "tier"] == "TIER_2B_RESCUE_CANDIDATE"


def test_empty_frame_emits_complete_schema_and_zero_summaries(tmp_path):
    master, arch, contigs = inputs()
    frame = build_profile_table(master, arch.head(0), contigs, TierPolicy())
    outputs = write_profile_outputs(frame, tmp_path, TierPolicy())
    assert outputs and all(p.is_file() for p in outputs)
    assert pd.read_csv(tmp_path/"tables/tiered_candidates.tsv",sep="\t").empty
    assert pd.read_csv(tmp_path/"tables/rescue_priority.tsv",sep="\t").empty
    summary = pd.read_csv(tmp_path/"summaries/universe_summary.tsv",sep="\t")
    assert summary["count"].sum() == 0 and summary["percent"].sum() == 0
    assert (tmp_path/"beds/nlr_candidates.bed").read_text() == ""


def test_rescue_policy_and_ties_are_explicit(tmp_path):
    master,arch,contigs=inputs()
    frame=build_profile_table(master,arch,contigs,TierPolicy())
    frame=pd.concat([frame.assign(protein_id="z",tier="TIER_1A_HIGH_CONFIDENCE",rescue_priority=30),
                     frame.assign(protein_id="b",tier="TIER_3B_ARCHITECTURAL_VARIANT",rescue_priority=30),
                     frame.assign(protein_id="a",tier="TIER_2B_RESCUE_CANDIDATE",rescue_priority=30)],ignore_index=True)
    write_profile_outputs(frame,tmp_path/"default",TierPolicy())
    assert pd.read_csv(tmp_path/"default/tables/rescue_priority.tsv",sep="\t").protein_id.tolist()==["a"]
    write_profile_outputs(frame,tmp_path/"all",TierPolicy(skip_rescue_for_tier1a=False,include_tier3_in_rescue=True))
    assert pd.read_csv(tmp_path/"all/tables/rescue_priority.tsv",sep="\t").protein_id.tolist()==["a","b","z"]


@pytest.mark.parametrize("columns,lines,expected", [
    ("protein_id\tnbd_ok\tnbd_confidence", ["p\tFalse\tHIGH","q\tTrue\tLOW"], ["q"]),
    ("protein_id\tnbd_confidence", ["p\tLOW","q\tmedium"], ["q"]),
    ("protein_id\ttier", ["p\tTIER_2A","q\tTIER_2B"], ["p","q"]),
])
def test_rescue_export_filter_precedence(tmp_path,columns,lines,expected):
    table, ids = tmp_path/"rescue.tsv",tmp_path/"rescue.ids"
    table.write_text(columns+"\n"+"\n".join(lines)+"\n")
    stage4_export.extract_rescue_ids(table,ids,True)
    assert ids.read_text().splitlines()==expected


def test_export_rejects_successful_but_incomplete_tool_output(tmp_path):
    class Runner:
        def resolve(self,key): return key
        def run(self,argv,**kwargs):
            assert kwargs["check"] is True and "-r" not in argv
            Path(kwargs["stdout_path"]).write_text(">p\nAAAA\n")
    class Context:
        runner=Runner()
    ids=tmp_path/"in.ids"
    ids.write_text("p\nq\n")
    with pytest.raises(StageError,match="missing=\\['q'\\]"):
        stage4_export.export_faa(Context(),ids,tmp_path/"prot.faa",tmp_path/"out.faa")


@pytest.mark.skipif(not os.environ.get("FUNLR_REFERENCE_DIR"), reason="requires separate reference data")
def test_reference_stage3_all_columns_and_rescue_membership(tmp_path):
    reference=Path(os.environ["FUNLR_REFERENCE_DIR"])
    master=pd.read_csv(reference/"stage0/master_table.tsv",sep="\t")
    arch=pd.read_csv(reference/"stage2/architecture_summary.tsv",sep="\t")
    contigs=pd.read_csv(reference/"stage0/contig_lengths.tsv",sep="\t",header=None,names=["scaffold","contig_len"])
    config=FunlrConfig()
    frame=build_profile_table(master,arch,contigs,TierPolicy.from_settings(config.scientific_settings()),settings=config.scientific_settings())
    write_profile_outputs(frame,tmp_path,TierPolicy())
    for name in ("tiered_candidates.tsv","rescue_priority.tsv"):
        actual=pd.read_csv(tmp_path/"tables"/name,sep="\t").set_index("protein_id").sort_index().sort_index(axis=1)
        expected=pd.read_csv(reference/"stage3/tables"/name,sep="\t").set_index("protein_id").sort_index().sort_index(axis=1)
        pd.testing.assert_frame_equal(actual,expected,check_dtype=False)
    for name in ("tier_summary.tsv","flag_counts.tsv","universe_summary.tsv"):
        assert (tmp_path/"summaries"/name).read_bytes()==(reference/"stage3/summaries"/name).read_bytes()
    assert (tmp_path/"beds/nlr_candidates.bed").read_bytes()==(reference/"stage3/beds/nlr_candidates.bed").read_bytes()
