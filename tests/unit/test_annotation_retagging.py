"""Retagging replaces FuNLR evidence without retaining obsolete calls."""
from pathlib import Path

import pytest

from funlr.stages.annotation_export import NLR_FIELDS, tagged_gff


def mrna(attributes):
    return f"chr1\tannotation\tmRNA\t12\t345\t.\t-\t.\t{attributes}\n"


@pytest.mark.parametrize("parent", ["Parent=gone;", ""])
def test_stale_calls_removed_without_damaging_other_attributes(tmp_path, parent):
    retained = (
        f"ID=tx;{parent}user_NLR_tier=keep;NLR_tier_extra=keep;"
        "Note=prior%3BNLR_tier%3Dold;product=immune%20candidate"
    )
    old_tags = ";".join(f"{key}=STALE" for key in NLR_FIELDS)
    source = tmp_path / "source.gff3"
    source.write_text("##gff-version 3\n" + mrna(retained + ";" + old_tags))
    target = tmp_path / "retagged.gff3"
    tagged_gff(source, target, {})
    assert target.read_text() == "##gff-version 3\n" + mrna(retained)
    assert source.read_text().endswith(mrna(retained + ";" + old_tags))


def test_retagging_replaces_calls_and_is_byte_idempotent(tmp_path):
    source = tmp_path / "source.gff3"
    source.write_text(mrna("ID=tx;Parent=gene%3A1;NLR_tier=OLD;NLR_flags=OLD"))
    calls = {"gene:1": {"NLR_tier": "TIER_1B_NEEDS_REVIEW", "NLR_flags": "A;B"}}
    first, second = tmp_path / "first.gff3", tmp_path / "second.gff3"
    tagged_gff(source, first, calls)
    tagged_gff(first, second, calls)
    assert first.read_bytes() == second.read_bytes()
    assert first.read_text() == mrna(
        "ID=tx;Parent=gene%3A1;NLR_tier=TIER_1B_NEEDS_REVIEW;NLR_flags=A%3BB"
    )


def test_untagged_unmatched_and_non_mrna_records_remain_unchanged(tmp_path):
    text = (
        "##gff-version 3\n# NLR_tier=comment\n"
        + mrna("ID=unmatched;Parent=other;;Note=untouched;")
        + "chr1\tannotation\tgene\t12\t345\t.\t-\t.\tID=other;NLR_tier=external\n"
        + "chr1\tannotation\tCDS\t12\t345\t.\t-\t0\tParent=unmatched;NLR_flags=external\n"
    )
    source, target = tmp_path / "source.gff3", tmp_path / "target.gff3"
    source.write_text(text)
    tagged_gff(source, target, {})
    assert target.read_text() == text
