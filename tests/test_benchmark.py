"""Reference coverage cannot silently become whole-proteome precision."""
import hashlib
import json
from pathlib import Path

import pytest

from funlr.benchmark import compare_reference


@pytest.fixture
def panel(tmp_path):
    source = Path(__file__).parents[1] / "examples/benchmark"
    for name in ("input.faa", "called.tsv", "reference.json"):
        (tmp_path / name).write_bytes((source / name).read_bytes())
    return tmp_path


def compare(panel, **kwargs):
    return compare_reference(panel / "reference.json", panel / "called.tsv", panel / "input.faa", **kwargs)


def edit(panel, function):
    path = panel / "reference.json"
    data = json.loads(path.read_text())
    function(data)
    path.write_text(json.dumps(data))


def test_partial_synthetic_panel_has_recovery_but_no_precision_and_is_read_only(panel):
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in panel.iterdir()}
    result = compare(panel)
    assert result["metrics"] == {"positive_recovery": 0.5, "precision": None, "specificity": None}
    assert result["counts"] == {"input_units": 4, "called_units": 2, "positive_present": 2,
                                "positive_recovered": 1, "negative_present": 1, "negative_called": 0,
                                "unavailable_references": 1, "unresolved_or_unlisted_calls": 1}
    assert result["unknown_called_ids"] == ["unresolved"]
    assert result["independence"] == "SYNTHETIC"
    assert result["strata"]["training_overlap"]["yes"]["positive_recovered"] == 1
    assert result["strata"]["training_overlap"]["no"]["positive_recovery"] == 0
    assert result["provenance"]["input"]["sha256"] == hashlib.sha256((panel / "input.faa").read_bytes()).hexdigest()
    assert {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in panel.iterdir()} == before


def test_unlisted_calls_are_unknown_even_with_a_negative_reference(panel):
    edit(panel, lambda d: d.update(records=[r for r in d["records"] if r["protein_id"] != "unresolved"]))
    result = compare(panel)
    assert result["unknown_called_ids"] == ["unresolved"]
    assert result["counts"]["negative_called"] == 0
    assert result["metrics"]["precision"] is None


def test_positive_only_zero_recovery_does_not_crash_or_invent_negatives(panel):
    edit(panel, lambda d: d.update(records=[r for r in d["records"] if r["expected_call"] == "positive"]))
    (panel / "called.tsv").write_text("protein_id\nnear_neighbor\n")
    result = compare(panel)
    assert result["metrics"]["positive_recovery"] == 0
    assert result["metrics"]["precision"] is None
    assert result["counts"]["negative_called"] == 0
    assert result["unknown_called_ids"] == ["near_neighbor"]


def test_empty_called_set_has_valid_zero_recovery(panel):
    (panel / "called.tsv").write_text("protein_id\n")
    assert compare(panel)["metrics"] == {"positive_recovery": 0, "precision": None, "specificity": None}


def test_exhaustive_panel_precision_is_explicit_and_has_known_denominator(panel):
    def make_complete(d):
        d["exhaustive_labels"] = True
        next(r for r in d["records"] if r["protein_id"] == "unresolved")["expected_call"] = "negative"
    edit(panel, make_complete)
    result = compare(panel)
    assert result["metrics"] == {"positive_recovery": 0.5, "precision": 0.5, "specificity": 0.5}
    assert result["precision_scope"] == "Entire declared evaluation FASTA"
    assert result["counts"]["negative_called"] == 1


def test_false_exhaustive_claim_is_rejected(panel):
    edit(panel, lambda d: d.update(exhaustive_labels=True))
    with pytest.raises(ValueError, match="every evaluation FASTA ID"):
        compare(panel)


def test_locus_mode_requires_explicit_same_unit_and_never_falls_back_to_genes(panel):
    def loci(d):
        d["unit"] = "locus"
        for row in d["records"]:
            row["locus_id"] = row.pop("protein_id")
    edit(panel, loci)
    (panel / "called.tsv").write_text("locus_id\npositive_a\n")
    with pytest.raises(ValueError, match="--id-column locus_id"):
        compare(panel)
    result = compare(panel, id_column="locus_id")
    assert result["metrics"]["positive_recovery"] == 0.5
    assert "curator assertions" in result["presence_check"]


@pytest.mark.parametrize("data,match", [
    ("protein_id\ttier\npositive_a\tTIER1\n", "exactly one column"),
    ("gene_id\npositive_a\n", "exactly one column"),
    ("protein_id\npositive_a\npositive_a\n", "Duplicate called"),
    ("protein_id\nmade_up_fusion\n", "absent from evaluation FASTA"),
    ("protein_id\n\"positive_a\nextra\"\n", "single-line"),
    ("", "empty"),
])
def test_ambiguous_or_inconsistent_prediction_sets_are_rejected(panel, data, match):
    (panel / "called.tsv").write_text(data)
    with pytest.raises(ValueError, match=match):
        compare(panel)


@pytest.mark.parametrize("mutation,match", [
    (lambda d: d.update(input_sha256="0" * 64), "input_sha256"),
    (lambda d: d["records"][0].update(reference_present=False), "reference_present disagrees"),
    (lambda d: d["records"][0].update(evidence_source=None), "requires evidence_source"),
    (lambda d: d["records"].append(d["records"][0]), "Duplicate reference"),
    (lambda d: d["records"][0].update(training_overlap="maybe"), "training_overlap"),
    (lambda d: d.update(synthetic=False), "Synthetic reference"),
    (lambda d: d.update(exhaustive_labels="false"), "JSON boolean"),
    (lambda d: d.update(schema_version=True), "schema_version"),
    (lambda d: d.update(extra="typo"), "fields must be exactly"),
    (lambda d: d.update(unit=[]), "Reference unit"),
    (lambda d: d["records"][0].update(expected_call=[]), "expected_call"),
    (lambda d: d["records"][0].update(evidence_category={}), "evidence_category"),
    (lambda d: d["records"][0].update(training_overlap=[]), "training_overlap"),
])
def test_reference_claims_are_validated_before_metrics(panel, mutation, match):
    edit(panel, mutation)
    with pytest.raises(ValueError, match=match):
        compare(panel)


def test_changed_input_at_same_path_cannot_reuse_reference_claim(panel):
    path = panel / "input.faa"
    path.write_bytes(path.read_bytes().replace(b"MALWTRK", b"AAAAAAA"))
    with pytest.raises(ValueError, match="input_sha256"):
        compare(panel)


def test_duplicate_json_keys_do_not_silently_change_a_label(panel):
    path = panel / "reference.json"
    path.write_text(path.read_text().replace('"expected_call": "positive"',
                                           '"expected_call": "negative", "expected_call": "positive"', 1))
    with pytest.raises(ValueError, match="Duplicate reference JSON field"):
        compare(panel)


@pytest.mark.parametrize("fasta,match", [
    (b">same\nMA\n>same\nMT\n", "Duplicate input"),
    (b">empty\n>next\nMA\n", "Empty input FASTA sequence"),
    (b"MA\n>next\nMT\n", "before first"),
])
def test_invalid_fasta_universe_cannot_produce_scores(panel, fasta, match):
    (panel / "input.faa").write_bytes(fasta)
    edit(panel, lambda d: d.update(input_sha256=hashlib.sha256(fasta).hexdigest()))
    with pytest.raises(ValueError, match=match):
        compare(panel)


def test_export_helper_excludes_fusions_and_refuses_overwrite(tmp_path):
    import runpy
    helper = Path(__file__).parents[1] / "examples/benchmark/export_strict_proteins.py"
    export = runpy.run_path(str(helper))["export_strict_proteins"]
    source = tmp_path / "nlr_strict_candidates.tsv"
    source.write_text("protein_id\tis_fusion_model\noriginal\tFalse\nfusion_1\tTrue\n")
    target = tmp_path / "called.tsv"
    assert export(source, target) == 1
    assert target.read_bytes() == b"protein_id\noriginal\n"
    with pytest.raises(FileExistsError):
        export(source, target)
    assert target.read_bytes() == b"protein_id\noriginal\n"


def test_export_helper_refuses_all_candidate_report_and_invalid_flags_before_write(tmp_path):
    import runpy
    helper = Path(__file__).parents[1] / "examples/benchmark/export_strict_proteins.py"
    export = runpy.run_path(str(helper))["export_strict_proteins"]
    source = tmp_path / "nlr_final_report.tsv"
    source.write_text("protein_id\tis_fusion_model\noriginal\tFalse\n")
    target = tmp_path / "called.tsv"
    with pytest.raises(ValueError, match="all-candidate"):
        export(source, target)
    assert not target.exists()
    strict = tmp_path / "nlr_strict_candidates.tsv"
    strict.write_text("protein_id\tis_fusion_model\noriginal\tunknown\n")
    with pytest.raises(ValueError, match="invalid is_fusion_model"):
        export(strict, target)
    assert not target.exists()
