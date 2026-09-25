"""Input measurements preserve bytes and describe evidence without grading it."""
import hashlib
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from funlr.input_qc import collect_input_qc, write_input_qc


def _inputs(tmp_path, genome=">chr1\nACGTacgtNNRY\n>chr2\nAAA\n", proteins=">p1\nMAA*\n>p2\nM*A\n>p3\nMA**\n>p4\nmAA\n"):
    paths = [tmp_path / name for name in ("genome.fa", "proteins.fa", "genes.gff3")]
    paths[0].write_text(genome)
    paths[1].write_text(proteins)
    paths[2].write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t1\t8\t.\t+\t.\tID=g%2C1\n"
        "chr1\tx\tmRNA\t1\t6\t.\t+\t.\tParent=g%2C1;ID=p1\n"
        "chr1\tx\ttranscript\t1\t8\t.\t+\t.\tID=p2;Parent=g%2C1\n"
        "chr1\tx\tmRNA\t2\t8\t.\t-\t.\tID=p3;Parent=g2,g3\n"
        "chr1\tx\tmRNA\t3\t8\t.\t+\t.\tID=p4\n"
        "chr1\tx\tCDS\t1\t6\t.\t+\t0\tParent=p1\n"
    )
    return paths


def test_measurements_have_explicit_denominators_and_no_quality_grade(tmp_path):
    inputs = _inputs(tmp_path)
    report = collect_input_qc(*inputs)
    assert report["schema_version"] == 1
    genome = report["genome"]
    assert genome["n_sequences"] == 2
    assert genome["total_sequence_characters"] == 15
    assert genome["n50_bp"] == 12
    assert genome["n_bases"] == 2
    assert genome["pct_n_bases"] == pytest.approx(100 * 2 / 15)
    assert genome["lowercase_acgt_bases"] == 4
    assert genome["pct_lowercase_acgt"] == pytest.approx(100 * 4 / 11)
    assert genome["other_sequence_characters"] == 2
    proteins = report["proteins"]
    assert proteins["n_with_internal_stop"] == 2  # M*A and MA**; MAA* is terminal only.
    assert proteins["pct_with_internal_stop"] == 50
    assert proteins["n_with_terminal_stop"] == 2
    assert proteins["n_starting_with_m"] == 4
    assert proteins["mean_length"] == 3.5
    assert not {"status", "grade", "passed", "orf_complete"} & report.keys()
    assert not {"orf_complete", "hard_masked", "repeat_masked"} & genome.keys()
    assert any("does not establish why" in note for note in report["notes"])
    assert "not an ORF-completeness assessment" in report["definitions"]["starting_with_m"]


def test_gff_parent_counts_do_not_infer_loci_or_collapse_encoded_ids(tmp_path):
    report = collect_input_qc(*_inputs(tmp_path))["gff3"]
    assert report["n_features"] == 6
    assert report["n_gene_features"] == report["n_distinct_gene_ids"] == 1
    assert report["n_transcript_features"] == report["n_distinct_transcript_ids"] == 4
    assert report["n_transcript_features_without_parent"] == 1
    assert report["n_transcript_features_with_multiple_parents"] == 1
    assert report["n_parent_ids_referenced_by_transcripts"] == 3  # g,1 is one escaped ID.
    assert report["n_parent_ids_with_multiple_transcripts"] == 1
    assert report["n_declared_genes_with_multiple_transcripts"] == 1
    assert "n_loci" not in report


def test_read_only_deterministic_report_and_actual_byte_hashes(tmp_path):
    inputs = _inputs(tmp_path)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    first = collect_input_qc(*inputs)
    assert collect_input_qc(*inputs) == first
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
    for section, path in zip(("genome", "proteins", "gff3"), inputs):
        assert first[section]["source"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_wrapping_case_and_filename_do_not_change_noncase_measurements(tmp_path):
    inputs = _inputs(tmp_path)
    first = collect_input_qc(*inputs)
    inputs[0].write_text(">chr1 description\r\nACGT\r\nacgt\r\nNNRY\r\n>chr2\r\nAAA\r\n")
    inputs[1].write_text(">p1\nM\nAA*\n>p2\nM*\nA\n>p3\nMA*\n*\n>p4\nm\nAA\n")
    renamed = tmp_path / "unmasked-but-name-does-not-classify.fa"
    inputs[0].rename(renamed)
    second = collect_input_qc(renamed, *inputs[1:])
    for section in ("genome", "proteins", "gff3"):
        assert {k: v for k, v in first[section].items() if k != "source"} == {
            k: v for k, v in second[section].items() if k != "source"
        }
    assert first["notes"] == second["notes"]


def test_all_ambiguous_genome_has_unknown_lowercase_fraction(tmp_path):
    paths = _inputs(tmp_path, genome=">chr1\nNNNNNNNNNNNN\n")
    report = collect_input_qc(*paths)
    assert report["genome"]["pct_n_bases"] == 100
    assert report["genome"]["pct_lowercase_acgt"] is None
    assert report["genome"]["n50_bp"] == 12
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("bad_proteins, message", [
    (">p1\nAAA\n>p1\nMAA\n", "Duplicate FASTA"),
    (">p1\nAA-A\n", "Invalid unaligned FASTA"),
    ("MAA\n", "Expected uncompressed FASTA"),
    (">p1\n", "nonempty sequences"),
])
def test_standalone_reuses_fasta_validation_without_models(tmp_path, bad_proteins, message):
    paths = _inputs(tmp_path, proteins=bad_proteins)
    with pytest.raises(ValueError, match=message):
        collect_input_qc(*paths)


@pytest.mark.parametrize("line, message", [
    ("chr1\tx\tgene\t1\t99\t.\t+\t.\tID=g\n", "outside genome"),
    ("chr1\tx\tgene\tone\t8\t.\t+\t.\tID=g\n", "Non-integer"),
    ("chr1\tx\tmRNA\t1\t8\t.\t.\t.\tID=t\n", "unique ID and"),
    ("chr1\tx\tgene\n", "9 fields"),
    ("##FASTA\n>chr1\nAA\n", "embedded FASTA"),
])
def test_standalone_gff_errors_are_explicit(tmp_path, line, message):
    paths = _inputs(tmp_path)
    paths[2].write_text(line)
    with pytest.raises(ValueError, match=message):
        collect_input_qc(*paths)


def test_new_json_writer_does_not_overwrite_inputs_or_follow_symlinks(tmp_path):
    paths = _inputs(tmp_path)
    report = collect_input_qc(*paths)
    out = tmp_path / "input_qc.json"
    write_input_qc(report, out)
    assert json.loads(out.read_text()) == report
    with pytest.raises(FileExistsError):
        write_input_qc(report, out)
    with pytest.raises(FileExistsError):
        write_input_qc(report, paths[0])
    link = tmp_path / "linked.json"
    link.symlink_to(out)
    with pytest.raises(ValueError, match="symlinks"):
        write_input_qc(report, link)
    with pytest.raises(FileNotFoundError):
        write_input_qc(report, tmp_path / "absent" / "input_qc.json")


def test_stage0_adds_checked_qc_without_changing_primary_exports(tmp_path):
    from funlr.config import FunlrConfig
    from funlr.core.context import RunPaths
    from funlr.stages import stage0_validate

    genome, proteins, gff = _inputs(tmp_path)
    Path(str(genome) + ".fai").write_text("chr1\t12\t0\t12\t13\nchr2\t3\t0\t3\t4\n")
    config = FunlrConfig({"inputs": {"genome": str(genome), "proteins": str(proteins), "gff3": str(gff)}})
    runner = SimpleNamespace(inventory={}, tools={"samtools": "samtools"}, has_tool=lambda name: True,
                             resolve=lambda name: name, tool_version=lambda name: "fixture")
    ctx = SimpleNamespace(config=config, runner=runner, paths=RunPaths(tmp_path / "run"), dry_run=False,
                          stage_logger=lambda *args: logging.getLogger("input-qc-stage-test"))
    outputs = stage0_validate.run(ctx)
    stage0 = ctx.paths.stage_dir(0)
    qc_path = stage0 / "input_qc.json"
    assert qc_path in outputs
    assert json.loads(qc_path.read_text())["proteins"]["n_with_internal_stop"] == 2
    assert (stage0 / "proteins_clean.faa").read_text() == proteins.read_text()
    assert (stage0 / "protein_lengths.tsv").read_text() == "protein_id\tprotein_length\np1\t4\np2\t3\np3\t4\np4\t3\n"
    assert (stage0 / "contig_lengths.tsv").read_text() == "chr1\t12\nchr2\t3\n"
