"""Stage 0: post-OrthoFiller master table and input provenance.

Indexes the runner's private genome copy, reconciles transcript/CDS protein
IDs, and retains header-defined eggNOG PFAM annotations for enrichment.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError
from funlr.core.tools import require_tools
from funlr.parsers.eggnog import parse_eggnog_light
from funlr.parsers.fasta import fasta_iter, write_fasta
from funlr.parsers.gff import read_transcript_coords
from funlr.stages._inputs import reconcile_coordinates, build_updated_master_table

STAGE_NUMBER = 0
STAGE_NAME = "validate"
REQUIRED_TOOLS = ["samtools"]


# ---------------------------------------------------------------------------
# Pure helpers (plain arguments; unit-testable without a RunContext)
# ---------------------------------------------------------------------------


def build_contig_lengths(fai_path: str | Path, out_path: str | Path) -> int:
    """cut -f1,2 of a .fai index -> 2-column TSV, NO header (legacy 0.1)."""
    n = 0
    with open(fai_path) as f, open(out_path, "w") as o:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            o.write(f"{fields[0]}\t{fields[1]}\n")
            n += 1
    return n


def build_protein_files(
    proteins_path: str | Path, out_faa: str | Path, out_len: str | Path
) -> int:
    """Clean protein FASTA (IDs only, 60-col wrap) + lengths TSV (legacy 0.2)."""
    n = 0
    with open(proteins_path) as f, open(out_faa, "w") as o, open(out_len, "w") as l:
        l.write("protein_id\tprotein_length\n")
        for pid, seq in fasta_iter(f):
            write_fasta(iter([(pid, seq)]), o, wrap=60)
            l.write(f"{pid}\t{len(seq)}\n")
            n += 1
    return n


def build_master_table(
    coords_path: str | Path,
    lens_path: str | Path,
    eggnog_path: str | Path | None,
    eggnog_light_out: str | Path,
    master_out: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """eggnog light table + coords.merge(lens).merge(eggnog) (legacy 0.4)."""
    lens = pd.read_csv(lens_path, sep="\t")
    coords = pd.read_csv(coords_path, sep="\t")
    eg = parse_eggnog_light(eggnog_path)
    eg.to_csv(eggnog_light_out, sep="\t", index=False)
    master = coords.merge(lens, on="protein_id", how="left").merge(
        eg, on="protein_id", how="left"
    )
    master.to_csv(master_out, sep="\t", index=False)
    return master, eg


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------


def _require_nonempty(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise StageError(f"Required input missing/empty: {path} ({label})")


def run(ctx: RunContext) -> list[Path]:
    log = ctx.stage_logger(STAGE_NUMBER, STAGE_NAME)
    log.info("=== Stage 0: Build master table ===")
    log.info("Start: %s", time.ctime())
    require_tools(ctx.runner, REQUIRED_TOOLS, log)

    cfg = ctx.config
    genome = Path(cfg.get("inputs", "genome") or "")
    proteins = Path(cfg.get("inputs", "proteins") or "")
    gff3 = Path(cfg.get("inputs", "gff3") or "")
    eggnog = cfg.get("inputs", "eggnog")
    annotations = cfg.get("inputs", "annotations")

    for p, label in ((genome, "inputs.genome"), (proteins, "inputs.proteins"), (gff3, "inputs.gff3")):
        _require_nonempty(p, label)
    for path, label in ((eggnog, "inputs.eggnog"), (annotations, "inputs.annotations")):
        if path:
            _require_nonempty(Path(path), label)

    stage0 = ctx.paths.stage_dir(STAGE_NUMBER)
    outputs: list[Path] = []

    # The orchestrator supplies a private genome copy under work/. Indexing
    # there keeps source inputs read-only and avoids symlinks in saved outputs.
    fai = Path(f"{genome}.fai")
    if not fai.is_file():
        log.info("[0.1] samtools faidx %s", genome)
        ctx.runner.run([ctx.runner.resolve("samtools"), "faidx", str(genome)], check=True)

    if ctx.dry_run:
        log.info("PLANNED: stage 0 outputs would be built from %s.fai", genome)
        return []

    contig_lengths = stage0 / "contig_lengths.tsv"
    n_contigs = build_contig_lengths(fai, contig_lengths)
    log.info("[0.1] contig_lengths.tsv: %d contigs", n_contigs)
    outputs.append(contig_lengths)

    # 0.2) Clean protein FASTA + lengths
    log.info("[0.2] Protein FASTA clean + lengths")
    out_faa = stage0 / "proteins_clean.faa"
    out_len = stage0 / "protein_lengths.tsv"
    n_prot = build_protein_files(proteins, out_faa, out_len)
    log.info("[0.2] %d proteins", n_prot)
    outputs += [out_faa, out_len]

    # 0.3) Reconcile mRNA/CDS attributes against the actual FASTA universe.
    log.info("[0.3] Reconcile transcript/CDS protein IDs")
    coords, mapping = reconcile_coordinates(gff3, out_faa)
    coords_out = stage0 / "id_coords.tsv"
    coords.to_csv(coords_out, sep="\t", index=False)
    log.info("Parsed %d transcript entries", len(coords))
    outputs.append(coords_out)
    mapping_out = stage0 / "id_mapping_report.tsv"
    pd.DataFrame(mapping.items(), columns=["metric", "value"]).to_csv(mapping_out, sep="\t", index=False)
    outputs.append(mapping_out)

    # 0.4) Parse eggNOG + merge (fail-soft if missing)
    log.info("[0.4] Parse eggNOG + merge (fail-soft if missing)")
    if not eggnog or not (Path(eggnog).is_file() and Path(eggnog).stat().st_size > 0):
        log.warning("eggNOG annotations missing/empty (continuing): %s", eggnog)
    eggnog_light = stage0 / "eggnog_light.tsv"
    master_out = stage0 / "master_table.tsv"
    master, eg = build_updated_master_table(coords_out, out_len, eggnog, eggnog_light, master_out)
    log.info("Master table rows: %d", master.shape[0])
    log.info("eggNOG light rows: %d", eg.shape[0])
    outputs += [eggnog_light, master_out]

    # Paths plus dereferenced file sizes; root provenance records content hashes.
    manifest_rows = []
    for label, path in (("GENOME_FILE", genome), ("PROTEINS_FILE", proteins),
                        ("GFF_FILE", gff3), ("ANNOTATIONS_FILE", annotations), ("EGGNOG_FILE", eggnog)):
        if path:
            p = Path(path)
            manifest_rows.append((label, str(p), str(p.resolve()), p.stat().st_size))
    input_manifest = stage0 / "input_manifest.tsv"
    pd.DataFrame(manifest_rows, columns=["label", "path", "resolved_path", "size_bytes"]).to_csv(input_manifest, sep="\t", index=False)
    software_rows = []
    for key, info in ctx.runner.inventory.items():
        software_rows.append(("tool", key, info.get("version_output", "unknown"),
                              info.get("path", ctx.runner.tools.get(key, key)), "detected"))
    for key in ("pfam", "nbd_hmms", "custom_hmms", "asm_hmms", "effector_hmms", "sensor_hmms"):
        path = cfg.get("databases", key)
        if path:
            software_rows.append(("hmm_db", key, str(path), str(Path(path).resolve()), "configured"))
    software_manifest = stage0 / "software_manifest.tsv"
    pd.DataFrame(software_rows, columns=["kind", "name", "value", "resolved_value", "status"]).to_csv(software_manifest, sep="\t", index=False)
    outputs += [input_manifest, software_manifest]

    # run_info.txt (stage root; additive per SPEC deviation 3 — legacy stage 0
    # wrote none, so fields mirror the stage-1 layout plus input paths).
    run_info = stage0 / "run_info.txt"
    run_info.write_text(
        f"stage: {STAGE_NUMBER}\n"
        f"date: {time.ctime()}\n"
        f"genome: {cfg.get('sample', 'species_id') or cfg.get('sample', 'sample_id') or 'unknown'}\n"
        f"assembly_version: {cfg.get('sample', 'assembly_version') or 'unknown'}\n"
        f"genome_file: {genome}\n"
        f"proteins_file: {proteins}\n"
        f"gff_file: {gff3}\n"
        f"eggnog_file: {eggnog or ''}\n"
    )
    outputs.append(run_info)

    # Descriptive measurements are additive and do not change master-table
    # values or select candidates. Original protein stop characters are kept.
    from funlr.input_qc import collect_input_qc, write_input_qc
    input_qc = collect_input_qc(genome, proteins, gff3, validated=True)
    qc_path = write_input_qc(input_qc, stage0 / "input_qc.json")
    outputs.append(qc_path)
    for note in input_qc["notes"]:
        log.info("Input measurements: %s", note)

    log.info("Done: %s", time.ctime())
    return outputs
