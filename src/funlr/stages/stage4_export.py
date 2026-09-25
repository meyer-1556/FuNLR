"""Stage 4: reference tier-table autodiscovery and checked literal FASTA export.

TIER_SPECS retains the earlier naming map for baseline tests; production exports
use each current tier filename plus tiered_candidates and rescue_priority.
"""

from __future__ import annotations

import csv
import os
import re
from collections import Counter
import time
from collections import OrderedDict
from pathlib import Path

from funlr.core.context import RunContext
from funlr.core.errors import StageError
from funlr.core.tools import require_tools

STAGE_NUMBER = 4
STAGE_NAME = "export"
REQUIRED_TOOLS = ["seqkit"]

# Tier label -> (stage3 tables TSV name, output FASTA name), in legacy order
# (legacy lines 1333-1358, 1371).
TIER_SPECS: "OrderedDict[str, tuple[str, str]]" = OrderedDict(
    [
        ("tier1A", ("TIER_1A_HIGH_CONFIDENCE.tsv", "tier1A_high_confidence.faa")),
        ("tier1B", ("TIER_1B_NEEDS_REVIEW.tsv", "tier1B_needs_review.faa")),
        ("tier2A", ("TIER_2A_HIGH_PRIORITY_RESCUE.tsv", "tier2A_high_priority_rescue.faa")),
        ("tier2B", ("TIER_2B_RESCUE_CANDIDATE.tsv", "tier2B_rescue_candidate.faa")),
        ("tier3", ("TIER_3_ARCHITECTURAL_VARIANT.tsv", "tier3_architectural_variant.faa")),
        ("tier4A", ("TIER_4A_REPEAT_ONLY_NO_NBD.tsv", "tier4A_repeat_only_no_nbd.faa")),
        ("tier4B", ("TIER_4B_LIKELY_STAND_HOUSEKEEPING.tsv", "tier4B_likely_stand_housekeeping.faa")),
        ("tier4C", ("TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL.tsv", "tier4C_no_nbd_no_repeat_low_signal.faa")),
    ]
)

ALL_TIERED_TSV = "tiered_candidates.tsv"
ALL_TIERED_OUT = "all_tiered_candidates.faa"


# ---------------------------------------------------------------------------
# Pure helpers (legacy extract_ids / export_faa / summary)
# ---------------------------------------------------------------------------


def _nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def extract_ids(in_tsv: str | Path, out_ids: str | Path) -> Path:
    """Read first-column IDs, skipping only an actual protein_id header."""
    in_tsv, out_ids = Path(in_tsv), Path(out_ids)
    ids = []
    if _nonempty(in_tsv):
        with in_tsv.open() as stream:
            for line_number, line in enumerate(stream):
                field = line.rstrip("\n\r").split("\t")[0]
                if line_number == 0 and field == "protein_id":
                    continue
                if field.strip():
                    ids.append(field)
    out_ids.write_text("".join(pid + "\n" for pid in ids))
    return out_ids


def extract_rescue_ids(in_tsv, out_ids, require_nbd_ok=False):
    """Honor nbd_ok when supplied; otherwise use confidence, then keep all."""
    in_tsv, out_ids = Path(in_tsv), Path(out_ids)
    if not require_nbd_ok or not _nonempty(in_tsv):
        return extract_ids(in_tsv, out_ids)
    ids = []
    with in_tsv.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = reader.fieldnames or []
        if "protein_id" not in columns:
            raise StageError(f"Rescue priority table lacks protein_id header: {in_tsv}")
        for row in reader:
            pid = row.get("protein_id", "")
            if not pid:
                continue
            if "nbd_ok" in columns:
                keep = str(row.get("nbd_ok", "")).strip() in {"True", "TRUE", "1"}
            elif "nbd_confidence" in columns:
                keep = str(row.get("nbd_confidence", "")).strip().upper() in {"HIGH", "MEDIUM"}
            else:
                keep = True
            if keep:
                ids.append(pid)
    out_ids.write_text("".join(pid + "\n" for pid in ids))
    return out_ids


def export_faa(ctx: RunContext, ids: str | Path, prot: str | Path, out_faa: str | Path) -> Path:
    """seqkit grep of an ID list; empty ID list -> empty FASTA, no call.

    Empty ID lists produce an empty FASTA without invoking SeqKit. Real tool
    failures propagate, so a failed export cannot become a successful stage.
    """
    ids = Path(ids)
    out_faa = Path(out_faa)
    out_faa.write_text("")  # : > "${out_faa}"

    if not _nonempty(ids):
        return out_faa

    ctx.runner.run(
        [ctx.runner.resolve("seqkit"), "grep", "-f", str(ids), str(prot)],
        stdout_path=out_faa,
        check=True,
    )
    requested = set(ids.read_text().splitlines())
    observed = Counter()
    with out_faa.open() as stream:
        for line in stream:
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields:
                    raise StageError(f"Malformed empty FASTA identifier in {out_faa}")
                observed[fields[0]] += 1
    missing = requested - set(observed)
    unexpected = set(observed) - requested
    duplicates = sorted(pid for pid, count in observed.items() if count != 1)
    if missing or unexpected or duplicates:
        raise StageError(
            f"FASTA extraction did not return exactly the requested IDs: "
            f"missing={sorted(missing)[:10]}, unexpected={sorted(unexpected)[:10]}, "
            f"duplicates={duplicates[:10]}; output={out_faa}"
        )
    return out_faa


def write_summary(stage4: str | Path, out_path: str | Path) -> Path:
    """fasta_export_summary.tsv: file/bytes/seqs for every .faa (legacy 1403-1412)."""
    stage4 = Path(stage4)
    rows: list[tuple[str, int, int]] = []
    for f in sorted(stage4.glob("*.faa"), key=lambda path: path.name.casefold()):  # bash glob order = sorted
        bytes_ = os.path.getsize(f)
        seqs = 0
        with open(f) as fh:
            for line in fh:
                if line.startswith(">"):
                    seqs += 1
        rows.append((f.name, bytes_, seqs))
    with open(out_path, "w") as out:
        out.write("file\tbytes\tseqs\n")
        for name, bytes_, seqs in rows:
            out.write(f"{name}\t{bytes_}\t{seqs}\n")
    return Path(out_path)


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> list[Path]:
    """Export all discovered reference tier tables and the two combined handoffs."""
    log = ctx.stage_logger(STAGE_NUMBER, STAGE_NAME)
    if ctx.dry_run:
        log.info("PLANNED: dynamic tier FASTA exports and confidence-filtered rescue handoff")
        return []
    require_tools(ctx.runner, REQUIRED_TOOLS, log)
    stage0 = ctx.paths.results_dir / "stage0"
    stage3 = ctx.paths.results_dir / "stage3"
    stage4 = ctx.paths.stage_dir(STAGE_NUMBER)
    proteins = stage0 / "proteins_clean.faa"
    tables = stage3 / "tables"
    if not _nonempty(proteins):
        raise StageError(f"Missing/empty proteins FASTA: {proteins}")
    if not tables.is_dir():
        raise StageError(f"Stage 3 tables directory not found: {tables}")
    required = [stage0 / "input_manifest.tsv", stage0 / "software_manifest.tsv",
                tables / "tiered_candidates.tsv", tables / "rescue_priority.tsv"]
    for path in required:
        if not _nonempty(path):
            raise StageError(f"Required Stage 4 input missing/empty: {path}")
    settings = ctx.config.scientific_settings()
    require_nbd_ok = bool(settings.get("REQUIRE_NBD_OK_IN_RESCUE_FASTA", 0))
    manifest = stage4 / "fasta_export_manifest.tsv"
    outputs = [manifest]
    def natural_key(path):
        return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", path.name)]
    tiers = sorted(tables.glob("TIER_*.tsv"), key=natural_key)
    sources = [(path.stem, path, False) for path in tiers if _nonempty(path)]
    sources.extend([
        ("tiered_candidates", tables / "tiered_candidates.tsv", False),
        ("rescue_priority", tables / "rescue_priority.tsv", require_nbd_ok),
    ])
    with manifest.open("w") as stream:
        stream.write("label\tsource_tsv\tids_file\tfasta_file\n")
        for label, source, filtered in sources:
            ids, fasta = stage4 / f"{label}.ids", stage4 / f"{label}.faa"
            if label == "rescue_priority":
                extract_rescue_ids(source, ids, filtered)
            else:
                extract_ids(source, ids)
            export_faa(ctx, ids, proteins, fasta)
            stream.write(f"{label}\t{source}\t{ids}\t{fasta}\n")
            outputs.extend([ids, fasta])
    summary = write_summary(stage4, stage4 / "fasta_export_summary.tsv")
    outputs.append(summary)
    run_info = stage4 / "run_info.txt"
    info = {
        "stage": STAGE_NUMBER, "date": time.ctime(),
        "genome": ctx.config.get("sample", "species_id") or ctx.sample_id,
        "assembly_version": ctx.config.get("sample", "assembly_version") or "",
        "nlr_profile": settings.get("NLR_PROFILE", "ILLUMINA"),
        "input_manifest": stage0 / "input_manifest.tsv", "software_manifest": stage0 / "software_manifest.tsv",
        "stage3_tables_dir": tables, "proteins_fasta": proteins,
        "seqkit_executable": ctx.runner.resolve("seqkit"),
        "require_nbd_ok_in_rescue_fasta": int(require_nbd_ok),
        "manifest": manifest, "summary": summary,
        "notes": "Tier TSV autodiscovery; checked literal ID extraction; missing requested IDs fail the stage",
    }
    run_info.write_text("".join(f"{key}: {value}\n" for key, value in info.items()))
    outputs.append(run_info)
    log.info("Exported %d FASTA subsets", len(sources))
    return outputs
