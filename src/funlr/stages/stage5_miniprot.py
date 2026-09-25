"""Stage 5: rescue candidates with miniprot (port of 05_rescue_miniprot.sh).

Builds the rescue query FASTA from stage 4 tier exports, runs miniprot
against the genome, parses the resulting GFF3 into feature/locus tables, and
emits the list of queries still needing exonerate refinement.
"""

from __future__ import annotations

from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError
from funlr.core.tools import require_tools
from funlr.parsers.fasta import dedup_keep_first, fasta_ids
from funlr.parsers.gff import MINIPROT_FEATURE_COLUMNS, miniprot_locus_summary, read_miniprot_features

STAGE_NUMBER = 5
STAGE_NAME = "miniprot"
REQUIRED_TOOLS = ["miniprot"]

# Empty-frame column sets (legacy stage 5.2, lines 1692-1694).
SUMMARY_COLUMNS = ["protein_id", "scaffold", "start", "end", "strand", "n_features"]
COUNTS_COLUMNS = ["feature", "count"]

# Manifest rows, in legacy order (lines 1750-1762).
_MANIFEST_LABELS = [
    "query_faa_raw",
    "query_faa_dedup",
    "miniprot_gff3",
    "miniprot_log",
    "miniprot_features",
    "miniprot_summary",
    "miniprot_hit_counts",
    "exonerate_refine_ids",
    "query_counts",
    "run_info",
]


# ---------------------------------------------------------------------------
# Pure logic (no RunContext; unit-testable)
# ---------------------------------------------------------------------------

def concat_query_fastas(
    tier_fastas: "OrderedDict[str, Path]", out_path: Path, include_tier3: bool
) -> Path:
    """Concatenate tier FASTAs in legacy order (append only if exists+nonempty).

    Order: tier2A, tier2B, tier1B, then tier3 iff ``include_tier3``.
    Pure Python text append (legacy ``cat >>``, lines 1492-1508).
    """
    order = ["tier2A", "tier2B", "tier1B"] + (["tier3"] if include_tier3 else [])
    out_path = Path(out_path)
    with out_path.open("w") as out:
        for label in order:
            src = Path(tier_fastas[label])
            if src.is_file() and src.stat().st_size > 0:
                with src.open() as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), ""):
                        out.write(chunk)
    return out_path


def count_faa(path: Path) -> int:
    """Count FASTA records; missing file counts as 0 (legacy count_faa)."""
    n = 0
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
    except FileNotFoundError:
        return 0
    return n


def write_query_counts(files: "OrderedDict[str, Path]", out_path: Path) -> None:
    """Write summaries/query_counts.tsv (label/file/seqs rows, legacy lines 1576-1606)."""
    with open(out_path, "w") as out:
        out.write("label\tfile\tseqs\n")
        for label, path in files.items():
            out.write(f"{label}\t{path}\t{count_faa(path)}\n")


def compute_refine_ids(all_ids: Iterable[str], hit_ids: Iterable[str]) -> list[str]:
    """Sorted set of query ids with NO miniprot locus (legacy lines 1716-1746)."""
    return sorted(set(all_ids) - set(hit_ids))


def write_refine_ids(ids: Iterable[str], out_path: Path) -> None:
    with open(out_path, "w") as out:
        for pid in ids:
            out.write(pid + "\n")


def hit_ids_from_summary(summary_path: Path) -> set[str]:
    """protein_id set from miniprot_summary.tsv; empty on any problem (legacy 1726-1736)."""
    hit: set[str] = set()
    try:
        with open(summary_path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if "protein_id" in header:
                pid_i = header.index("protein_id")
                for line in fh:
                    if not line.strip():
                        continue
                    hit.add(line.rstrip("\n").split("\t")[pid_i])
    except Exception:  # noqa: BLE001 - legacy swallows all errors here
        hit = set()
    return hit


def write_manifest(paths: "OrderedDict[str, Path]", out_path: Path) -> None:
    with open(out_path, "w") as out:
        out.write("label\tpath\n")
        for label, path in paths.items():
            out.write(f"{label}\t{path}\n")


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(ctx: RunContext) -> list[Path]:
    from funlr.stages.rescue_v3 import run_miniprot
    return run_miniprot(ctx)
