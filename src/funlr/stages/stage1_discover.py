"""Stage 1 (discover): NBD hmmscan + candidate-id union.

The stage entry point runs updated strict/relaxed NBD discovery, Pfam NBD
scanning, enrichment, and a configured description union. Earlier pure parser
helpers remain available for compatibility with existing callers.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path

import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError
from funlr.core.tools import ensure_pressed, is_pressed, require_tools
from funlr.parsers.hmmer import parse_domtblout


def _ensure_pressed(ctx: RunContext, log, hmm_path: str | Path) -> Path:
    """ensure_pressed, but dry-run safe: log the hmmpress plan and move on."""
    if ctx.dry_run and not is_pressed(hmm_path):
        ctx.runner.run([ctx.runner.resolve("hmmpress"), "-f", str(hmm_path)], check=True)
        return Path(hmm_path)
    return ensure_pressed(hmm_path, ctx.paths.pressed_dir, ctx.runner, log)

STAGE_NUMBER = 1
STAGE_NAME = "discover"
REQUIRED_TOOLS = ["hmmscan", "hmmpress", "seqkit", "hmmfetch"]

# Legacy keyword list for the eggNOG priority list (exact, legacy lines 450-455).
EGGNOG_PRIORITY_KEYWORDS = [
    "NACHT", "NB-ARC", "STAND", "NOD-like",
    "WD40", "ANK", "TPR", "HEAT", "Kelch",
    "HET", "gasdermin", "HeLo", "MLKL", "Goodbye",
    "inflammasome", "immune",
]

_EGGNOG_PRIORITY_PATTERN = re.compile(
    "|".join(re.escape(k) for k in EGGNOG_PRIORITY_KEYWORDS), re.IGNORECASE
)


# ---------------------------------------------------------------------------
# Pure helpers (plain arguments; unit-testable without a RunContext)
# ---------------------------------------------------------------------------


def write_nbd_parse_outputs(
    domtbl_path: str | Path,
    eval_nbd: float,
    min_nbd_ali_len: int,
    hits_out: str | Path,
    ids_out: str | Path,
    report_out: str | Path,
) -> tuple[pd.DataFrame, dict]:
    """Parse NBD domtblout -> hits TSV, candidate-id list, parse report.

    Legacy 1.2: hard filter on i-Evalue and aligned length; report key order
    is exactly as the legacy script wrote it.
    """
    df, rep = parse_domtblout(domtbl_path, eval_nbd, min_nbd_ali_len, source="nbd")
    hits = df.rename(
        columns={"domain": "nbd_hmm", "start": "ali_start", "end": "ali_end"}
    )[["protein_id", "nbd_hmm", "dom_i_evalue", "ali_start", "ali_end", "ali_len"]]
    hits.to_csv(hits_out, sep="\t", index=False)

    cand = hits[["protein_id"]].drop_duplicates()
    cand.to_csv(ids_out, index=False, header=False)

    with open(report_out, "w") as out:
        out.write(f"total_domtbl_lines\t{rep['total_domtbl_lines']}\n")
        out.write(f"kept_after_filters\t{rep['kept_after_filters']}\n")
        out.write(f"filtered_high_i_eval\t{rep['filtered_high_i_eval']}\n")
        out.write(f"filtered_short_align\t{rep['filtered_short_align']}\n")
        out.write(f"unique_candidate_proteins\t{cand.shape[0]}\n")
    return hits, rep


def write_eggnog_priority_ids(
    master_path: str | Path, out_path: str | Path
) -> pd.DataFrame:
    """eggNOG keyword priority list (legacy 1.3; non-gating, fail-soft).

    Case-insensitive keyword regex over Description/Preferred_name/GOs columns
    of master_table.tsv. Missing/empty master (or missing protein_id column)
    yields an empty output file.
    """
    out_path = Path(out_path)
    if (not os.path.exists(master_path)) or os.path.getsize(master_path) == 0:
        out_path.write_text("")
        return pd.DataFrame(columns=["protein_id"])

    master = pd.read_csv(master_path, sep="\t")
    if "protein_id" not in master.columns:
        out_path.write_text("")
        return pd.DataFrame(columns=["protein_id"])

    cols = [c for c in ["Description", "Preferred_name", "GOs"] if c in master.columns]

    def row_hit(r) -> bool:
        for c in cols:
            v = r.get(c, "")
            if isinstance(v, str) and _EGGNOG_PRIORITY_PATTERN.search(v):
                return True
        return False

    if cols:
        prio = master[master.apply(row_hit, axis=1)][["protein_id"]].dropna().drop_duplicates()
    else:
        prio = master[["protein_id"]].head(0)

    prio.to_csv(out_path, index=False, header=False)
    return prio


def write_union_ids(
    id_files: list[str | Path], out_path: str | Path
) -> list[str]:
    """cat id files | drop empty lines | sort -u (legacy 1.4)."""
    ids: set[str] = set()
    for p in id_files:
        p = Path(p)
        if not p.is_file():
            continue
        with p.open() as f:
            for line in f:
                line = line.strip()
                if line:
                    ids.add(line)
    ordered = sorted(ids)
    with open(out_path, "w") as o:
        for pid in ordered:
            o.write(pid + "\n")
    return ordered


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------


from ._discovery_run import run  # Native orchestration; compatibility helpers remain importable.
