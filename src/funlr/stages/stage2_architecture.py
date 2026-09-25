"""Stage 2 (architecture): Pfam/custom domain scans + architecture summary.

Faithful port of legacy ``02_domain_architecture.sh`` (legacy lines 499-756):
Pfam hmmscan of union candidates, optional custom HMM scan, hard-filtered
parse of both, per-protein domain architecture summary over ALL proteins in
the union FASTA (including zero-hit proteins).
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
from funlr.parsers.fasta import fasta_ids
from funlr.parsers.hmmer import parse_domtblout


def _ensure_pressed(ctx: RunContext, log, hmm_path: str | Path) -> Path:
    """ensure_pressed, but dry-run safe: log the hmmpress plan and move on."""
    if ctx.dry_run and not is_pressed(hmm_path):
        ctx.runner.run([ctx.runner.resolve("hmmpress"), "-f", str(hmm_path)], check=True)
        return Path(hmm_path)
    return ensure_pressed(hmm_path, ctx.paths.pressed_dir, ctx.runner, log)

STAGE_NUMBER = 2
STAGE_NAME = "architecture"
REQUIRED_TOOLS = ["hmmscan", "hmmpress"]

NBD_GROUPS = {"NACHT", "NB-ARC"}
REPEAT_GROUPS = {"WD40", "Ank", "TPR", "HEAT", "Kelch"}

# Legacy FP_PATTERNS (exact, legacy lines 683-687).
FP_PATTERNS = [
    r"\bSMC\b", r"ABC_tran|ABC_transporter", r"Septin", r"Dynamin",
    r"Helicase_C|DEAD|DEXDc", r"ResIII|Restriction", r"Toprim",
    r"FtsK_SpoIIIE", r"\bPkinase\b",
]
FP_RE = re.compile("|".join(FP_PATTERNS), re.I)


def group_domain(d: str) -> str:
    """Map a raw domain name to its architecture group (legacy, exact regexes)."""
    if re.search(r"\bNACHT\b", d, re.I):
        return "NACHT"
    if re.search(r"\bNB-ARC\b", d, re.I):
        return "NB-ARC"
    if re.search(r"Ank", d, re.I):
        return "Ank"
    if re.search(r"\bTPR\b", d, re.I):
        return "TPR"
    if re.search(r"WD40|WD_40|WD-repeat", d, re.I):
        return "WD40"
    if re.search(r"\bHEAT\b", d, re.I):
        return "HEAT"
    if re.search(r"Kelch", d, re.I):
        return "Kelch"
    if re.search(r"LRR", d, re.I):
        return "LRR"
    return d


# ---------------------------------------------------------------------------
# Pure helpers (plain arguments; unit-testable without a RunContext)
# ---------------------------------------------------------------------------


def summarize_architecture(
    all_prots: list[str],
    hits: pd.DataFrame,
    stage1_nbd: set[str],
) -> pd.DataFrame:
    """Per-protein architecture summary over ALL proteins (legacy 2.3).

    Zero-hit proteins are kept; ``has_nbd = has_nbd_stage1 | has_nbd_pfam``.
    Column names/order exactly as legacy.
    """
    df = hits.copy()
    if not df.empty:
        df["domain_group"] = df["domain"].apply(group_domain)
    else:
        df["domain_group"] = pd.Series(dtype=str)

    hit_map: dict[str, pd.DataFrame] = {}
    if not df.empty:
        for pid, g in df.groupby("protein_id"):
            hit_map[pid] = g.sort_values(
                ["start", "end", "dom_i_evalue"], ascending=[True, True, True]
            )

    summ = []
    for pid in all_prots:
        if pid in hit_map:
            g = hit_map[pid]
            doms_raw = g["domain"].tolist()
            doms_grp = g["domain_group"].tolist()
            has_nbd_pfam = any(d in NBD_GROUPS for d in doms_grp)
            has_rep = any(d in REPEAT_GROUPS for d in doms_grp)
            has_lrr = any(d == "LRR" for d in doms_grp)
            has_fp = any(FP_RE.search(d) for d in doms_raw)
            order_ok = True
            if has_nbd_pfam and has_rep:
                nbd_min = g[g["domain_group"].isin(NBD_GROUPS)]["start"].min()
                rep_min = g[g["domain_group"].isin(REPEAT_GROUPS)]["start"].min()
                order_ok = nbd_min <= rep_min
            summ.append(
                {
                    "protein_id": pid,
                    "domains_raw": ";".join(doms_raw),
                    "domains_grouped": ";".join(doms_grp),
                    "has_nbd_pfam": bool(has_nbd_pfam),
                    "has_repeat": bool(has_rep),
                    "order_ok": bool(order_ok),
                    "flag_fp_domains": bool(has_fp),
                    "flag_rare_lrr": bool(has_lrr),
                    "has_nbd_stage1": bool(pid in stage1_nbd),
                }
            )
        else:
            # No Pfam/custom domains; keep row so Stage 3 never silently drops it
            summ.append(
                {
                    "protein_id": pid,
                    "domains_raw": "",
                    "domains_grouped": "",
                    "has_nbd_pfam": False,
                    "has_repeat": False,
                    "order_ok": True,
                    "flag_fp_domains": False,
                    "flag_rare_lrr": False,
                    "has_nbd_stage1": bool(pid in stage1_nbd),
                }
            )

    out = pd.DataFrame(summ, columns=[
        "protein_id", "domains_raw", "domains_grouped", "has_nbd_pfam",
        "has_repeat", "order_ok", "flag_fp_domains", "flag_rare_lrr",
        "has_nbd_stage1",
    ])
    # Final NBD call used downstream: stage1 OR pfam
    out["has_nbd"] = out["has_nbd_stage1"] | out["has_nbd_pfam"]
    return out


def parse_domain_hits(
    domtbl_paths: dict[str, str | Path],
    eval_use: float,
    min_dom_ali_len: int,
    allowed_proteins: set[str],
    all_hits_out: str | Path,
    report_out: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse pfam+custom domtblout files -> all_domain_hits + parse report.

    Legacy 2.3: concat both sources, restrict to proteins in the union FASTA,
    and write a two-row parse_filter_report.tsv.
    """
    frames: list[pd.DataFrame] = []
    reports: list[dict] = []
    for source, path in domtbl_paths.items():
        df, rep = parse_domtblout(path, eval_use, min_dom_ali_len, source=source)
        frames.append(df)
        reports.append(rep)

    nonempty = [f for f in frames if not f.empty]
    if nonempty:
        df = pd.concat(nonempty, ignore_index=True)
    else:
        df = pd.DataFrame(columns=frames[0].columns)
    df = df[df["protein_id"].isin(allowed_proteins)].copy()  # safety
    df.to_csv(all_hits_out, sep="\t", index=False)

    rep_df = pd.DataFrame(reports)
    rep_df.to_csv(report_out, sep="\t", index=False)
    return df, rep_df


def read_id_set(path: str | Path) -> set[str]:
    """Read a one-id-per-line file into a set (legacy stage1_nbd.ids reader)."""
    ids: set[str] = set()
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as f:
            for line in f:
                ids.add(line.strip())
    return ids


def write_sorted_ids(ids: set[str], out_path: str | Path) -> None:
    """sort -u of an id file (legacy stage1_nbd.ids writer)."""
    with open(out_path, "w") as o:
        for pid in sorted(ids):
            o.write(pid + "\n")


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------


from ._architecture_run import run  # Native orchestration; compatibility helpers remain importable.
