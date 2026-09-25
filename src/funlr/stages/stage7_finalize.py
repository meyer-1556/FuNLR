"""Stage 7: final report (port of 07_finalize_report.sh).

Merges miniprot/exonerate rescue loci onto the tiered candidate table, writes
the final report, strict candidate set, tier/QC summaries, and BED tracks —
each into results/stage7/ AND final_results/.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError

STAGE_NUMBER = 7
STAGE_NAME = "finalize"
REQUIRED_TOOLS: list[str] = []

# Strict tiers (legacy lines 2276-2282).
STRICT_TIERS = {
    "TIER_1A_HIGH_CONFIDENCE",
    "TIER_1B_NEEDS_REVIEW",
    "TIER_2A_HIGH_PRIORITY_RESCUE",
    "TIER_2B_RESCUE_CANDIDATE",
    "TIER_3_ARCHITECTURAL_VARIANT",
}

# Preferred output columns, kept if present (legacy lines 2253-2261).
WANT_COLUMNS = [
    "protein_id", "gene_id", "transcript_id",
    "scaffold", "start", "end", "strand", "protein_length",
    "tier", "rescue_priority", "flags",
    "domains_raw", "domains_grouped",
    "Description", "Preferred_name", "GOs", "KEGG_ko",
    "mp_scaffold", "mp_start", "mp_end", "mp_strand",
    "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features",
]

MP_NEED = ["protein_id", "scaffold", "start", "end", "strand"]
EXO_COLUMNS = ["protein_id", "exo_scaffold", "exo_start", "exo_end", "exo_strand", "n_features"]


# ---------------------------------------------------------------------------
# Pure logic (no RunContext; unit-testable)
# ---------------------------------------------------------------------------

def load_miniprot_summary(path: Path) -> pd.DataFrame:
    """Load stage-5 miniprot loci; missing/empty -> empty frame (legacy 2160-2178)."""
    path = Path(path)
    if path.exists() and path.stat().st_size > 0:
        mp = pd.read_csv(path, sep="\t")
        for c in MP_NEED:
            if c not in mp.columns:
                mp[c] = pd.NA
        mp = mp[MP_NEED].copy()
    else:
        mp = pd.DataFrame(columns=MP_NEED)
    return mp.rename(
        columns={"scaffold": "mp_scaffold", "start": "mp_start",
                 "end": "mp_end", "strand": "mp_strand"}
    )


_TARGET_RE = re.compile(r"(Target=)([^;\s]+)")
_PARENT_RE = re.compile(r"(Parent=)([^;\s]+)")
_FUNLR_QUERY_RE = re.compile(r"(?:^|;)FUNLR_query=([^;]+)")


def _parse_exonerate_merged_gff(gff_path: Path) -> pd.DataFrame:
    """Fallback: recover loci from merged GFF via Target=/Parent= (legacy 2200-2235)."""
    rows = []
    gff_path = Path(gff_path)
    if gff_path.exists() and gff_path.stat().st_size > 0:
        with gff_path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue
                chrom, _src, _ft, start, end, _score, strand, _phase, attr = parts
                pid = ""
                known = _FUNLR_QUERY_RE.search(attr)
                m = _TARGET_RE.search(attr)
                if known:
                    pid = unquote(known.group(1))
                elif m:
                    pid = m.group(2)
                else:
                    m2 = _PARENT_RE.search(attr)
                    if m2:
                        pid = m2.group(2)
                if not pid:
                    continue
                try:
                    rows.append((pid, chrom, int(start), int(end), strand))
                except ValueError:
                    continue
    if rows:
        ex = pd.DataFrame(
            rows, columns=["protein_id", "exo_scaffold", "exo_start", "exo_end", "exo_strand"]
        )
        ex = ex.groupby("protein_id", as_index=False).agg(
            exo_scaffold=("exo_scaffold", "first"),
            exo_start=("exo_start", "min"),
            exo_end=("exo_end", "max"),
            exo_strand=("exo_strand", "first"),
        )
        ex["n_features"] = pd.NA
        return ex
    return pd.DataFrame(columns=EXO_COLUMNS)


def load_exonerate_summary(summary_path: Path, merged_gff_path: Path) -> pd.DataFrame:
    """Load stage-6 exonerate loci; fall back to merged-GFF parse (legacy 2183-2238)."""
    summary_path = Path(summary_path)
    if summary_path.exists() and summary_path.stat().st_size > 0:
        ex = pd.read_csv(summary_path, sep="\t")
        for old, new in (("scaffold", "exo_scaffold"), ("start", "exo_start"),
                         ("end", "exo_end"), ("strand", "exo_strand")):
            if old in ex.columns:
                ex = ex.rename(columns={old: new})
        keep = [c for c in EXO_COLUMNS if c in ex.columns]
        return ex[keep].copy()
    return _parse_exonerate_merged_gff(merged_gff_path)


def merge_final_report(
    tiered: pd.DataFrame, mp: pd.DataFrame, ex: pd.DataFrame
) -> pd.DataFrame:
    """Left-merge loci onto tiered candidates; keep WANT columns (legacy 2243-2263)."""
    if "protein_id" not in tiered.columns:
        raise StageError("tiered_candidates.tsv missing required column: protein_id")
    out = tiered.copy()
    out = out.merge(mp, on="protein_id", how="left")
    out = out.merge(ex, on="protein_id", how="left")
    cols = [c for c in WANT_COLUMNS if c in out.columns]
    return out[cols].copy()


def has_non_nlr(x) -> bool:
    """Verbatim legacy NON_NLR_ANNOT flag test (lines 2290-2294)."""
    if pd.isna(x):
        return False
    s = str(x)
    return (
        "NON_NLR_ANNOT" in s.split("|")
        or "NON_NLR_ANNOT" in s.split(";")
        or "NON_NLR_ANNOT" in s.split(",")
        or "NON_NLR_ANNOT" in s
    )


def strict_candidates(out: pd.DataFrame) -> pd.DataFrame:
    """Strict set: STRICT_TIERS minus anything flagged NON_NLR_ANNOT (legacy 2283-2295)."""
    strict = out.copy()
    if "tier" in strict.columns:
        strict = strict[strict["tier"].isin(STRICT_TIERS)].copy()
    if "flags" in strict.columns:
        strict = strict.loc[~strict["flags"].apply(has_non_nlr).astype(bool)].copy()
    return strict


def tier_summary(out: pd.DataFrame) -> pd.DataFrame:
    """Tier counts (legacy lines 2305-2309)."""
    if "tier" in out.columns:
        summ = out["tier"].value_counts(dropna=False).reset_index()
        summ.columns = ["tier", "count"]
        return summ
    return pd.DataFrame({"tier": [], "count": []})


# -- QC helpers (legacy lines 2319-2439) ------------------------------------

def normalize_bool_series(s: pd.Series):
    """Accepts bool/int/str-ish; returns boolean series with NA -> False."""
    if s is None:
        return None
    if s.dtype == bool:
        return s.fillna(False)
    if np.issubdtype(s.dtype, np.number):
        return s.fillna(0) != 0
    # strings like TRUE/FALSE/yes/no/1/0
    x = s.astype(str).str.strip().str.lower()
    return x.isin(["true", "t", "yes", "y", "1", "ok", "pass"])


def find_col(df: pd.DataFrame, candidates: list[str]):
    """First matching column, with case-insensitive fallback."""
    for c in candidates:
        if c in df.columns:
            return c
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def infer_has_nbd_from_domains(df: pd.DataFrame):
    col = find_col(df, ["domains_grouped", "domains_raw"])
    if not col:
        return None
    d = df[col].fillna("").astype(str)
    return d.str.contains(r"\bNACHT\b|\bNB-ARC\b", regex=True)


def infer_has_repeat_from_domains(df: pd.DataFrame):
    col = find_col(df, ["domains_grouped", "domains_raw"])
    if not col:
        return None
    d = df[col].fillna("").astype(str)
    return d.str.contains(r"WD40|ANK|TPR|HEAT", regex=True)


def count_flag_tokens(flags: pd.Series) -> dict[str, int]:
    """Split flags on [;,|], drop PASS/empty tokens, count (legacy 2391-2398)."""
    flag_counts: dict[str, int] = {}
    for val in flags.fillna("").astype(str):
        toks = re.split(r"[;,|]+", val.strip())
        toks = [t.strip() for t in toks if t.strip() and t.strip().upper() != "PASS"]
        for t in toks:
            flag_counts[t] = flag_counts.get(t, 0) + 1
    return flag_counts


def build_qc_summary(tiered: pd.DataFrame, strict: pd.DataFrame) -> pd.DataFrame:
    """QC metric rows exactly as legacy (lines 2342-2434)."""
    nbd_col = find_col(tiered, ["has_nbd", "Has_NBD", "nbd_present", "NBD_present", "NBD", "hasNBD"])
    rep_col = find_col(tiered, ["has_repeat", "Has_repeat", "repeat_present", "Repeat_present", "HAS_REPEAT"])
    order_col = find_col(tiered, ["order_ok", "Order_OK", "orderOK", "NBD_repeat_order_ok", "ORDER_OK"])

    has_nbd = None
    has_rep = None
    order_ok = None

    if nbd_col:
        has_nbd = normalize_bool_series(tiered[nbd_col])
    if rep_col:
        has_rep = normalize_bool_series(tiered[rep_col])
    if order_col:
        order_ok = normalize_bool_series(tiered[order_col])

    if has_nbd is None:
        tmp = infer_has_nbd_from_domains(tiered)
        if tmp is not None:
            has_nbd = tmp.fillna(False)
    if has_rep is None:
        tmp = infer_has_repeat_from_domains(tiered)
        if tmp is not None:
            has_rep = tmp.fillna(False)

    # Order_OK cannot be reliably inferred from domains without positional info.
    if order_ok is None:
        order_ok = pd.Series([pd.NA] * len(tiered))

    flag_counts = count_flag_tokens(tiered["flags"]) if "flags" in tiered.columns else {}
    top_flags = sorted(flag_counts.items(), key=lambda x: x[1], reverse=True)[:20]

    qc_lines: list[tuple[str, object]] = []
    qc_lines.append(("n_total_candidates", len(tiered)))

    if "tier" in tiered.columns:
        for tier_name, count in tiered["tier"].value_counts(dropna=False).items():
            qc_lines.append((f"tier_count::{tier_name}", int(count)))

    if has_nbd is not None:
        qc_lines.append(("has_nbd_true", int(has_nbd.sum())))
    else:
        qc_lines.append(("has_nbd_true", "NA"))

    if has_rep is not None:
        qc_lines.append(("has_repeat_true", int(has_rep.sum())))
    else:
        qc_lines.append(("has_repeat_true", "NA"))

    # order_ok may have NA; count only True.
    try:
        qc_lines.append(("order_ok_true", int(pd.Series(order_ok).fillna(False).sum())))
    except Exception:  # noqa: BLE001 - legacy try/except
        qc_lines.append(("order_ok_true", "NA"))

    qc_lines.append(("n_strict_candidates", int(len(strict))))

    for k, v in top_flags:
        qc_lines.append((f"flag::{k}", int(v)))

    return pd.DataFrame(qc_lines, columns=["metric", "value"])


# -- BED writers (legacy lines 2444-2480) ------------------------------------

def bed_lines_original(out: pd.DataFrame) -> list[str]:
    sub = out.dropna(subset=["scaffold", "start", "end", "strand", "protein_id", "tier"])
    return [
        f"{r['scaffold']}\t{int(r['start']) - 1}\t{int(r['end'])}\t"
        f"{r['protein_id']}|{r['tier']}\t0\t{r['strand']}\n"
        for _, r in sub.iterrows()
    ]


def bed_lines_miniprot(out: pd.DataFrame) -> list[str]:
    sub = out.dropna(subset=["mp_scaffold", "mp_start", "mp_end", "mp_strand", "protein_id"])
    return [
        f"{r['mp_scaffold']}\t{int(r['mp_start']) - 1}\t{int(r['mp_end'])}\t"
        f"{r['protein_id']}|MINIPROT\t0\t{r['mp_strand']}\n"
        for _, r in sub.iterrows()
    ]


def bed_lines_exonerate(out: pd.DataFrame) -> list[str]:
    sub = out.dropna(subset=["exo_scaffold", "exo_start", "exo_end", "exo_strand", "protein_id"])
    return [
        f"{r['exo_scaffold']}\t{int(r['exo_start']) - 1}\t{int(r['exo_end'])}\t"
        f"{r['protein_id']}|EXONERATE\t0\t{r['exo_strand']}\n"
        for _, r in sub.iterrows()
    ]


def write_bed_pair(lines: list[str], stage_path: Path, final_path: Path) -> None:
    text = "".join(lines)
    stage_path.write_text(text)
    final_path.write_text(text)


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(ctx: RunContext) -> list[Path]:
    import json
    from funlr.stages.reporting_v3 import run_final
    from funlr.stages.rescue_v3 import manifest
    from funlr.evidence import export_evidence
    outputs = run_final(ctx)
    if not ctx.dry_run:
        outputs += export_evidence(
            ctx.paths.output_dir,
            [ctx.paths.stage_dir(7) / "evidence", ctx.paths.final_dir / "evidence"],
            ctx.config.scientific_settings(),
        )
        stage7 = ctx.paths.stage_dir(7)
        params = json.loads((stage7 / "metadata/run_info.json").read_text())
        manifest(stage7, "stage7_manifest.tsv", params)
    return outputs
