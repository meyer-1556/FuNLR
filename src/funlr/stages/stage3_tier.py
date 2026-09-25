"""Stage 3 entry point for the updated reference profile-aware tiering workflow.

The pure helpers below retain the earlier prototype rules for regression
comparisons. Production run() delegates to the updated :mod:`tiering` rules.
"""

from __future__ import annotations

import logging
import re
import time
from pathlib import Path

import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError

STAGE_NUMBER = 3
STAGE_NAME = "tier"
REQUIRED_TOOLS: list[str] = []

# Legacy required-column sets (exact, legacy lines 898-909).
REQUIRED_ARCH_COLUMNS = {
    "protein_id", "has_nbd", "has_repeat", "order_ok", "flag_fp_domains",
    "flag_rare_lrr", "domains_raw", "domains_grouped",
}
REQUIRED_MASTER_COLUMNS = {
    "protein_id", "scaffold", "start", "end", "strand", "protein_length",
}

# Tokens that should count as "NBD itself" (legacy lines 980-985).
NBD_TOKENS = {
    # common labels people use in summaries
    "NACHT", "NB-ARC", "NBD", "NACHT_NTPase", "STAND", "STAND_NTPase",
    # Pfam-ish names that sometimes appear
    "NACHT_N", "NB_ARC",
}


# ---------------------------------------------------------------------------
# Pure helpers (ported verbatim from the legacy heredoc; unit-testable)
# ---------------------------------------------------------------------------


def safe_int(x):
    try:
        if pd.isna(x):
            return None
        return int(x)
    except Exception:
        return None


def contig_end_flag(r, end_bp: int) -> bool:
    clen_ = safe_int(r.get("contig_len", None))
    s = safe_int(r.get("start", None))
    e = safe_int(r.get("end", None))
    if clen_ is None or s is None or e is None:
        return False
    return (s <= end_bp) or ((clen_ - e) <= end_bp)


def rescue_possible_flag(r, min_rescue_flank_bp: int) -> bool:
    clen_ = safe_int(r.get("contig_len", None))
    s = safe_int(r.get("start", None))
    e = safe_int(r.get("end", None))
    if clen_ is None or s is None or e is None:
        return False
    left_flank = s - 1
    right_flank = clen_ - e
    return max(left_flank, right_flank) >= min_rescue_flank_bp


def length_flags(r) -> list[str]:
    L = r.get("protein_length", None)
    flags = []
    try:
        if pd.isna(L):
            return flags
        L = int(L)
    except Exception:
        return flags
    if L < 500:
        flags.append("SHORT_PROTEIN")
    if L > 3000:
        flags.append("VERY_LONG_PROTEIN")
    return flags


# ---------- Domain parsing for Tier 3 logic ----------
# We interpret domains_grouped / domains_raw as a delimiter-separated set.
def _split_domains(s) -> list[str]:
    s = str(s or "")
    if not s:
        return []
    # common separators: ; , | whitespace
    parts = re.split(r"[;,\|\s]+", s)
    return [p.strip() for p in parts if p and p.strip()]


def has_other_domains_besides_nbd(r) -> bool:
    doms = _split_domains(r.get("domains_grouped", ""))
    if not doms:
        doms = _split_domains(r.get("domains_raw", ""))
    # If we still have nothing, be conservative: treat as unknown -> no "other domains"
    if not doms:
        return False
    # If ANY domain token is not an NBD token, we consider "other domains present"
    for d in doms:
        if d not in NBD_TOKENS:
            return True
    return False


# ---------- Annotation-based "likely non-NLR" flagging (soft) ----------
def _compile_or_none(pat, warn=None):
    pat = (pat or "").strip()
    if not pat:
        return None
    try:
        return re.compile(pat, re.IGNORECASE)
    except re.error as e:
        msg = f"WARNING: bad regex pattern ignored: {pat!r} ({e})"
        if warn is not None:
            warn(msg)
        else:
            logging.getLogger("funlr.stage3").warning(msg)
        return None


def non_nlr_annot_flag(r, re_kw, re_go, re_dom) -> bool:
    txt = " ".join([str(r.get("Preferred_name", "")), str(r.get("Description", ""))])
    gos = str(r.get("GOs", ""))
    dom = str(r.get("domains_raw", ""))

    hit = False
    if re_kw and re_kw.search(txt):
        hit = True
    if re_go and re_go.search(gos):
        hit = True
    if re_dom and re_dom.search(dom):
        hit = True
    return bool(hit)


# ---------- Repeat proximity heuristic (soft; verbatim port of legacy loop) ----------
def compute_repeat_nearby(
    df: pd.DataFrame, repeat_prox_bp: int, repeat_prox_neighbors: int
) -> pd.Series:
    """Scaffold-grouped +-window / +-neighbor-count repeat-proximity scan.

    Only candidates with has_nbd and not has_repeat are considered; any
    repeat-bearing neighbor within ``repeat_prox_neighbors`` positions whose
    span overlaps the +-repeat_prox_bp window marks the candidate.
    """
    repeat_nearby = pd.Series(False, index=df.index)

    for scaf, sub in df.dropna(subset=["scaffold", "_start_i", "_end_i"]).groupby(
        "scaffold", sort=False
    ):
        sub = sub.sort_values("_start_i").copy()
        idxs = sub.index.tolist()

        is_repeat = sub["has_repeat"].fillna(False).astype(bool)
        starts = sub["_start_i"].astype(int).tolist()
        ends = sub["_end_i"].astype(int).tolist()

        for j, idx in enumerate(idxs):
            has_nbd = bool(sub.loc[idx, "has_nbd"])
            has_rep = bool(sub.loc[idx, "has_repeat"])
            if (not has_nbd) or has_rep:
                continue

            s = starts[j]
            e = ends[j]
            win_s = s - repeat_prox_bp
            win_e = e + repeat_prox_bp

            j0 = max(0, j - repeat_prox_neighbors)
            j1 = min(len(idxs) - 1, j + repeat_prox_neighbors)

            nearby = False
            for k in range(j0, j1 + 1):
                if k == j:
                    continue
                if not bool(is_repeat.iloc[k]):
                    continue
                if ends[k] < win_s or starts[k] > win_e:
                    continue
                nearby = True
                break

            repeat_nearby.at[idx] = nearby

    return repeat_nearby.fillna(False).astype(bool)


# ---------- Flags ----------
def illumina_arch_flags(r, end_bp: int, min_rescue_flank_bp: int) -> list[str]:
    flags = []
    has_nbd = bool(r.get("has_nbd", False))
    has_rep = bool(r.get("has_repeat", False))

    # NBD_ONLY is now *strict*: NBD present, no repeat, and no other domains besides NBD
    if has_nbd and (not has_rep) and (not has_other_domains_besides_nbd(r)):
        flags.append("NBD_ONLY")

    if not bool(r.get("order_ok", True)):
        flags.append("DOMAIN_ORDER_WEIRD")

    if bool(r.get("flag_fp_domains", False)):
        flags.append("FP_DOMAIN_PRESENT")

    if bool(r.get("flag_rare_lrr", False)):
        flags.append("FLAG_RARE_LRR")

    if bool(r.get("non_nlr_annot", False)):
        flags.append("NON_NLR_ANNOT")

    if contig_end_flag(r, end_bp):
        flags.append("CONTIG_END")
        if not rescue_possible_flag(r, min_rescue_flank_bp):
            flags.append("RESCUE_NO_FLANK")
        else:
            flags.append("RESCUE_POSSIBLE")

    if bool(r.get("repeat_nearby", False)):
        flags.append("REPEAT_NEARBY")

    flags.extend(length_flags(r))
    return flags


# ---------- Tiering (new scheme + Tier 3) ----------
def assign_tier(r) -> str:
    has_nbd = bool(r.get("has_nbd", False))
    has_rep = bool(r.get("has_repeat", False))
    flags = r.get("flags", "") or ""
    non_nlr = "NON_NLR_ANNOT" in flags

    # Tier 4 split: NO NBD cases
    if not has_nbd:
        if has_rep:
            return "TIER_4A_REPEAT_ONLY_NO_NBD"
        return "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL"

    # Tier 4 split: NBD present but strong housekeeping signature
    if non_nlr:
        return "TIER_4B_LIKELY_STAND_HOUSEKEEPING"

    # NBD + repeat (canonical-ish)
    if has_rep and ("FP_DOMAIN_PRESENT" not in flags):
        if ("CONTIG_END" in flags) or ("SHORT_PROTEIN" in flags):
            return "TIER_1B_NEEDS_REVIEW"
        return "TIER_1A_HIGH_CONFIDENCE"

    # NBD present, no repeat
    # Rescue-worthy first:
    if ("CONTIG_END" in flags) or ("REPEAT_NEARBY" in flags):
        if "RESCUE_NO_FLANK" in flags:
            # still rescue-ish, but weaker
            return "TIER_2B_RESCUE_CANDIDATE"
        return (
            "TIER_2A_HIGH_PRIORITY_RESCUE"
            if ("CONTIG_END" in flags)
            else "TIER_2B_RESCUE_CANDIDATE"
        )

    # True "NBD only" (strict) but not near contig ends / repeats -> lower priority rescue candidate
    if "NBD_ONLY" in flags:
        return "TIER_2B_RESCUE_CANDIDATE"

    # Otherwise: architectural variants (NBD + other domains, but no repeat)
    return "TIER_3_ARCHITECTURAL_VARIANT"


# ---------- Rescue priority scoring ----------
def rescue_score(r) -> int:
    score = 0
    t = r.get("tier", "")
    flags = r.get("flags", "") or ""

    if t == "TIER_2A_HIGH_PRIORITY_RESCUE":
        score += 90
    elif t == "TIER_2B_RESCUE_CANDIDATE":
        score += 70
    elif t == "TIER_1B_NEEDS_REVIEW":
        score += 60

    if "CONTIG_END" in flags:
        score += 15
    if "RESCUE_POSSIBLE" in flags:
        score += 10
    if "RESCUE_NO_FLANK" in flags:
        score -= 30

    if "REPEAT_NEARBY" in flags:
        score += 10

    if "SHORT_PROTEIN" in flags:
        score += 10
    if "FP_DOMAIN_PRESENT" in flags:
        score -= 25
    if "NON_NLR_ANNOT" in flags:
        score -= 35
    if "VERY_LONG_PROTEIN" in flags:
        score -= 10

    return max(0, min(100, score))


# ---------- Full table build (legacy merge + flag/tier/score pipeline) ----------


def build_tiered_table(
    master: pd.DataFrame,
    arch: pd.DataFrame,
    clen: pd.DataFrame,
    *,
    contig_end_bp: int,
    min_rescue_flank_bp: int,
    repeat_prox_bp: int,
    repeat_prox_neighbors: int,
    non_nlr_keywords_regex: str = "",
    non_nlr_domain_regex: str = "",
    non_nlr_go_regex: str = "",
    log: logging.Logger | None = None,
) -> pd.DataFrame:
    """Merge arch+master+contig lengths and compute flags/tier/rescue_priority.

    Column validation, dedup, NA-fraction warning and every flag/tier/score
    step match the legacy heredoc exactly (legacy lines 897-1185).
    """
    warn = log.warning if log is not None else logging.getLogger("funlr.stage3").warning

    # Ensure key columns exist
    missing_arch = REQUIRED_ARCH_COLUMNS - set(arch.columns)
    if missing_arch:
        raise StageError(
            f"ERROR: architecture_summary.tsv missing columns: {sorted(missing_arch)}"
        )

    missing_master = REQUIRED_MASTER_COLUMNS - set(master.columns)
    if missing_master:
        raise StageError(
            f"ERROR: master_table.tsv missing columns: {sorted(missing_master)}"
        )

    # Defensive: protein_id should be unique in arch
    if arch["protein_id"].duplicated().any():
        dups = int(arch["protein_id"].duplicated().sum())
        warn(
            "WARNING: architecture_summary.tsv has %d duplicate protein_id rows; "
            "keeping first occurrence.",
            dups,
        )
        arch = arch.drop_duplicates(subset=["protein_id"], keep="first")

    # Merge data
    try:
        df = arch.merge(master, on="protein_id", how="left", validate="one_to_one").merge(
            clen, on="scaffold", how="left"
        )
    except pd.errors.MergeError as exc:
        raise StageError(f"ERROR: merge failed (one_to_one): {exc}") from exc

    # A valid proteome can yield no candidates. Preserve the ordinary output
    # schema without relying on pandas' ambiguous empty DataFrame.apply shape.
    if df.empty:
        for name, dtype in (
            ("non_nlr_annot", bool), ("_start_i", float), ("_end_i", float),
            ("repeat_nearby", bool), ("flags_list", object), ("flags", object),
            ("tier", object), ("rescue_priority", int),
        ):
            df[name] = pd.Series(index=df.index, dtype=dtype)
        return df

    na_frac = df["contig_len"].isna().mean()
    if na_frac > 0.25:
        warn(
            "WARNING: contig_len is NA for %.1f%% of rows. Check scaffold names "
            "between master_table and contig_lengths.",
            na_frac * 100,
        )

    re_kw = _compile_or_none(non_nlr_keywords_regex, warn)
    re_dom = _compile_or_none(non_nlr_domain_regex, warn)
    re_go = _compile_or_none(non_nlr_go_regex, warn)

    df["non_nlr_annot"] = df.apply(
        lambda r: non_nlr_annot_flag(r, re_kw, re_go, re_dom), axis=1
    )

    # Repeat proximity heuristic
    df["_start_i"] = pd.to_numeric(df["start"], errors="coerce")
    df["_end_i"] = pd.to_numeric(df["end"], errors="coerce")
    df["repeat_nearby"] = compute_repeat_nearby(df, repeat_prox_bp, repeat_prox_neighbors)

    df["flags_list"] = df.apply(
        lambda r: illumina_arch_flags(r, contig_end_bp, min_rescue_flank_bp), axis=1
    )
    df["flags"] = df["flags_list"].apply(lambda xs: ";".join(xs) if xs else "PASS")
    df["tier"] = df.apply(assign_tier, axis=1)
    df["rescue_priority"] = df.apply(rescue_score, axis=1)
    return df


# ---------- Outputs ----------


def write_bed(path, frame: pd.DataFrame) -> None:
    """Legacy BED writer: 0-based start, name=protein_id|tier, score=priority*10."""
    with open(path, "w") as o:
        for _, r in frame.iterrows():
            o.write(
                f"{r['scaffold']}\t{int(r['start']) - 1}\t{int(r['end'])}"
                f"\t{r['protein_id']}|{r['tier']}\t{int(r['rescue_priority'] * 10)}\t{r['strand']}\n"
            )


def write_outputs(
    df: pd.DataFrame, tables_dir: Path, beds_dir: Path, summ_dir: Path
) -> list[Path]:
    """Write every legacy stage-3 output file (legacy lines 1187-1239)."""
    tables_dir = Path(tables_dir)
    beds_dir = Path(beds_dir)
    summ_dir = Path(summ_dir)
    outputs: list[Path] = []

    tiered_path = tables_dir / "tiered_candidates.tsv"
    df.drop(columns=["flags_list"], errors="ignore").to_csv(
        tiered_path, sep="\t", index=False
    )
    outputs.append(tiered_path)

    for tier, sub in df.groupby("tier", dropna=False):
        out = tables_dir / f"{tier}.tsv"
        sub.drop(columns=["flags_list"], errors="ignore").to_csv(out, sep="\t", index=False)
        outputs.append(out)

    tier_series = df["tier"].fillna("")
    resc_mask = tier_series.str.contains("RESCUE", regex=False) | (
        tier_series == "TIER_1B_NEEDS_REVIEW"
    )
    resc = df[resc_mask].sort_values("rescue_priority", ascending=False)
    resc_path = tables_dir / "rescue_priority.tsv"
    resc.drop(columns=["flags_list"], errors="ignore").to_csv(
        resc_path, sep="\t", index=False
    )
    outputs.append(resc_path)

    bed = beds_dir / "nlr_candidates.bed"
    bed_df = df.dropna(subset=["scaffold", "start", "end", "strand"]).copy()
    write_bed(bed, bed_df)
    outputs.append(bed)

    by_tier_dir = beds_dir / "by_tier"
    by_tier_dir.mkdir(parents=True, exist_ok=True)
    for tier, sub in bed_df.groupby("tier", dropna=False):
        safe_tier = str(tier).replace("/", "_")
        out = by_tier_dir / f"{safe_tier}.bed"
        write_bed(out, sub)
        outputs.append(out)

    summary = df["tier"].value_counts().to_frame("count")
    tier_summary = summ_dir / "tier_summary.tsv"
    summary.to_csv(tier_summary, sep="\t")
    outputs.append(tier_summary)

    flag_counts = (
        df["flags"]
        .fillna("")
        .str.split(";")
        .explode()
        .replace({"PASS": pd.NA, "": pd.NA})
        .dropna()
        .value_counts()
    )
    flag_counts_path = summ_dir / "flag_counts.tsv"
    flag_counts.to_frame("count").to_csv(flag_counts_path, sep="\t")
    outputs.append(flag_counts_path)

    return outputs


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------


def run(ctx: RunContext) -> list[Path]:
    """Tier the updated reference architecture table using explicit profile settings."""
    from .tiering import TierPolicy, build_profile_table, write_profile_outputs

    log = ctx.stage_logger(STAGE_NUMBER, STAGE_NAME)
    if ctx.dry_run:
        log.info("PLANNED: profile-aware tiering, rescue selection, and BED exports")
        return []
    settings = ctx.config.scientific_settings()
    policy = TierPolicy.from_settings(settings)
    stage0 = ctx.paths.results_dir / "stage0"
    stage2 = ctx.paths.results_dir / "stage2"
    stage3 = ctx.paths.stage_dir(STAGE_NUMBER)
    inputs = {
        "master": stage0 / "master_table.tsv",
        "contig": stage0 / "contig_lengths.tsv",
        "architecture": stage2 / "architecture_summary.tsv",
        "input_manifest": stage0 / "input_manifest.tsv",
        "software_manifest": stage0 / "software_manifest.tsv",
        "stage2_parse_filter_report": stage2 / "parse_filter_report.tsv",
        "stage2_run_info": stage2 / "run_info.txt",
    }
    for path in inputs.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise StageError(f"Required input missing/empty: {path}")
    try:
        master = pd.read_csv(inputs["master"], sep="\t")
        arch = pd.read_csv(inputs["architecture"], sep="\t")
        contigs = pd.read_csv(inputs["contig"], sep="\t", header=None, names=["scaffold", "contig_len"])
    except (pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise StageError(f"Malformed Stage 3 input table: {exc}") from exc
    frame = build_profile_table(master, arch, contigs, policy, settings=settings, log=log)
    outputs = write_profile_outputs(frame, stage3, policy)
    metadata = stage3 / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    run_info = metadata / "run_info.txt"
    info = {
        "stage": STAGE_NUMBER, "date": time.ctime(),
        "genome": ctx.config.get("sample", "species_id") or ctx.sample_id,
        "assembly_version": ctx.config.get("sample", "assembly_version") or "",
        "nlr_profile": policy.profile, **inputs,
        **{name: getattr(policy, name) for name in policy.__dataclass_fields__ if name != "profile"},
        "non_nlr_keywords_regex": settings.get("NON_NLR_KEYWORDS_REGEX", ""),
        "non_nlr_domain_regex": settings.get("NON_NLR_DOMAIN_REGEX", ""),
        "non_nlr_go_regex": settings.get("NON_NLR_GO_REGEX", ""),
        "confidence_policy": "Use supplied confidence; infer legacy evidence only if confidence column is absent",
        "stage2_header_columns": len(arch.columns),
    }
    run_info.write_text("".join(f"{key}: {value}\n" for key, value in info.items()))
    outputs.append(run_info)
    log.info("Tiered %d candidates with profile %s", len(frame), policy.profile)
    log.info("Tier counts: %s", frame["tier"].value_counts().to_dict())
    return outputs
