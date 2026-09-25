"""Profile-aware tiering from the completed April 2026 analysis policy.

The scientific decision order and score terms follow 03_tier_and_flags.sh.
An explicitly supplied LOW confidence is authoritative; legacy inference is
used only when the confidence column is absent. Empty candidate tables retain
their schema. The old stage3 helpers remain available for baseline comparisons.
"""
from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
import re

import pandas as pd

from funlr.core.errors import StageError


REQUIRED_ARCH = {
    "protein_id", "has_nbd", "has_sensor", "has_ssfr_sensor",
    "has_non_ssfr_sensor", "domains_raw", "domains_grouped",
}
REQUIRED_MASTER = {"protein_id", "scaffold", "start", "end", "strand", "protein_length"}
OPTIONAL_BOOL = (
    "has_lrr", "order_invalid", "flag_fp_domains", "flag_fp_overlaps_nbd",
    "integrated_domain_cterm", "asm_present", "has_any_repeat_region",
    "has_effector", "has_effector_nterm", "effector_is_nterm_of_nbd",
)
NBD_TOKENS = {"NACHT", "NB-ARC", "STAND_LIKE", "NBD_LIKE"}
TEMP_COLUMNS = ["flags_list", "_start_i", "_end_i"]


def safe_int(value):
    try:
        return None if pd.isna(value) else int(value)
    except (TypeError, ValueError, OverflowError):
        return None


def as_bool(value) -> bool:
    """Treat missing calls as false, without treating text 'False' as true."""
    if value is None or pd.isna(value):
        return False
    if isinstance(value, str):
        if value.strip().lower() in {"true", "1"}:
            return True
        if value.strip().lower() in {"false", "0", "", "na", "nan"}:
            return False
        raise StageError(f"Invalid boolean architecture value: {value!r}")
    return bool(value)


def other_domains(row) -> bool:
    # Keep the source's delimiter and token rules, including its treatment of
    # an explicit textual 'nan' domain token; missing values are not domains.
    for column in ("domains_grouped", "domains_raw"):
        value = row.get(column, "")
        text = "" if value is None or pd.isna(value) else str(value)
        domains = [p for p in re.split(r"[;,|\s]+", text) if p]
        if domains:
            return any(domain not in NBD_TOKENS for domain in domains)
    return False


@dataclass(frozen=True)
class TierPolicy:
    profile: str = "ILLUMINA"
    hard_end_bp: int = 10000
    subtelo_bp: int = 0
    min_rescue_flank_bp: int = 5000
    repeat_prox_bp: int = 10000
    repeat_prox_neighbors: int = 10
    very_short_protein_aa: int = 333
    short_protein_aa: int = 500
    skip_rescue_for_tier1a: bool = True
    include_tier3_in_rescue: bool = False

    @classmethod
    def from_settings(cls, settings):
        fields = {name: settings[name.upper()] for name in cls.__dataclass_fields__
                  if name != "profile" and name.upper() in settings}
        fields["profile"] = str(settings.get("NLR_PROFILE", "ILLUMINA")).upper()
        return cls(**fields)

    def distance_to_end(self, row):
        length, start, end = [safe_int(row.get(k)) for k in ("contig_len", "start", "end")]
        if None in (length, start, end):
            return None
        return min(start - 1, length - end)

    def rescue_possible(self, row) -> bool:
        length, start, end = [safe_int(row.get(k)) for k in ("contig_len", "start", "end")]
        return None not in (length, start, end) and max(start - 1, length - end) >= self.min_rescue_flank_bp

    def length_flags(self, row) -> list[str]:
        length = safe_int(row.get("protein_length"))
        if length is None:
            return []
        flags = []
        if length < self.very_short_protein_aa:
            flags.append("VERY_SHORT_PROTEIN")
        elif length < self.short_protein_aa:
            flags.append("SHORT_PROTEIN")
        if length > 3000:
            flags.append("VERY_LONG_PROTEIN")
        return flags

    def build_flags(self, row) -> list[str]:
        flags = []
        sensor = as_bool(row.get("has_sensor", False))
        nbd_ok = as_bool(row.get("nbd_ok", False))
        if as_bool(row.get("has_effector", False)):
            flags.append("HAS_EFFECTOR_ANY")
        if as_bool(row.get("has_effector_nterm", row.get("effector_is_nterm_of_nbd", False))):
            flags.append("HAS_EFFECTOR_NTERM")
        if nbd_ok and not sensor and not other_domains(row):
            flags.append("NBD_ONLY_UNKNOWN_SENSOR")
        if as_bool(row.get("order_invalid", False)):
            flags.append("ORIENTATION_INVALID")
        else:
            orientation = row.get("sensor_is_cterm_of_nbd")
            if nbd_ok and sensor and not pd.isna(orientation) and orientation is not None and not as_bool(orientation):
                flags.append("DOMAIN_ORDER_WEIRD")
        for column, flag in (
            ("flag_fp_domains", "FP_DOMAIN_PRESENT"),
            ("flag_fp_overlaps_nbd", "FP_DOMAIN_OVERLAPS_NBD"),
            ("has_lrr", "EXCLUDE_LRR"), ("non_nlr_annot", "NON_NLR_ANNOT"),
            ("integrated_domain_cterm", "INTEGRATED_DOMAIN_CTERM"),
            ("asm_present", "ASM_PRESENT"),
        ):
            if as_bool(row.get(column, False)):
                flags.append(flag)
        distance = self.distance_to_end(row)
        if distance is not None and distance <= self.hard_end_bp:
            flags.extend(["HARD_END", "RESCUE_POSSIBLE" if self.rescue_possible(row) else "RESCUE_NO_FLANK"])
        if self.subtelo_bp > 0 and distance is not None and distance <= self.subtelo_bp:
            flags.append("SUBTELO")
        if as_bool(row.get("repeat_nearby", False)):
            flags.append("REPEAT_NEARBY")
        flags.extend(self.length_flags(row))
        if "NBD_ONLY_UNKNOWN_SENSOR" in flags and any(
            flag in flags for flag in ("HARD_END", "SHORT_PROTEIN", "VERY_SHORT_PROTEIN")
        ):
            flags.append("NBD_ONLY_TRUNCATED")
        return flags

    def assign_tier(self, row) -> str:
        raw_nbd = as_bool(row.get("has_nbd_any", row.get("has_nbd", False)))
        confidence = str(row.get("nbd_confidence", "LOW")).strip().upper()
        nbd_ok = as_bool(row.get("nbd_ok", confidence in {"HIGH", "MEDIUM"}))
        sensor = as_bool(row.get("has_sensor", False))
        flags = row.get("flags", "") or ""
        hard_end, nearby = "HARD_END" in flags, "REPEAT_NEARBY" in flags
        no_flank = "RESCUE_NO_FLANK" in flags
        # Preserve the reference substring rule (VERY_SHORT_PROTEIN also contains
        # SHORT_PROTEIN); changing that rule is a separate scientific decision.
        very_short, short = "VERY_SHORT_PROTEIN" in flags, "SHORT_PROTEIN" in flags
        repeat_region = as_bool(row.get("has_any_repeat_region", False))
        asm = as_bool(row.get("asm_present", False))
        if not raw_nbd:
            return "TIER_4A_REPEAT_ONLY_NO_NBD" if sensor else "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL"
        if as_bool(row.get("has_lrr", False)):
            return "TIER_4D_EXCLUDED_LRR"
        if as_bool(row.get("order_invalid", False)) and not (not sensor and repeat_region):
            return "TIER_4E_ORIENTATION_INVALID"
        if as_bool(row.get("non_nlr_annot", False)) and confidence != "HIGH":
            return "TIER_4B_LIKELY_STAND_HOUSEKEEPING"
        if very_short and not (hard_end or nearby or asm):
            return "TIER_4F_VERY_SHORT_FRAGMENT"
        if nbd_ok:
            if sensor and as_bool(row.get("has_ssfr_sensor", False)) and not as_bool(row.get("flag_fp_overlaps_nbd", False)):
                return "TIER_1B_NEEDS_REVIEW" if hard_end or short else "TIER_1A_HIGH_CONFIDENCE"
            if sensor and as_bool(row.get("has_non_ssfr_sensor", False)) and not as_bool(row.get("flag_fp_overlaps_nbd", False)):
                return "TIER_1C_CANONICAL_NON_SSFR_SENSOR"
        if nbd_ok and not sensor and repeat_region:
            if self.profile == "ILLUMINA" and hard_end and not no_flank:
                return "TIER_2A_HIGH_PRIORITY_RESCUE"
            return "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE"
        rescue_context = hard_end or nearby or asm
        if not nbd_ok:
            if rescue_context:
                return "TIER_2C_LOW_PRIORITY_FRAGMENT" if no_flank else "TIER_2B_RESCUE_CANDIDATE"
            return "TIER_3B_ARCHITECTURAL_VARIANT"
        if rescue_context:
            if no_flank:
                return "TIER_2B_RESCUE_CANDIDATE"
            if self.profile == "ILLUMINA" and hard_end:
                return "TIER_2A_HIGH_PRIORITY_RESCUE"
            return "TIER_2B_RESCUE_CANDIDATE"
        if "NBD_ONLY_UNKNOWN_SENSOR" in flags:
            return "TIER_4G_SUSPECTED_PSEUDOGENE" if self.profile == "HIFI" else "TIER_2B_RESCUE_CANDIDATE"
        if short:
            return "TIER_2C_LOW_PRIORITY_FRAGMENT"
        if as_bool(row.get("integrated_domain_cterm", False)):
            return "TIER_3A_INTEGRATED_DECOY"
        return "TIER_3B_ARCHITECTURAL_VARIANT"

    def rescue_score(self, row) -> int:
        tier = row.get("tier", "")
        flags = row.get("flags", "") or ""
        confidence = str(row.get("nbd_confidence", "LOW")).strip().upper()
        bases = {"TIER_2A_HIGH_PRIORITY_RESCUE": 90, "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE": 75,
                 "TIER_2B_RESCUE_CANDIDATE": 70, "TIER_1B_NEEDS_REVIEW": 60,
                 "TIER_1C_CANONICAL_NON_SSFR_SENSOR": 50, "TIER_2C_LOW_PRIORITY_FRAGMENT": 40,
                 "TIER_3A_INTEGRATED_DECOY": 30, "TIER_3B_ARCHITECTURAL_VARIANT": 20}
        score = bases.get(tier, 0)
        for flag, amount in (
            ("HAS_EFFECTOR_NTERM", 20), ("HARD_END", 15 if self.profile == "ILLUMINA" else 5),
            ("RESCUE_POSSIBLE", 10), ("RESCUE_NO_FLANK", -30), ("REPEAT_NEARBY", 10), ("SUBTELO", 5),
            ("ASM_PRESENT", 25), ("INTEGRATED_DOMAIN_CTERM", -10), ("SHORT_PROTEIN", 10),
            ("VERY_SHORT_PROTEIN", -40), ("VERY_LONG_PROTEIN", -10),
            ("FP_DOMAIN_PRESENT", -25), ("FP_DOMAIN_OVERLAPS_NBD", -50),
            ("NON_NLR_ANNOT", -15 if confidence == "HIGH" else -45),
        ):
            if flag in flags:
                score += amount
        if confidence == "HIGH":
            score += 10
            if "NBD_ONLY_UNKNOWN_SENSOR" in flags or "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE" in tier:
                score += 5
        elif confidence == "MEDIUM":
            score += 5
        else:
            score -= 10
        if "EXCLUDE_LRR" in flags or "ORIENTATION_INVALID" in flags:
            score = 0
        return max(0, min(100, score))


def repeat_nearby(frame, policy):
    result = pd.Series(False, index=frame.index)
    for _, group in frame.dropna(subset=["scaffold", "_start_i", "_end_i"]).groupby("scaffold", sort=False):
        group = group.sort_values("_start_i")
        indices = group.index.tolist()
        starts, ends = group["_start_i"].astype(int).tolist(), group["_end_i"].astype(int).tolist()
        repeats = group["has_ssfr_sensor"].map(as_bool).tolist()
        for position, index in enumerate(indices):
            if not as_bool(group.loc[index, "nbd_ok"]) or as_bool(group.loc[index, "has_sensor"]):
                continue
            first = max(0, position - policy.repeat_prox_neighbors)
            last = min(len(indices) - 1, position + policy.repeat_prox_neighbors)
            result.at[index] = any(
                k != position and repeats[k] and ends[k] >= starts[position] - policy.repeat_prox_bp
                and starts[k] <= ends[position] + policy.repeat_prox_bp for k in range(first, last + 1)
            )
    return result


def build_profile_table(master, arch, contigs, policy, *, settings=None, log=None):
    """Return tiered rows without modifying any caller-owned input table."""
    log = log or logging.getLogger("funlr.stage3")
    settings = settings or {}
    arch = arch.copy()
    malformed = [c for c in arch if "sensor_is_cterm_of_nbdorder_invalid" in str(c)]
    if malformed:
        raise StageError(f"Malformed Stage 2 header (missing tab): {malformed[0]}")
    for frame, required, label in ((arch, REQUIRED_ARCH, "architecture_summary.tsv"), (master, REQUIRED_MASTER, "master_table.tsv")):
        missing = required - set(frame.columns)
        if missing:
            raise StageError(f"{label} missing required columns: {sorted(missing)}")
    if {"scaffold", "contig_len"} - set(contigs.columns):
        raise StageError("contig_lengths.tsv requires scaffold and contig_len")
    supplied_confidence = "nbd_confidence" in arch
    duplicate_ids = arch.loc[arch["protein_id"].duplicated(keep=False), ["protein_id"]].drop_duplicates()
    if not duplicate_ids.empty:
        log.warning("Duplicate architecture IDs: keeping first occurrence for %d IDs", len(duplicate_ids))
        arch = arch.drop_duplicates("protein_id", keep="first")
    for name in OPTIONAL_BOOL:
        if name not in arch:
            arch[name] = False
    if "sensor_is_cterm_of_nbd" not in arch:
        arch["sensor_is_cterm_of_nbd"] = pd.NA
    if not supplied_confidence:
        arch["nbd_confidence"] = "LOW"
    for name in ("has_nbd", "has_sensor", "has_ssfr_sensor", "has_non_ssfr_sensor", *OPTIONAL_BOOL, "has_nbd_stage1", "has_nbd_pfam"):
        if name in arch:
            arch[name] = arch[name].map(as_bool).astype(bool)
    try:
        frame = arch.merge(master, on="protein_id", how="left", validate="one_to_one").merge(contigs, on="scaffold", how="left")
    except pd.errors.MergeError as exc:
        raise StageError(f"Stage 3 one-to-one merge failed: {exc}") from exc
    frame["nbd_confidence"] = frame["nbd_confidence"].fillna("LOW").astype(str).str.strip().str.upper().replace({"NA": "LOW", "NONE": "LOW", "": "LOW"})
    if not supplied_confidence:
        log.warning("Stage 2 omitted nbd_confidence; inferring only from legacy NBD evidence columns")
        if "has_nbd_stage1" in frame:
            frame.loc[frame["has_nbd_stage1"], "nbd_confidence"] = "HIGH"
        if "has_nbd_pfam" in frame:
            frame.loc[frame["has_nbd_pfam"] & frame["nbd_confidence"].eq("LOW"), "nbd_confidence"] = "MEDIUM"
    frame["nbd_ok"] = frame["nbd_confidence"].isin(["HIGH", "MEDIUM"])
    frame["has_nbd_any"] = frame["has_nbd"].map(as_bool).astype(bool)
    expressions = []
    for key in ("NON_NLR_KEYWORDS_REGEX", "NON_NLR_GO_REGEX", "NON_NLR_DOMAIN_REGEX"):
        text = settings.get(key, "") or ""
        try:
            expressions.append(re.compile(text, re.IGNORECASE) if text else None)
        except re.error as exc:
            raise StageError(f"Invalid {key}: {exc}") from exc
    def non_nlr(row):
        texts = [" ".join(str(row.get(k, "")) for k in ("Preferred_name", "Description")), str(row.get("GOs", "")), str(row.get("domains_raw", ""))]
        return any(regex is not None and regex.search(text) is not None for regex, text in zip(expressions, texts))
    frame["non_nlr_annot"] = pd.Series([non_nlr(row) for _, row in frame.iterrows()], index=frame.index, dtype=bool)
    frame["_start_i"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["_end_i"] = pd.to_numeric(frame["end"], errors="coerce")
    frame["repeat_nearby"] = repeat_nearby(frame, policy)
    frame["flags_list"] = pd.Series([policy.build_flags(row) for _, row in frame.iterrows()], index=frame.index, dtype=object)
    frame["flags"] = frame["flags_list"].map(lambda values: ";".join(values) if values else "PASS").astype(object)
    frame["tier"] = pd.Series([policy.assign_tier(row) for _, row in frame.iterrows()], index=frame.index, dtype=object)
    frame["rescue_priority"] = pd.Series([policy.rescue_score(row) for _, row in frame.iterrows()], index=frame.index, dtype=int)
    if not frame.empty and frame["contig_len"].isna().mean() > 0.25:
        log.warning("More than 25%% of candidates have no contig length; check scaffold names")
    frame.attrs["duplicate_arch_ids"] = duplicate_ids["protein_id"].tolist()
    return frame


def write_profile_outputs(frame, stage3, policy):
    """Write the reference layout and clean schemas, including a zero-candidate run."""
    stage3 = Path(stage3)
    tables, beds, summaries, qc = [stage3 / name for name in ("tables", "beds", "summaries", "qc")]
    for directory in (tables, beds, summaries, qc, beds / "by_tier"):
        directory.mkdir(parents=True, exist_ok=True)
    outputs = []
    def tsv(data, path, *, index=False):
        data.to_csv(path, sep="\t", index=index)
        outputs.append(path)
    duplicate_arch = frame.attrs.get("duplicate_arch_ids")
    if duplicate_arch:
        tsv(pd.DataFrame({"protein_id": duplicate_arch}), qc / "stage2_arch_duplicate_protein_ids.tsv")
    # Preserve the defined universe counts before coordinate-join deduplication.
    metrics = {
        "total_candidates": len(frame), "has_nbd_any": frame["has_nbd_any"].sum(), "has_nbd_ok": frame["nbd_ok"].sum(),
        "has_sensor": frame["has_sensor"].sum(), "has_ssfr_sensor": frame["has_ssfr_sensor"].sum(),
        "has_non_ssfr_sensor": frame["has_non_ssfr_sensor"].sum(),
        "has_nbd_any_and_sensor": (frame["has_nbd_any"] & frame["has_sensor"]).sum(),
        "has_nbd_ok_and_sensor": (frame["nbd_ok"] & frame["has_sensor"]).sum(),
        "has_nbd_ok_and_ssfr_sensor": (frame["nbd_ok"] & frame["has_ssfr_sensor"]).sum(),
        "non_nlr_annot": frame["non_nlr_annot"].sum(),
        "has_nbd_ok_and_non_nlr_annot": (frame["nbd_ok"] & frame["non_nlr_annot"]).sum(), "has_lrr": frame["has_lrr"].sum(),
    }
    tsv(pd.DataFrame([{"metric": key, "count": int(count), "percent": count / len(frame) * 100 if len(frame) else 0.0} for key, count in metrics.items()]), summaries / "universe_summary.tsv")
    if frame["protein_id"].duplicated().any():
        frame = frame.copy()
        frame["_tier_rank"] = frame["tier"].str.extract(r"^TIER_([1-4])", expand=False).fillna(99).astype(int)
        frame["_conf_rank"] = frame["nbd_confidence"].map({"HIGH": 3, "MEDIUM": 2, "LOW": 1}).fillna(0).astype(int)
        order = ["protein_id", "has_nbd_any", "rescue_priority", "_tier_rank", "_conf_rank"]
        frame = frame.sort_values(order, ascending=[True, False, False, True, False])
        duplicates = frame[frame["protein_id"].duplicated(keep=False)].copy()
        tsv(duplicates, qc / "duplicate_protein_ids.full.tsv")
        columns = [c for c in ("tier", "rescue_priority", "nbd_confidence", "flags", "domains_grouped", "scaffold", "start", "end") if c in duplicates]
        duplicates["_sig"] = duplicates[columns].astype(str).agg("|".join, axis=1)
        tsv(duplicates.groupby("protein_id", as_index=False).agg(duplicate_rows=("protein_id", "size"), unique_signatures=("_sig", "nunique")), qc / "duplicate_protein_ids.summary.tsv")
        frame = frame.drop_duplicates("protein_id", keep="first").drop(columns=["_tier_rank", "_conf_rank"])
    else:
        tsv(pd.DataFrame(columns=["protein_id"]), qc / "duplicate_protein_ids.full.tsv")
        tsv(pd.DataFrame(columns=["protein_id", "duplicate_rows", "unique_signatures"]), qc / "duplicate_protein_ids.summary.tsv")
    clean = frame.drop(columns=TEMP_COLUMNS, errors="ignore")
    ordered = ["protein_id"] + sorted(set(clean.columns) - {"protein_id"})
    tsv(clean[ordered], tables / "tiered_candidates.tsv")
    for tier, group in clean.groupby("tier", dropna=False):
        tsv(group, tables / (str(tier).replace("/", "_").replace("\\", "_") + ".tsv"))
    prefixes = ["TIER_1B", "TIER_1C", "TIER_2A", "TIER_2B", "TIER_2C", "TIER_2D"]
    if not policy.skip_rescue_for_tier1a:
        prefixes.insert(0, "TIER_1A")
    if policy.include_tier3_in_rescue:
        prefixes.extend(["TIER_3A", "TIER_3B"])
    # Explicit ties make cap-boundary selection reproducible across pandas /
    # NumPy builds; the defined single-column quicksort did not define tie ordering.
    rescue = clean[clean["tier"].fillna("").str.contains("|".join(prefixes), regex=True)].sort_values(
        ["rescue_priority", "protein_id"], ascending=[False, True], kind="stable"
    )
    tsv(rescue, tables / "rescue_priority.tsv")
    def bed(data, path):
        with path.open("w") as stream:
            for _, row in data.iterrows():
                stream.write(f"{row['scaffold']}\t{int(row['start'])-1}\t{int(row['end'])}\t{row['protein_id']}|{row['tier']}\t{int(float(row['rescue_priority'])*10)}\t{row['strand']}\n")
        outputs.append(path)
    located = clean.dropna(subset=["scaffold", "start", "end", "strand"])
    bed(located, beds / "nlr_candidates.bed")
    for tier, group in located.groupby("tier", dropna=False):
        bed(group, beds / "by_tier" / (str(tier).replace("/", "_").replace("\\", "_") + ".bed"))
    tsv(clean["tier"].value_counts().to_frame("count"), summaries / "tier_summary.tsv", index=True)
    flags = clean["flags"].fillna("").str.split(";").explode().replace({"PASS": pd.NA, "": pd.NA}).dropna().value_counts()
    tsv(flags.to_frame("count"), summaries / "flag_counts.tsv", index=True)
    return outputs
