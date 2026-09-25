"""Domain architecture plots (native matplotlib port of the v3 R/ggplot2 block).

This module replaces the embedded ``Rscript`` plotting section of the v3 bash
stage 7 (``07_finalize_report.sh``) and the standalone
``plot_nlr_domains.sh`` Slurm script. The standalone script is an older,
simpler R variant (dplyr/ggplot2/RColorBrewer only, 14x10 facet figure,
newline-joined architecture summary); it is **superseded** by the richer
stage-7 embedded block ported here — Effector/NBD/Sensor/ASM/Other domain
categories, HCL-derived per-category fills, five summary charts, a 16x12
tier-faceted overview, per-protein track PNGs, and a counted
``domain_architecture_summary.tsv``.

Includes contributed rendering code under BSD-3-Clause; see
THIRD_PARTY_LICENSE.txt. Domain classification and architecture counts retain
the documented rules. The renderer uses the true Stage 3 protein length in the
scatter, preserves both legends, and validates plot coordinates and names.
The default renderer needs no R installation; the original R backend remains
available for comparison. Images are not claimed to be pixel-identical.

Public API (consumed by ``funlr.stages.stage7_finalize``)::

    generate_domain_plots(all_domain_hits_tsv, tiered_candidates_tsv, outdir, log) -> list[Path]
"""

from __future__ import annotations

import logging
import hashlib
import platform
import math
import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering; must precede pyplot import

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import PercentFormatter

__all__ = ["generate_domain_plots", "classify_domain_type", "plotting_provenance"]

DPI = 300

# Columns the module actually needs from each input table.
HITS_REQUIRED_COLUMNS = ("protein_id", "domain", "start", "end")
TIERED_REQUIRED_COLUMNS = ("protein_id", "tier", "protein_length")

# awk '$tiercol ~ /^TIER_[123]/' in the bash: Tier 1-3 candidates only.
# FUSION_RESCUE rows never match this pattern (and fusion models have no
# domain hits), so fusions are excluded by construction.
TIER_PATTERN = re.compile(r"^TIER_[123]")

# R `tier_order` vector: canonical display order for tiers, intersected with
# the tiers actually present in the data.
TIER_ORDER = [
    "TIER_1A_HIGH_CONFIDENCE",
    "TIER_1B_NEEDS_REVIEW",
    "TIER_1C_CANONICAL_NON_SSFR_SENSOR",
    "TIER_2A_HIGH_PRIORITY_RESCUE",
    "TIER_2B_RESCUE_CANDIDATE",
    "TIER_2C_LOW_PRIORITY_FRAGMENT",
    "TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE",
    "TIER_3A_INTEGRATED_DECOY",
    "TIER_3B_ARCHITECTURAL_VARIANT",
]

# Truthy strings for the optional is_fusion_model column in tiered_candidates.
_TRUTHY = {"true", "1", "yes", "t"}

# ---------------------------------------------------------------------------
# Domain token cascade (R `case_when`, first match wins, ignore.case = TRUE)
# ---------------------------------------------------------------------------

# gsub("__.*", "", domain): strip taxonomic suffixes like
# "NB-ARC__class__Agaricomycetes" before classification.
_SUFFIX_RE = re.compile(r"__.*")

# (label, include regex, exclude regex or None) — evaluated in order.
_DOMAIN_RULES: list[tuple[str, re.Pattern[str], re.Pattern[str] | None]] = [
    # ---------- Effector domains ----------
    ("Goodbye", r"Goodbye", None),
    ("HeLo", r"HeLo", r"HeLo-like"),
    ("HeLo-like", r"HeLo-like|HELL", None),
    ("HET", r"HET", r"HET-s"),
    ("HET-s", r"HET-s", None),
    ("TIR", r"TIR", None),
    ("Patatin", r"Patatin", None),
    ("PNP_UDP", r"PNP_UDP|PUP", None),
    ("RelA_SpoT", r"RelA_SpoT", None),
    ("Ses", r"Ses[AB]|ses[ab]", None),
    ("CHAT", r"CHAT", None),
    ("Crinkler", r"Crinkler", None),
    ("SAM", r"SAM", None),
    ("C2", r"C2[ -]domain|PF00168", None),
    ("Peptidase_S8", r"Peptidase_S8", None),
    ("CARD/Pyrin", r"CARD|Pyrin", None),
    ("CC", r"CC[_-]|Coiled[_-]coil", None),
    # ---------- NBD domains ----------
    ("NACHT", r"NACHT", None),
    ("NB-ARC", r"NB-ARC", None),
    ("AAA", r"AAA_16|AAA_22|AAA", None),  # collapse AAA variants
    # ---------- Sensor domains ----------
    ("WD40", r"WD40|WD_40|WD-repeat", None),
    ("ANK", r"ANK|Ank|ANKYRIN", None),
    ("TPR", r"TPR", None),
    ("HEAT", r"HEAT", None),
    ("LRR", r"LRR|Leucine[ -]rich", None),
    ("SPRY", r"SPRY", None),
    ("Kinase", r"PKinase|Kinase", None),
    ("Zinc_finger", r"C2H2|Zinc finger|ZINC_FINGER", None),
    ("ZZ", r"ZZ[_-]|ZZ-type", None),
    ("HMA/WRKY/LIM", r"HMA|WRKY|LIM", None),
    # ---------- ASM domains ----------
    ("HRAM", r"HRAM", None),
    ("PP", r"PP[ -]motif|NLR07|NLR39", None),
    ("sigma", r"sigma", None),  # also effector
    ("PUASM", r"PUASM|NLR32", None),
    ("Basidio_ASM", r"NLR05|NLR08|NLR22|NLR29|NLR44", None),
    ("BASS", r"BASS", None),
    ("Lineage_ASM", r"NLR17|NLR19|NLR34", None),
]
_DOMAIN_RULES = [
    (label, re.compile(pat, re.IGNORECASE), re.compile(neg, re.IGNORECASE) if neg else None)
    for label, pat, neg in _DOMAIN_RULES
]

# R `domain_category` case_when: domain type -> category (border color).
DOMAIN_CATEGORY: dict[str, str] = {
    **{t: "Effector" for t in (
        "Goodbye", "HeLo", "HeLo-like", "HET", "HET-s", "TIR", "Patatin",
        "PNP_UDP", "RelA_SpoT", "Ses", "CHAT", "Crinkler", "SAM", "C2",
        "Peptidase_S8", "CARD/Pyrin", "CC",
    )},
    **{t: "NBD" for t in ("NACHT", "NB-ARC", "AAA")},
    **{t: "Sensor" for t in (
        "WD40", "ANK", "TPR", "HEAT", "LRR", "SPRY", "Kinase",
        "Zinc_finger", "ZZ", "HMA/WRKY/LIM",
    )},
    **{t: "ASM" for t in (
        "HRAM", "PP", "sigma", "PUASM", "Basidio_ASM", "BASS", "Lineage_ASM",
    )},
}

# R `border_colors`.
BORDER_COLORS: dict[str, str] = {
    "Effector": "#E41A1C",
    "NBD": "#377EB8",
    "Sensor": "#4DAF4A",
    "ASM": "#984EA3",
    "Other": "#999999",
}

# R `hue_ranges` for the HCL-derived per-category fill colors.
_HUE_RANGES: dict[str, tuple[float, float]] = {
    "Effector": (0.0, 30.0),
    "NBD": (200.0, 260.0),
    "Sensor": (80.0, 140.0),
    "ASM": (260.0, 320.0),
    "Other": (0.0, 0.0),
}

_MARKER_CYCLE = ["o", "s", "^", "D", "v", "<", ">", "P", "X", "*"]


# ---------------------------------------------------------------------------
# Domain classification
# ---------------------------------------------------------------------------

def classify_domain_type(domain: str) -> str:
    """Simplify one domain token via the ordered regex cascade.

    Port of the R ``case_when`` block: taxonomic suffixes (``__...``) are
    stripped first, then the first matching rule wins; unmatched tokens map
    to ``"Other"``.
    """
    clean = _SUFFIX_RE.sub("", str(domain))
    for label, pat, neg in _DOMAIN_RULES:
        if pat.search(clean) and not (neg is not None and neg.search(clean)):
            return label
    return "Other"


def domain_category(domain_type: str) -> str:
    """Map a simplified domain type to its category (R ``domain_category``)."""
    return DOMAIN_CATEGORY.get(domain_type, "Other")


# ---------------------------------------------------------------------------
# HCL -> sRGB (R grDevices::hcl, polar CIELUV, D65, fixup = TRUE)
# ---------------------------------------------------------------------------

def _hcl(h: float, c: float, lum: float) -> str:
    """Convert polar-CIELUV (h, c, l) to a clamped sRGB hex color."""
    if lum <= 0:
        return "#000000"
    hrad = math.radians(h % 360.0)
    u = c * math.cos(hrad)
    v = c * math.sin(hrad)
    # D65 reference white, Yn = 100
    xn, yn, zn = 95.047, 100.0, 108.883
    un = 4.0 * xn / (xn + 15.0 * yn + 3.0 * zn)
    vn = 9.0 * yn / (xn + 15.0 * yn + 3.0 * zn)
    up = u / (13.0 * lum) + un
    vp = v / (13.0 * lum) + vn
    if lum > 8.0:
        y = yn * ((lum + 16.0) / 116.0) ** 3
    else:
        y = yn * lum * (3.0 / 29.0) ** 3
    x = y * 9.0 * up / (4.0 * vp)
    z = y * (12.0 - 3.0 * up - 20.0 * vp) / (4.0 * vp)
    x, y, z = x / 100.0, y / 100.0, z / 100.0
    # XYZ -> linear sRGB
    r = 3.2406 * x - 1.5372 * y - 0.4986 * z
    g = -0.9689 * x + 1.8758 * y + 0.0415 * z
    b = 0.0557 * x - 0.2040 * y + 1.0570 * z

    def _gamma(chan: float) -> float:
        if chan <= 0.0031308:
            return 12.92 * chan
        return 1.055 * chan ** (1.0 / 2.4) - 0.055

    # fixup = TRUE: clamp out-of-gamut channels
    rgb = [min(1.0, max(0.0, _gamma(chan))) for chan in (r, g, b)]
    return "#" + "".join(f"{int(round(chan * 255)):02X}" for chan in rgb)


def _fill_colors(domain_types: list[str]) -> dict[str, str]:
    """HCL-derived fill color per domain type (R ``get_category_colors``).

    Types are grouped by category (first-occurrence order, matching the R
    ``split(unique(...))`` semantics); within a category hues sweep the
    category range and luminance alternates 40/70. ``Other`` is flat grey.
    """
    by_cat: dict[str, list[str]] = {}
    for dtype in dict.fromkeys(domain_types):
        by_cat.setdefault(domain_category(dtype), []).append(dtype)
    colors: dict[str, str] = {}
    for cat, types in by_cat.items():
        if cat == "Other":
            colors.update({t: "#CCCCCC" for t in types})
            continue
        lo, hi = _HUE_RANGES[cat]
        n = len(types)
        # R seq(lo, hi, length.out = n): n == 1 yields lo
        hues = [lo + (hi - lo) * i / (n - 1) if n > 1 else lo for i in range(n)]
        for i, dtype in enumerate(types):
            lum = 40.0 if i % 2 == 0 else 70.0  # rep(c(40, 70), length.out = n)
            colors[dtype] = _hcl(hues[i], 80.0, lum)
    return colors


# ---------------------------------------------------------------------------
# Input loading / validation
# ---------------------------------------------------------------------------

def _read_tsv(path: Path, what: str, required: tuple[str, ...]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"{what} not found: {path}. Domain plots require this stage output; "
            f"run the upstream stages first."
        )
    if path.stat().st_size == 0:
        raise ValueError(f"{what} is empty: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{what} ({path}) is missing required column(s): "
            f"{', '.join(missing)}; found: {', '.join(df.columns)}"
        )
    return df


def _atomic_save_fig(fig: plt.Figure, path: Path) -> Path:
    """Render a PNG to a temp sibling then rename — never a partial PNG."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / (path.name + ".tmp")
    try:
        fig.savefig(tmp, format="png", dpi=DPI)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    finally:
        plt.close(fig)
    return path


def _atomic_write_text(text: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / (path.name + ".tmp")
    try:
        tmp.write_text(text)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return path


def _minimal_axes(ax: plt.Axes, rotate_x: bool = True) -> None:
    """theme_minimal()-ish styling: light grid, no top/right spines."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.set_axisbelow(True)
    if rotate_x:
        ax.tick_params(axis="x", labelrotation=45)
        for label in ax.get_xticklabels():
            label.set_ha("right")


def _cmap_colors(name: str) -> list[str]:
    return [matplotlib.colors.to_hex(c) for c in matplotlib.colormaps[name].colors]


# ---------------------------------------------------------------------------
# Per-protein feature aggregation (R section 5)
# ---------------------------------------------------------------------------

def _protein_features(domains: pd.DataFrame, protein_lengths: dict[str, float]) -> pd.DataFrame:
    """Aggregate domain hits to per-protein features (R ``protein_features``)."""
    rows = []
    for (pid, tier), g in domains.groupby(["protein_id", "tier"], sort=False):
        cat = g["domain_category"]
        nbd = g[cat == "NBD"]
        sensor = g[cat == "Sensor"]
        nbd_types = ";".join(dict.fromkeys(nbd["domain_type"]))
        sensor_types = ";".join(dict.fromkeys(sensor["domain_type"]))
        row = {
            "protein_id": pid,
            "tier": tier,
            "nbd_types": nbd_types,
            "nbd_start": float(nbd["start"].min()) if not nbd.empty else math.nan,
            "nbd_end": float(nbd["end"].max()) if not nbd.empty else math.nan,
            "sensor_types": sensor_types,
            "sensor_start": (
                float(sensor["start"].min()) if not sensor.empty else math.nan
            ),
            "has_asm": bool((cat == "ASM").any()),
            "has_effector": bool((cat == "Effector").any()),
            "protein_length": protein_lengths[pid],
        }
        # nbd_simple: prefer NACHT or NB-ARC, else AAA or "None"
        if "NACHT" in nbd_types:
            row["nbd_simple"] = "NACHT"
        elif "NB-ARC" in nbd_types:
            row["nbd_simple"] = "NB-ARC"
        elif "AAA" in nbd_types:
            row["nbd_simple"] = "AAA"
        else:
            row["nbd_simple"] = "None"
        # sensor_class: first matching class, else None/Other
        sensor_class = "Other"
        for klass in ("WD40", "ANK", "TPR", "HEAT", "LRR", "SPRY"):
            if klass in sensor_types:
                sensor_class = klass
                break
        else:
            sensor_class = "None" if sensor_types == "" else "Other"
        row["sensor_class"] = sensor_class
        # order_ok: NBD start < sensor start (with buffer 10); NA if either absent
        if not math.isnan(row["nbd_start"]) and not math.isnan(row["sensor_start"]):
            row["order_ok"] = bool(row["nbd_start"] + 10 < row["sensor_start"])
        else:
            row["order_ok"] = None
        rows.append(row)
    columns = [
        "protein_id", "tier", "nbd_types", "nbd_start", "nbd_end",
        "sensor_types", "sensor_start", "has_asm", "has_effector",
        "protein_length", "nbd_simple", "sensor_class", "order_ok",
    ]
    return pd.DataFrame(rows, columns=columns)


# ---------------------------------------------------------------------------
# Summary plots (R section 7)
# ---------------------------------------------------------------------------

def _stacked_proportion_bar(
    path: Path,
    features: pd.DataFrame,
    value_col: str,
    title: str,
    legend_title: str,
    color_map: dict[str, str],
    level_order: list[str],
    tier_levels: list[str],
) -> Path:
    """geom_bar(position = "fill") port: per-tier stacked proportions."""
    fig, ax = plt.subplots(figsize=(8, 5))
    tiers = [t for t in tier_levels if t in set(features["tier"])]
    x = list(range(len(tiers)))
    totals = {
        t: float((features["tier"] == t).sum()) for t in tiers
    }
    bottoms = [0.0] * len(tiers)
    for level in level_order:
        vals = []
        for t in tiers:
            n = float(((features["tier"] == t) & (features[value_col] == level)).sum())
            vals.append(n / totals[t] if totals[t] else 0.0)
        ax.bar(x, vals, bottom=bottoms, color=color_map[level], label=level,
               width=0.9, edgecolor="white", linewidth=0.5)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_xticks(x)
    ax.set_xticklabels(tiers)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.set_xlabel("Tier")
    ax.set_ylabel("Proportion")
    if level_order:
        ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.02, 0.5),
                  frameon=False)
    _minimal_axes(ax)
    fig.tight_layout()
    return _atomic_save_fig(fig, path)


def _sensor_plot(path: Path, features: pd.DataFrame, tier_levels: list[str]) -> Path:
    sub = features[features["sensor_class"] != "None"]
    levels = sorted(sub["sensor_class"].unique())  # ggplot alphabetical levels
    colors = dict(zip(levels, _cmap_colors("Set2")))  # scale_fill_brewer Set2
    return _stacked_proportion_bar(
        path, sub, "sensor_class",
        title="Sensor class distribution by tier",
        legend_title="Sensor class", color_map=colors, level_order=levels,
        tier_levels=tier_levels,
    )


def _nbd_plot(path: Path, features: pd.DataFrame, tier_levels: list[str]) -> Path:
    sub = features[features["nbd_simple"] != "None"]
    levels = sorted(sub["nbd_simple"].unique())
    colors = dict(zip(levels, _cmap_colors("Dark2")))  # scale_fill_brewer Dark2
    return _stacked_proportion_bar(
        path, sub, "nbd_simple",
        title="NBD subtype distribution by tier",
        legend_title="NBD type", color_map=colors, level_order=levels,
        tier_levels=tier_levels,
    )


def _asm_plot(path: Path, features: pd.DataFrame, tier_levels: list[str]) -> Path:
    sub = features.copy()
    sub["asm_present"] = sub["has_asm"].map({True: "ASM present", False: "No ASM"})
    fixed = {"ASM present": "#984EA3", "No ASM": "#CCCCCC"}
    # alphabetical levels, restricted to observed
    levels = [lvl for lvl in sorted(fixed) if lvl in set(sub["asm_present"])]
    colors = {lvl: fixed[lvl] for lvl in levels}
    return _stacked_proportion_bar(
        path, sub, "asm_present",
        title="ASM presence by tier",
        legend_title="", color_map=colors, level_order=levels,
        tier_levels=tier_levels,
    )


def _scatter_plot(path: Path, features: pd.DataFrame, tier_levels: list[str]) -> Path:
    sub = features[features["rescue_priority"].notna()]
    fig, ax = plt.subplots(figsize=(8, 5))
    tier_colors = dict(zip(tier_levels, _cmap_colors("Set1")))
    sensor_levels = sorted(sub["sensor_class"].unique())
    markers = {lvl: _MARKER_CYCLE[i % len(_MARKER_CYCLE)]
               for i, lvl in enumerate(sensor_levels)}
    for tier in tier_levels:
        tdf = sub[sub["tier"] == tier]
        for lvl in sensor_levels:
            pts = tdf[tdf["sensor_class"] == lvl]
            if pts.empty:
                continue
            ax.scatter(
                pts["protein_length"], pts["rescue_priority"],
                s=16, alpha=0.7, color=tier_colors.get(tier, "#333333"),
                marker=markers[lvl], edgecolors="none",
            )
    ax.set_title("Protein length vs rescue priority")
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Rescue priority")
    _minimal_axes(ax, rotate_x=False)
    tier_handles = [
        Line2D([], [], marker="o", linestyle="", color=tier_colors[t], label=t)
        for t in tier_levels if t in set(sub["tier"])
    ]
    marker_handles = [
        Line2D([], [], marker=markers[lvl], linestyle="", color="0.4", label=lvl)
        for lvl in sensor_levels
    ]
    # Explicitly retain the first legend when adding the second. Reserve space
    # inside the canvas: legends outside an unconstrained axes can be clipped.
    if tier_handles:
        tier_legend = ax.legend(
            handles=tier_handles, title="Tier", loc="upper left",
            bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=6.5,
            title_fontsize=8,
        )
        ax.add_artist(tier_legend)
    if marker_handles:
        ax.legend(handles=marker_handles, title="Sensor class",
                  loc="lower left", bbox_to_anchor=(1.02, 0.0),
                  frameon=False, fontsize=8, title_fontsize=8)
    fig.subplots_adjust(left=0.10, right=0.56, bottom=0.14, top=0.90)
    return _atomic_save_fig(fig, path)


def _order_plot(path: Path, features: pd.DataFrame, tier_levels: list[str]) -> Path:
    sub = features[features["order_ok"].notna()]
    tiers, fractions, labels = [], [], []
    for tier in tier_levels:
        tdf = sub[sub["tier"] == tier]
        if tdf.empty:
            continue
        total = len(tdf)
        ok = int(tdf["order_ok"].sum())
        tiers.append(tier)
        fractions.append(ok / total)
        labels.append(f"{ok}/{total}")
    fig, ax = plt.subplots(figsize=(8, 5))
    x = list(range(len(tiers)))
    ax.bar(x, fractions, color="#377EB8", width=0.9)
    for xi, frac, lab in zip(x, fractions, labels):
        ax.text(xi, frac + 0.02, lab, ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(tiers)
    ax.set_ylim(0, max([1.0] + [f * 1.15 for f in fractions]))
    ax.set_title("Order sanity (NBD before sensor)")
    ax.set_xlabel("Tier")
    ax.set_ylabel("Fraction with correct order")
    _minimal_axes(ax)
    fig.tight_layout()
    return _atomic_save_fig(fig, path)


# ---------------------------------------------------------------------------
# Per-protein track plots + tier-faceted overview
# ---------------------------------------------------------------------------

def _draw_domain_track(ax: plt.Axes, pdata: pd.DataFrame,
                       fill_colors: dict[str, str]) -> None:
    for _, hit in pdata.iterrows():
        ax.add_patch(Rectangle(
            (hit["start"], 0.0), hit["end"] - hit["start"], 1.0,
            facecolor=fill_colors.get(hit["domain_type"], "#CCCCCC"),
            edgecolor=BORDER_COLORS.get(hit["domain_category"], "#999999"),
            linewidth=1.0,
        ))
    ax.set_xlim(0, float(pdata["end"].max()) + 50)
    ax.set_ylim(0, 1)


def _individual_plots(outdir: Path, domains: pd.DataFrame,
                      fill_colors: dict[str, str], log: logging.Logger) -> list[Path]:
    """One 8x1.5in domain-track PNG per protein (R individual loop)."""
    ind_dir = outdir / "individual"
    ind_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    # unique(domains_filtered$protein_id): first-occurrence order
    for pid in dict.fromkeys(domains["protein_id"]):
        pdata = domains[domains["protein_id"] == pid]
        tier = pdata["tier"].iloc[0]
        fig, ax = plt.subplots(figsize=(8, 1.5))
        _draw_domain_track(ax, pdata, fill_colors)
        ax.axis("off")  # theme_void() + legend.position = "none"
        ax.set_title(f"{pid}\n{tier}", ha="center", fontsize=10)
        fig.tight_layout()
        paths.append(_atomic_save_fig(fig, ind_dir / f"{pid}.png"))
    log.info("Wrote %d individual protein plots to %s", len(paths), ind_dir)
    return paths


def _facet_plot(path: Path, domains: pd.DataFrame,
                fill_colors: dict[str, str], tier_levels: list[str]) -> Path:
    """16x12in tier-faceted overview: one row per protein (R ``p_summary``)."""
    facet = domains[domains["tier"].isin(tier_levels)]
    tiers = [t for t in tier_levels if t in set(facet["tier"])]
    present_types = [t for t in dict.fromkeys(facet["domain_type"])]
    present_cats = [c for c in ("Effector", "NBD", "Sensor", "ASM", "Other")
                    if c in set(facet["domain_category"])]

    if not tiers:
        fig, ax = plt.subplots(figsize=(16, 12))
        ax.text(0.5, 0.5, "No Tier 1-3 domain hits", ha="center", va="center",
                transform=ax.transAxes)
        ax.axis("off")
        return _atomic_save_fig(fig, path)

    fig, axes = plt.subplots(len(tiers), 1, figsize=(16, 12), squeeze=False)
    xmax = float(facet["end"].max()) * 1.05  # shared x scale (scales = "free_y")
    for ax, tier in zip(axes.ravel(), tiers):
        tdf = facet[facet["tier"] == tier]
        # arrange(tier, protein_id): proteins alphabetical within tier
        proteins = sorted(tdf["protein_id"].unique())
        for y, pid in enumerate(proteins, start=1):
            pdata = tdf[tdf["protein_id"] == pid]
            for _, hit in pdata.iterrows():
                ax.add_patch(Rectangle(
                    (hit["start"], y - 0.4), hit["end"] - hit["start"], 0.8,
                    facecolor=fill_colors.get(hit["domain_type"], "#CCCCCC"),
                    edgecolor=BORDER_COLORS.get(hit["domain_category"], "#999999"),
                    linewidth=0.6,
                ))
        ax.set_xlim(0, xmax)
        ax.set_yticks(range(1, len(proteins) + 1))
        ax.set_yticklabels(proteins, fontsize=7)
        ax.set_ylim(0.5, len(proteins) + 0.5)
        ax.margins(y=0.02)
        ax.set_title(tier, fontweight="bold", fontsize=12, loc="left")
        ax.grid(axis="y", visible=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes.ravel()[-1].set_xlabel("Amino Acid Position")
    for ax in axes.ravel():
        ax.set_ylabel("Protein ID")

    fill_handles = [
        Patch(facecolor=fill_colors.get(t, "#CCCCCC"), edgecolor="none", label=t)
        for t in present_types
    ]
    cat_handles = [
        Line2D([], [], marker="s", linestyle="", markersize=10,
               markerfacecolor="none", markeredgecolor=BORDER_COLORS[c], label=c)
        for c in present_cats
    ]
    leg1 = None
    if fill_handles:
        leg1 = fig.legend(handles=fill_handles, title="Domain Type", ncol=3,
                          loc="lower center", bbox_to_anchor=(0.5, 0.06),
                          frameon=False, fontsize="small", title_fontsize="medium")
        fig.add_artist(leg1)
    if cat_handles:
        fig.legend(handles=cat_handles, title="Category", ncol=5,
                   loc="lower center", bbox_to_anchor=(0.5, 0.01),
                   frameon=False, fontsize="small", title_fontsize="medium")
    fig.subplots_adjust(bottom=0.20 if fill_handles or cat_handles else 0.07,
                        hspace=0.5, left=0.15)
    return _atomic_save_fig(fig, path)


# ---------------------------------------------------------------------------
# domain_architecture_summary.tsv (R ``summary_stats``)
# ---------------------------------------------------------------------------

def _architecture_summary(domains: pd.DataFrame) -> pd.DataFrame:
    """Per-tier architecture counts: tier, total_proteins, architecture, count."""
    columns = ["tier", "total_proteins", "architecture", "count"]
    if domains.empty:
        return pd.DataFrame(columns=columns)
    arch = (
        domains.groupby(["tier", "protein_id"], sort=False)["domain_type"]
        .agg(lambda s: ";".join(sorted(set(s))))
        .reset_index(name="architecture")
    )
    counts = (
        arch.groupby(["tier", "architecture"], sort=False)
        .size().reset_index(name="count")
    )
    totals = (
        domains[["tier", "protein_id"]].drop_duplicates()
        .groupby("tier", sort=False).size().reset_index(name="total_proteins")
    )
    out = counts.merge(totals, on="tier", how="left")
    # R arrange(tier, desc(count), architecture) on character tier
    out = out.sort_values(
        ["tier", "count", "architecture"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    return out[columns].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def generate_domain_plots(
    all_domain_hits_tsv: Path,
    tiered_candidates_tsv: Path,
    outdir: Path,
    log: logging.Logger,
) -> list[Path]:
    """Generate the stage-7 domain architecture plot set.

    Parameters
    ----------
    all_domain_hits_tsv:
        ``stage2/all_domain_hits.tsv`` — needs ``protein_id, domain, start,
        end`` columns.
    tiered_candidates_tsv:
        ``stage3/tables/tiered_candidates.tsv`` — needs ``protein_id, tier,
        protein_length``;
        ``rescue_priority`` is used when present (skipped otherwise) and an
        optional truthy ``is_fusion_model`` column excludes fusion rows.
    outdir:
        Destination directory (e.g. ``results/stage7/domain_plots``); created
        if missing. PNGs/TSV are written atomically (temp file + rename), so
        a failure never leaves partial artifacts.
    log:
        Logger for progress/warnings.

    Returns
    -------
    list[Path]
        All files produced, in generation order. Empty when there are no
        Tier 1-3 candidates (plotting is skipped, mirroring the bash
        warn-and-skip behavior).

    Raises
    ------
    FileNotFoundError
        If an input TSV does not exist.
    ValueError
        If an input TSV is empty or lacks required columns.
    """
    all_domain_hits_tsv = Path(all_domain_hits_tsv)
    tiered_candidates_tsv = Path(tiered_candidates_tsv)
    outdir = Path(outdir)

    log.info("Reading domain hits...")
    hits = _read_tsv(all_domain_hits_tsv, "all_domain_hits.tsv", HITS_REQUIRED_COLUMNS)
    tiered = _read_tsv(
        tiered_candidates_tsv, "tiered_candidates.tsv", TIERED_REQUIRED_COLUMNS
    )

    # Tier 1-3 candidate IDs (bash awk '$tiercol ~ /^TIER_[123]/'); fusion
    # rows never match the tier pattern, and are additionally excluded via an
    # optional is_fusion_model column.
    is_tier123 = tiered["tier"].str.match(TIER_PATTERN)
    if "is_fusion_model" in tiered.columns:
        is_fusion = tiered["is_fusion_model"].str.strip().str.lower().isin(_TRUTHY)
        is_tier123 &= ~is_fusion
    candidates = tiered[is_tier123]
    for pid in candidates["protein_id"]:
        if (pid in {"", ".", ".."} or "/" in pid or "\\" in pid
                or any(ord(ch) < 32 or ord(ch) == 127 for ch in pid)):
            raise ValueError(f"Protein ID cannot be used as a plot filename: {pid!r}")
    if candidates["protein_id"].duplicated().any():
        raise ValueError("Tier 1-3 plot protein IDs must be unique")
    if candidates.empty:
        log.warning("No Tier 1-3 candidates found - skipping domain plots.")
        return []

    # tier_lookup <- setNames(tiered$tier, tiered$protein_id): first wins.
    tier_lookup: dict[str, str] = {}
    for pid, tier in zip(candidates["protein_id"], candidates["tier"]):
        tier_lookup.setdefault(pid, tier)

    domains = hits[hits["protein_id"].isin(tier_lookup)].copy()
    domains["start"] = pd.to_numeric(domains["start"], errors="coerce")
    domains["end"] = pd.to_numeric(domains["end"], errors="coerce")
    bad = (~domains["start"].map(math.isfinite)
           | ~domains["end"].map(math.isfinite)
           | (domains["start"] < 1)
           | (domains["end"] < domains["start"]))
    if bad.any():
        raise ValueError(
            f"{all_domain_hits_tsv}: {int(bad.sum())} row(s) have invalid "
            f"start/end coordinates"
        )
    domains["tier"] = domains["protein_id"].map(tier_lookup)
    domains["domain_type"] = domains["domain"].map(classify_domain_type)
    domains["domain_category"] = domains["domain_type"].map(domain_category)

    lengths = pd.to_numeric(candidates["protein_length"], errors="coerce")
    bad_length = ~lengths.map(math.isfinite) | (lengths <= 0)
    if bad_length.any():
        raise ValueError("Tier 1-3 protein_length must contain positive finite lengths")
    protein_lengths = dict(zip(candidates["protein_id"], lengths.astype(float)))
    features = _protein_features(domains, protein_lengths)
    # R factor(tier, levels = intersect(TIER_ORDER, present)): rows in tiers
    # outside TIER_ORDER are dropped from the tier-indexed plots.
    tier_levels = [t for t in TIER_ORDER if t in set(features["tier"])] \
        if not features.empty else []
    dropped = sorted(set(features["tier"]) - set(TIER_ORDER)) \
        if not features.empty else []
    if dropped:
        log.info("Tiers outside the canonical order are excluded from "
                 "tier-indexed plots: %s", ", ".join(dropped))
    features_plot = features[features["tier"].isin(tier_levels)].copy() \
        if not features.empty else features

    outdir.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []

    log.info("Generating summary plots...")
    produced.append(_sensor_plot(outdir / "sensor_distribution.png",
                                 features_plot, tier_levels))
    produced.append(_nbd_plot(outdir / "nbd_distribution.png",
                              features_plot, tier_levels))
    produced.append(_asm_plot(outdir / "asm_presence.png",
                              features_plot, tier_levels))

    if "rescue_priority" in tiered.columns:
        prio = tiered.drop_duplicates("protein_id")[
            ["protein_id", "rescue_priority"]
        ]
        prio = prio.assign(
            rescue_priority=pd.to_numeric(prio["rescue_priority"], errors="coerce")
        )
        features_plot = features_plot.merge(prio, on="protein_id", how="left")
        produced.append(_scatter_plot(outdir / "length_vs_priority.png",
                                      features_plot, tier_levels))
    else:
        log.warning(
            "tiered_candidates.tsv has no rescue_priority column - "
            "skipping length_vs_priority.png"
        )

    produced.append(_order_plot(outdir / "order_sanity.png",
                                features_plot, tier_levels))

    fill = _fill_colors(list(domains["domain_type"]))

    log.info("Generating individual protein plots...")
    produced.extend(_individual_plots(outdir, domains, fill, log))

    log.info("Generating tier summary figure...")
    produced.append(_facet_plot(outdir / "nlr_domains_by_tier.png",
                                domains, fill, tier_levels))

    summary = _architecture_summary(domains)
    produced.append(_atomic_write_text(
        summary.to_csv(sep="\t", index=False),
        outdir / "domain_architecture_summary.tsv",
    ))

    log.info("Domain plotting completed: %d artifacts in %s", len(produced), outdir)
    return produced


def plotting_provenance() -> dict:
    """Record the actual renderer and bundled font bytes used by the plots."""
    from matplotlib import font_manager
    font_files = {}
    for weight in ("normal", "bold"):
        properties = font_manager.FontProperties(
            family=matplotlib.rcParams["font.family"], weight=weight,
        )
        font = Path(font_manager.findfont(properties, fallback_to_default=True))
        font_files[weight] = {
            "name": font_manager.FontProperties(fname=str(font)).get_name(),
            "path": str(font.resolve()),
            "sha256": hashlib.sha256(font.read_bytes()).hexdigest(),
        }
    return {
        "backend": "matplotlib",
        "renderer": str(matplotlib.get_backend()),
        "matplotlib_version": matplotlib.__version__,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "dpi": DPI,
        "protein_length_source": "Stage3 protein_length",
        "font_family": list(matplotlib.rcParams["font.family"]),
        "font_sans_serif": list(matplotlib.rcParams["font.sans-serif"]),
        "font_files": font_files,
        "configured_cache": os.environ.get("MPLCONFIGDIR"),
        "actual_cache": matplotlib.get_cachedir(),
        "images_pixel_identical_to_R": False,
    }
