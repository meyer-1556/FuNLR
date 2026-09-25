"""GFF3 parsers ported from the legacy pipeline (stages 0, 5, 6)."""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

TRANSCRIPT_COLUMNS = ["protein_id", "transcript_id", "gene_id", "scaffold", "start", "end", "strand"]
MINIPROT_FEATURE_COLUMNS = [
    "protein_id", "scaffold", "start", "end", "strand", "feature", "source", "score", "phase", "attrs_raw"
]
EXONERATE_FEATURE_COLUMNS = [
    "protein_id", "scaffold", "start", "end", "strand", "feature", "score", "source"
]


def parse_gff3_attributes(attr_str: str) -> dict[str, str]:
    """Legacy attribute parser: split on ';', then 'key=value'."""
    d: dict[str, str] = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def read_transcript_coords(gff_path: str | Path) -> pd.DataFrame:
    """Parse mRNA/transcript features -> coordinates table (legacy stage 0.3)."""
    rows = []
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype not in ("mRNA", "transcript"):
                continue
            a = parse_gff3_attributes(attrs)
            tid = a.get("ID", "")
            gid = a.get("Parent", "")
            if not tid:
                continue
            rows.append(
                {
                    "protein_id": tid,
                    "transcript_id": tid,
                    "gene_id": gid,
                    "scaffold": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                }
            )
    return pd.DataFrame(rows, columns=TRANSCRIPT_COLUMNS).drop_duplicates()


_ATTR_RE = re.compile(r"([A-Za-z0-9_]+)=([^;]+)")


def _attrs_regex(s: str) -> dict[str, str]:
    """Regex-based attribute scan used by the legacy miniprot parser."""
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(s)}


def read_miniprot_features(gff_path: str | Path) -> pd.DataFrame:
    """Parse miniprot GFF3 into a feature table (legacy stage 5.2)."""
    rows = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, src, ft, start, end, score, strand, phase, attr = parts
            a = _attrs_regex(attr)
            q = ""
            if "Target" in a:
                q = a["Target"].split()[0]
            elif "query" in a:
                q = a["query"]
            elif "Q" in a:
                q = a["Q"]
            elif "ID" in a:
                q = a["ID"]
            rows.append(
                {
                    "protein_id": q,
                    "scaffold": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "feature": ft,
                    "source": src,
                    "score": score,
                    "phase": phase,
                    "attrs_raw": attr,
                }
            )
    return pd.DataFrame(rows, columns=MINIPROT_FEATURE_COLUMNS)


def miniprot_locus_summary(features: pd.DataFrame) -> pd.DataFrame:
    """One locus per protein: min(start), max(end) across features (legacy 5.2)."""
    return features.groupby("protein_id", dropna=False).agg(
        scaffold=("scaffold", "first"),
        start=("start", "min"),
        end=("end", "max"),
        strand=("strand", "first"),
        n_features=("feature", "size"),
    ).reset_index()


_EXON_PID_RE = re.compile(r"(query|Query|protein|prot|sequence|seq)[ _:=]+([A-Za-z0-9_.-]+)")
_EXON_FALLBACK_RE = re.compile(r"([A-Za-z0-9_.-]+)")


def read_exonerate_features(gff_path: str | Path) -> pd.DataFrame:
    """Best-effort parse of exonerate --showtargetgff output (legacy stage 6.3)."""
    rows = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, _phase, attrs = parts
            pid = ""
            tagged = parse_gff3_attributes(attrs).get("FUNLR_query")
            m = _EXON_PID_RE.search(attrs)
            if tagged is not None:
                pid = unquote(tagged)
            elif m:
                pid = m.group(2)
            else:
                m2 = _EXON_FALLBACK_RE.search(attrs)
                pid = m2.group(1) if m2 else ""
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            rows.append(
                {
                    "protein_id": pid,
                    "scaffold": chrom,
                    "start": s,
                    "end": e,
                    "strand": strand,
                    "feature": feature,
                    "score": score,
                    "source": source,
                }
            )
    return pd.DataFrame(rows, columns=EXONERATE_FEATURE_COLUMNS)


def exonerate_locus_summary(features: pd.DataFrame) -> pd.DataFrame:
    """One locus per protein (legacy stage 6.3)."""
    return features.groupby("protein_id", dropna=False).agg(
        scaffold=("scaffold", "first"),
        start=("start", "min"),
        end=("end", "max"),
        strand=("strand", "first"),
        n_features=("feature", "size"),
    ).reset_index()
