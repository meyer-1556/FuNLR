"""HMMER domtblout parsing shared by stages 1 and 2 (legacy-faithful)."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

_HIT_COLUMNS = [
    "protein_id", "domain", "acc", "dom_i_evalue", "start", "end", "ali_len", "source",
]


def parse_domtblout(
    path: str | Path,
    max_i_evalue: float,
    min_ali_len: int,
    source: str,
) -> tuple[pd.DataFrame, dict]:
    """Parse a domtblout file with the legacy hard filters.

    Filters: domain i-Evalue <= max_i_evalue AND aligned length >= min_ali_len.
    Aligned length = abs(ali_end - ali_start) + 1 from domtblout columns 18/19
    (0-based indices 17/18). Returns (hits_df, report_dict) where report_dict has
    the legacy counter names: total_domtbl_lines, kept_after_filters,
    filtered_high_i_eval, filtered_short_align (+ source).
    """
    rows = []
    total = kept = bade = short = 0
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        df = pd.DataFrame(columns=_HIT_COLUMNS)
        return df, {
            "source": source,
            "total_domtbl_lines": 0,
            "kept_after_filters": 0,
            "filtered_high_i_eval": 0,
            "filtered_short_align": 0,
        }
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            x = line.split()
            # domtblout has fixed first 22 columns; description may extend beyond
            if len(x) < 22:
                continue
            dom_name = x[0]  # target name
            dom_acc = x[1]
            prot = x[3]      # query name
            try:
                i_eval = float(x[12])   # domain i-Evalue
                ali_start = int(x[17])
                ali_end = int(x[18])
            except (ValueError, IndexError):
                continue
            alilen = abs(ali_end - ali_start) + 1
            if i_eval > max_i_evalue:
                bade += 1
                continue
            if alilen < min_ali_len:
                short += 1
                continue
            kept += 1
            rows.append((prot, dom_name, dom_acc, i_eval, ali_start, ali_end, alilen, source))
    df = pd.DataFrame(rows, columns=_HIT_COLUMNS)
    rep = {
        "source": source,
        "total_domtbl_lines": total,
        "kept_after_filters": kept,
        "filtered_high_i_eval": bade,
        "filtered_short_align": short,
    }
    return df, rep
