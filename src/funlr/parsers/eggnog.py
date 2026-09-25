"""eggNOG-mapper light parser (legacy stage 0.4)."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

EGGNOG_LIGHT_COLUMNS = ["protein_id", "COG_category", "Description", "Preferred_name", "GOs", "KEGG_ko"]


def parse_eggnog_light(eggnog_path: str | Path | None) -> pd.DataFrame:
    """Extract light columns from an emapper.annotations file.

    Fail-soft (legacy behavior): missing/empty input yields an empty frame with
    the expected columns.
    """
    rows = []
    if (
        eggnog_path
        and os.path.exists(eggnog_path)
        and os.path.getsize(eggnog_path) > 0
    ):
        with open(eggnog_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                rec = {"protein_id": cols[0]}
                if len(cols) > 6:
                    rec["COG_category"] = cols[6]
                if len(cols) > 7:
                    rec["Description"] = cols[7]
                if len(cols) > 8:
                    rec["Preferred_name"] = cols[8]
                if len(cols) > 9:
                    rec["GOs"] = cols[9]
                if len(cols) > 11:
                    rec["KEGG_ko"] = cols[11]
                rows.append(rec)
    eg = pd.DataFrame(rows)
    if eg.empty:
        eg = pd.DataFrame(columns=EGGNOG_LIGHT_COLUMNS)
    return eg.drop_duplicates()
