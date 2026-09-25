"""Stage 6: exonerate refinement of hard rescue cases (port of 06_rescue_exonerate.sh).

Runs exonerate (protein2genome) per query for ids that miniprot did not place,
merges the per-query GFFs, and summarizes loci.
"""

from __future__ import annotations

from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from funlr.core.context import RunContext
from funlr.core.errors import StageError
from funlr.core.tools import require_tools
from funlr.parsers.fasta import fasta_ids
from funlr.parsers.gff import exonerate_locus_summary, read_exonerate_features

STAGE_NUMBER = 6
STAGE_NAME = "exonerate"
REQUIRED_TOOLS = ["seqkit", "exonerate"]

# Empty-frame column sets (legacy stage 6.3, lines 2038-2040).
SUMMARY_COLUMNS = ["protein_id", "scaffold", "start", "end", "strand", "n_features"]
COUNTS_COLUMNS = ["feature", "count"]


# ---------------------------------------------------------------------------
# Pure logic (no RunContext; unit-testable)
# ---------------------------------------------------------------------------

def apply_exon_maxn(ids: list[str], exon_maxn: int) -> list[str]:
    """Cap the refine list (legacy lines 1882-1891): 0 => all, else first N."""
    if exon_maxn == 0:
        return list(ids)
    return list(ids)[:exon_maxn]


def read_id_list(path: Path) -> list[str]:
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def write_id_list(ids: list[str], out_path: Path) -> None:
    with open(out_path, "w") as out:
        for pid in ids:
            out.write(pid + "\n")


def merge_gffs(per_query_gffs: list[tuple[str, Path]], out_path: Path) -> int:
    """Merge feature rows from this run and attach their known query identity.

    Raw Exonerate stdout remains in the per-query files. Alignment prose is not
    GFF; original attributes are escaped and preserved in ``LegacyAttrs``.
    """
    n_features = 0
    with open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for pid, gff in per_query_gffs:
            with open(gff) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 9:
                        continue
                    try:
                        int(fields[3])
                        int(fields[4])
                    except ValueError:
                        continue
                    fields[8] = "FUNLR_query=" + quote(pid, safe="") + ";LegacyAttrs=" + quote(fields[8], safe="")
                    out.write("\t".join(fields) + "\n")
                    n_features += 1
    return n_features


def write_manifest(paths: "OrderedDict[str, Path]", out_path: Path) -> None:
    with open(out_path, "w") as out:
        out.write("label\tpath\n")
        for label, path in paths.items():
            out.write(f"{label}\t{path}\n")


def _write_empty_parse_outputs(summary_tsv: Path, counts_tsv: Path) -> None:
    pd.DataFrame(columns=SUMMARY_COLUMNS).to_csv(summary_tsv, sep="\t", index=False)
    pd.DataFrame(columns=COUNTS_COLUMNS).to_csv(counts_tsv, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(ctx: RunContext) -> list[Path]:
    from funlr.stages.rescue_v3 import run_exonerate
    return run_exonerate(ctx)
