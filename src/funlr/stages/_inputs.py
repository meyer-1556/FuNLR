"""Protein/transcript ID and eggNOG reconciliation for Stage 0.

The source's first-transcript/first-CDS policy and coordinate-driven master
universe are retained. Headerless eggNOG is rejected explicitly rather than
failing through an uninitialized variable.
"""
from pathlib import Path

import pandas as pd

from funlr.core.errors import StageError

COORD_COLUMNS = ["protein_id", "transcript_id", "gene_id", "scaffold", "start",
                 "end", "strand", "mapping_source", "cds_protein_id"]
EGGNOG_COLUMNS = ["protein_id", "COG_category", "Description", "Preferred_name",
                  "GOs", "KEGG_ko", "PFAMs"]


def reconcile_coordinates(gff_path, proteins_path):
    with open(proteins_path) as handle:
        protein_ids = {line[1:].strip().split()[0] for line in handle if line.startswith(">")}
    tx_rows, tx_to_prot, seen_tx = [], {}, set()
    with open(gff_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end, _, strand, _, attrs = parts
            attr = dict(item.split("=", 1) for item in attrs.split(";") if "=" in item)
            if ftype in ("mRNA", "transcript"):
                tid = attr.get("ID", "")
                if not tid or tid in seen_tx:
                    continue
                seen_tx.add(tid)
                tx_rows.append({"transcript_id": tid, "gene_id": attr.get("Parent", ""),
                                "scaffold": seqid, "start": int(start), "end": int(end),
                                "strand": strand})
            elif ftype == "CDS":
                parent = attr.get("Parent", "")
                pid = attr.get("protein_id", "") or attr.get("orig_protein_id", "")
                if parent and pid:
                    tx_to_prot.setdefault(parent.split(",")[0], pid.split("|")[-1])
    counts = {"gff_transcripts": 0, "proteins_in_fasta": len(protein_ids),
              "cds_protein_id_mappings_raw": len(tx_to_prot),
              "assigned_via_cds_protein_id": 0, "assigned_via_transcript_id": 0,
              "assigned_via_gene_id": 0, "fallback_cds_protein_id_not_in_fasta": 0,
              "fallback_transcript_id_not_in_fasta": 0}
    rows = []
    for tx in tx_rows:
        tid, gid = tx["transcript_id"], tx["gene_id"]
        cds_pid = tx_to_prot.get(tid, "")
        if cds_pid and cds_pid in protein_ids:
            pid, source, counter = cds_pid, "CDS_PROTEIN_ID", "assigned_via_cds_protein_id"
        elif tid in protein_ids:
            pid, source, counter = tid, "TRANSCRIPT_ID_MATCH", "assigned_via_transcript_id"
        elif gid and gid in protein_ids:
            pid, source, counter = gid, "GENE_ID_MATCH", "assigned_via_gene_id"
        elif cds_pid:
            pid, source, counter = cds_pid, "CDS_PROTEIN_ID_NOT_IN_FASTA", "fallback_cds_protein_id_not_in_fasta"
        else:
            pid, source, counter = tid, "TRANSCRIPT_ID_NOT_IN_FASTA", "fallback_transcript_id_not_in_fasta"
        counts[counter] += 1
        rows.append({"protein_id": pid, **tx, "mapping_source": source, "cds_protein_id": cds_pid})
    coords = pd.DataFrame(rows, columns=COORD_COLUMNS)
    counts["gff_transcripts"] = len(coords)
    counts["rows_with_protein_id_not_in_fasta"] = len(set(coords["protein_id"].astype(str)) - protein_ids)
    return coords, counts


def parse_header_eggnog(path):
    """Read the original #query-bearing eggNOG table, retaining PFAM names."""
    if not path or not Path(path).is_file() or not Path(path).stat().st_size:
        return pd.DataFrame(columns=EGGNOG_COLUMNS)
    header = None
    with open(path) as handle:
        for line in handle:
            if line.startswith("#query"):
                header = {name.lstrip("#"): i for i, name in enumerate(line.strip().split("\t"))}
                break
    if header is None:
        raise StageError(f"eggNOG annotations need the original #query header: {path}")
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            row = {"protein_id": fields[0]}
            for name in EGGNOG_COLUMNS[1:]:
                i = header.get(name)
                row[name] = fields[i] if i is not None and i < len(fields) else ""
            if row["PFAMs"] in {"-", "", " "}:
                row["PFAMs"] = ""
            rows.append(row)
    return pd.DataFrame(rows, columns=EGGNOG_COLUMNS)


def build_updated_master_table(coords_path, lens_path, eggnog_path, eggnog_light_out, master_out):
    lens = pd.read_csv(lens_path, sep="\t")
    coords = pd.read_csv(coords_path, sep="\t")
    eggnog = parse_header_eggnog(eggnog_path)
    # reference writes deduplicated light rows but merges the original parsed rows.
    eggnog.drop_duplicates().to_csv(eggnog_light_out, sep="\t", index=False)
    master = coords.merge(lens, on="protein_id", how="left").merge(eggnog, on="protein_id", how="left")
    for name in EGGNOG_COLUMNS[1:]:
        if name not in master:
            master[name] = ""
    master["annotation_text"] = (
        master["Description"].fillna("").astype(str) + " " +
        master["Preferred_name"].fillna("").astype(str) + " " +
        master["GOs"].fillna("").astype(str) + " " +
        master["COG_category"].fillna("").astype(str) + " " +
        master["KEGG_ko"].fillna("").astype(str) + " " +
        master["PFAMs"].fillna("").astype(str)
    ).str.replace(r"\s+", " ", regex=True).str.strip()
    master["PFAMs"] = master["PFAMs"].fillna("").astype(str).replace(["", "-", " ", "--"], "")
    master.to_csv(master_out, sep="\t", index=False)
    return master, eggnog
