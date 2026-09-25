"""Non-destructive reference annotation exports using native gffread.

These copies annotate existing gene models; they do not insert rescue models.
"""
from __future__ import annotations

import csv
from pathlib import Path
import re
from urllib.parse import quote, unquote

import pandas as pd

from funlr.core.errors import StageError

NLR_FIELDS = ["NLR_tier", "NLR_nbd_confidence", "NLR_rescue_priority", "NLR_flags", "NLR_domains", "NLR_architecture"]
REPORT_FIELDS = ["tier", "nbd_confidence", "rescue_priority", "flags", "domains_grouped", "NLR_architecture"]


def attrs(text):
    return {key: unquote(value) for field in text.split(";") if "=" in field for key, value in [field.split("=", 1)]}


def clean(value):
    return re.sub(r"\s+", " ", str(value or "").replace("[", "(").replace("]", ")")).strip()


def tagged_gff(source, target, gene_calls):
    """Retain gene-level propagation; replace managed tags, clearing stale calls."""
    with open(source) as inp, open(target, "w") as out:
        for line in inp:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9 or fields[2] != "mRNA":
                out.write(line)
                continue
            parsed = attrs(fields[8])
            call = gene_calls.get(parsed.get("Parent"))
            # A previously tagged input may now have no call for this gene.
            # Clear only our exact keys; similarly named attributes and encoded
            # values belong to the input annotation and must remain untouched.
            if call or any(key in parsed for key in NLR_FIELDS):
                kept = [x for x in fields[8].split(";") if x and x.split("=", 1)[0] not in NLR_FIELDS]
                kept += [key + "=" + quote(value, safe="._:-") for key, value in (call or {}).items()]
                fields[8] = ";".join(kept)
            out.write("\t".join(fields) + "\n")


def rewrite_protein_headers(gff, raw, output):
    metadata = {}
    with open(gff) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9 or fields[2] != "mRNA":
                continue
            a = attrs(fields[8])
            if a.get("ID"):
                product = next((a[k] for k in ("product", "Product", "description", "Description", "Name") if a.get(k)), "hypothetical protein")
                metadata[a["ID"]] = {"product": clean(product), "tier": clean(a.get("NLR_tier")), "arch": clean(a.get("NLR_architecture"))}
    with open(raw) as inp, open(output, "w") as out:
        for line in inp:
            if not line.startswith(">"):
                out.write(line)
                continue
            tid = line[1:].split()[0]
            m = metadata.get(tid, {"product": "hypothetical protein", "tier": "", "arch": ""})
            header = f">{tid} [product={m['product']}]"
            if m["tier"]:
                header += f" [NLR_tier={m['tier']}]"
            if m["arch"]:
                header += f" [NLR_architecture={m['arch']}]"
            out.write(header + "\n")


def export_annotations(ctx, report, destination):
    destination.mkdir(parents=True, exist_ok=True)
    cfg = ctx.config
    protein_calls, gene_calls, conflicts = {}, {}, []
    t1, t2 = set(), set()
    for row in report.fillna("").to_dict("records"):
        pid, gid = str(row["protein_id"]), str(row.get("gene_id", ""))
        call = {key: str(row.get(source, "")) for key, source in zip(NLR_FIELDS, REPORT_FIELDS)}
        protein_calls[pid] = call
        if gid:
            if gid in gene_calls and gene_calls[gid] != call:
                conflicts.append((gid, pid, "FIRST_REPORT_ISOFORM_RETAINED"))
            gene_calls.setdefault(gid, call)
            if str(row["tier"]).startswith("TIER_1"):
                t1.add(gid)
            if str(row["tier"]).startswith("TIER_2"):
                t2.add(gid)
    pd.DataFrame([{"gene_id": gid, **value} for gid, value in sorted(gene_calls.items())], columns=["gene_id", *NLR_FIELDS]).to_csv(destination / "nlr_calls.tsv", sep="\t", index=False)
    pd.DataFrame(conflicts, columns=["gene_id", "protein_id", "policy"]).to_csv(destination / "isoform_conflicts.tsv", sep="\t", index=False)
    for suffix, ids in (("T1", t1), ("T2", t2), ("T1T2", t1 | t2)):
        (destination / f"nlr_genes_{suffix}.txt").write_text("".join(key + "\n" for key in sorted(ids)))
    report["tier"].value_counts().rename_axis("tier").reset_index(name="count").to_csv(destination / "nlr_tier_counts.tsv", sep="\t", index=False)
    status = []
    annotation = cfg.get("inputs", "annotations")
    if annotation and Path(annotation).is_file():
        path = Path(annotation)
        output = destination / (path.stem + ".withNLR.tsv")
        with path.open(newline="") as inp, output.open("w", newline="") as out:
            reader = csv.DictReader(inp, delimiter="\t")
            if not reader.fieldnames or "GeneID" not in reader.fieldnames:
                raise StageError("Funannotate annotation table requires a GeneID header")
            writer = csv.DictWriter(out, fieldnames=[*reader.fieldnames, *(key for key in NLR_FIELDS if key not in reader.fieldnames)], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                row.update(gene_calls.get(row.get("GeneID"), {key: "" for key in NLR_FIELDS}))
                writer.writerow(row)
        status.append(("annotations", "WRITTEN", str(output)))
    else:
        status.append(("annotations", "SKIPPED_NO_INPUT", ""))
    gff = cfg.get("inputs", "gff3")
    tagged = None
    if cfg.get("reporting", "write_tagged_gff3", True) and gff and Path(gff).is_file():
        tagged = destination / (Path(gff).stem + ".withNLR.gff3")
        tagged_gff(gff, tagged, gene_calls)
        status.append(("gff3", "WRITTEN_EXISTING_MODELS", str(tagged)))
    else:
        status.append(("gff3", "SKIPPED", ""))
    eggnog = cfg.get("inputs", "eggnog")
    if cfg.get("reporting", "write_tagged_eggnog", True) and eggnog and Path(eggnog).is_file():
        output = destination / (Path(eggnog).stem + ".withNLR.annotations")
        with open(eggnog) as inp, output.open("w") as out:
            for line in inp:
                if line.startswith("#query"):
                    out.write(line.rstrip("\n") + "\t" + "\t".join(NLR_FIELDS) + "\n")
                elif line.startswith("#") or not line.strip():
                    out.write(line)
                else:
                    pid = line.split("\t", 1)[0]
                    values = protein_calls.get(pid, {})
                    out.write(line.rstrip("\n") + "\t" + "\t".join(values.get(key, "") for key in NLR_FIELDS) + "\n")
        status.append(("eggnog", "WRITTEN", str(output)))
    else:
        status.append(("eggnog", "SKIPPED", ""))
    if cfg.get("reporting", "write_proteins", True) and tagged:
        raw = destination / ("." + tagged.stem + ".proteins.raw.fa")
        output = destination / (tagged.stem + ".proteins.fa")
        ctx.runner.run([ctx.runner.resolve("gffread"), tagged, "-g", cfg.get("inputs", "genome"), "-y", raw], stderr_path=destination / "gffread.log")
        if not raw.is_file() or not raw.stat().st_size:
            raise StageError("gffread produced no translated proteins")
        rewrite_protein_headers(tagged, raw, output)
        raw.unlink()
        status.append(("proteins", "TRANSLATED_EXISTING_MODELS", str(output)))
    else:
        status.append(("proteins", "SKIPPED", ""))
    pd.DataFrame(status, columns=["export", "status", "path"]).to_csv(destination / "export_status.tsv", sep="\t", index=False)
    return sorted(path for path in destination.iterdir() if path.is_file())
