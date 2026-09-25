"""Export original-protein IDs from the strict table of a verified completed run.

Usage: python export_strict_proteins.py RUN/final_results/nlr_strict_candidates.tsv called.tsv
The new output is created exclusively. Fusion models require a separate locus
evaluation and are excluded here. This does not verify the source run itself.
"""
import csv
from pathlib import Path
import sys


def export_strict_proteins(source, destination):
    source = Path(source)
    if source.name != "nlr_strict_candidates.tsv":
        raise ValueError("Use nlr_strict_candidates.tsv, not the all-candidate final report")
    selected = []
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not {"protein_id", "is_fusion_model"}.issubset(reader.fieldnames):
            raise ValueError("Strict table is missing protein_id or is_fusion_model")
        for row in reader:
            flag = row["is_fusion_model"].strip().lower()
            if flag not in {"true", "false", "1", "0"}:
                raise ValueError("Strict table has an invalid is_fusion_model value")
            if flag in {"false", "0"}:
                rid = row["protein_id"]
                if not rid or any(c.isspace() for c in rid) or rid in selected:
                    raise ValueError("Strict table contains missing, invalid or duplicate protein IDs")
                selected.append(rid)
    with Path(destination).open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["protein_id"])
        writer.writerows([rid] for rid in selected)
    return len(selected)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: export_strict_proteins.py nlr_strict_candidates.tsv NEW-called.tsv")
    try:
        count = export_strict_proteins(sys.argv[1], sys.argv[2])
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    print(f"Exported {count} original-protein strict calls; fusion models excluded")
