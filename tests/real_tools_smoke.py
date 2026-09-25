"""Opt-in real-executable integration smoke; synthetic sequences, not biology.

Run with the dependency environment active:
    python tests/real_tools_smoke.py --outdir /tmp/funlr-real-smoke
Requires hmmbuild in addition to the production executables. No downloads.
"""
import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root = Path(args.outdir).resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / "inputs"
    inputs.mkdir()
    rng = random.Random(717)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    peptide = lambda n: "".join(rng.choice(alphabet) for _ in range(n))
    nbd, repeat, spacer = peptide(220), peptide(60), peptide(80)
    canonical = nbd + spacer + repeat * 4
    rescue = nbd + peptide(90)
    absent = nbd + peptide(110)
    proteins = {"synthetic_canonical": canonical, "synthetic_rescue": rescue, "synthetic_absent": absent}
    codons = dict(zip(alphabet, "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))
    with (inputs / "proteins.faa").open("w") as fasta, (inputs / "genome.fa").open("w") as genome, (inputs / "genes.gff3").open("w") as gff:
        gff.write("##gff-version 3\n")
        for index, (name, sequence) in enumerate(proteins.items()):
            fasta.write(f">{name}\n{sequence}\n")
            start = 20001 if index == 0 else 501
            # The third supplied model deliberately lacks genomic homology;
            # valid coordinates do not imply a biologically valid gene model.
            bases = "".join(rng.choice("ACGT") for _ in range(60000))
            if index < 2:
                cds = "".join(codons[aa] for aa in sequence)
                bases = bases[:start - 1] + cds + bases[start - 1 + len(cds):]
            genome.write(f">chr{index}\n")
            for offset in range(0, len(bases), 60):
                genome.write(bases[offset:offset + 60] + "\n")
            gff.write(f"chr{index}\tSYNTHETIC\tmRNA\t{start}\t{start + len(sequence)*3-1}\t.\t+\t.\tID={name};Parent=gene{index}\n")
    for model, seq in (("NACHT", nbd), ("WD40", repeat)):
        alignment = inputs / (model + ".afa")
        alignment.write_text(f">synthetic_training_sequence\n{seq}\n")
        with (root / (model + "_hmmbuild.log")).open("w") as log:
            subprocess.run(["hmmbuild", "--amino", "-n", model, str(inputs / (model + ".hmm")), str(alignment)], stdout=log, stderr=subprocess.STDOUT, check=True)
    (inputs / "pfam_synthetic.hmm").write_bytes((inputs / "NACHT.hmm").read_bytes() + (inputs / "WD40.hmm").read_bytes())
    command = [sys.executable, "-m", "funlr", "run", "--genome", str(inputs / "genome.fa"),
               "--proteins", str(inputs / "proteins.faa"), "--gff3", str(inputs / "genes.gff3"),
               "--nbd-hmms", str(inputs / "NACHT.hmm"), "--pfam", str(inputs / "pfam_synthetic.hmm"),
               "--outdir", str(root / "run"), "--threads", "1", "--sample", "SYNTHETIC_REAL_TOOL_SMOKE"]
    subprocess.run(command, check=True)
    # Run real Exonerate refinement on a known matching query as a separate
    # positive format check; miniprot already finds it in the ordinary run.
    stage5 = root / "run/results/stage5"
    env = os.environ.copy()
    # Preserve full-run outputs and state. The direct-stage experiment receives
    # separate copies of only its required inputs.
    import shutil
    direct = root / "direct_exonerate"
    (direct / "results").mkdir(parents=True)
    shutil.copytree(root / "run/results/stage0", direct / "results/stage0")
    shutil.copytree(stage5, direct / "results/stage5")
    (direct / "results/stage5/queries/exonerate_refine_ids.txt").write_text("synthetic_rescue\n")
    runtime = json.loads((root / "run/work/runtime.json").read_text())
    runtime["config"]["execution"]["output_dir"] = str(direct)
    settings = direct / "runtime.json"
    settings.write_text(json.dumps(runtime))
    with (direct / "stage6.log").open("w") as log:
        subprocess.run([sys.executable, "-m", "funlr.stage_worker", "6", str(settings)], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    import csv
    with (root / "run/final_results/nlr_final_report.tsv").open() as handle:
        final = {row["protein_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if final["synthetic_canonical"]["tier"] != "TIER_1A_HIGH_CONFIDENCE":
        raise RuntimeError("Real HMMER scans did not yield the designed canonical architecture")
    if not final["synthetic_rescue"]["mp_scaffold"]:
        raise RuntimeError("Real miniprot did not map the designed rescue query")
    with (direct / "results/stage6/tables/exonerate_summary.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not any(row["protein_id"] == "synthetic_rescue" for row in rows):
        raise RuntimeError("Real Exonerate output did not preserve the synthetic query ID")
    subprocess.run([*command, "--resume"], check=True)
    evidence = {"kind": "synthetic real-executable smoke; NOT biological equivalence",
                "full_run_manifest": str(root / "run/manifest.json"), "positive_exonerate_rows": rows}
    (root / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(f"Real-tool synthetic smoke passed: {root}")


if __name__ == "__main__":
    main()
