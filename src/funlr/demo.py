"""Bundled real-tool installation test with reproducible synthetic inputs.

The synthetic HMMs are created locally with hmmbuild. They are test fixtures,
not fungal discovery databases, and must not be used for research analyses.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import random
import shutil
import sys

from . import __version__
from .doctor import check_installation
from .runner import invoke, sha256, write_json


_PURPOSE = "Synthetic real-tool installation test; NOT biological validation."


def _write_inputs(inputs: Path) -> list[tuple[str, Path]]:
    """Create deterministic sequences and alignments, using the original smoke fixture."""
    inputs.mkdir()
    rng = random.Random(717)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def peptide(length):
        return "".join(rng.choice(alphabet) for _ in range(length))

    nbd, repeat, spacer = peptide(220), peptide(60), peptide(80)
    proteins = {
        "synthetic_canonical": nbd + spacer + repeat * 4,
        "synthetic_rescue": nbd + peptide(90),
        "synthetic_absent": nbd + peptide(110),
    }
    codons = dict(zip(alphabet, "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))
    with (inputs / "proteins.faa").open("w") as fasta, (inputs / "genome.fa").open("w") as genome, (inputs / "genes.gff3").open("w") as gff:
        gff.write("##gff-version 3\n")
        for index, (name, sequence) in enumerate(proteins.items()):
            fasta.write(f">{name}\n{sequence}\n")
            start = 20001 if index == 0 else 501
            bases = "".join(rng.choice("ACGT") for _ in range(60000))
            # The third model intentionally has no matching genomic sequence.
            if index < 2:
                cds = "".join(codons[aa] for aa in sequence)
                bases = bases[:start - 1] + cds + bases[start - 1 + len(cds):]
            genome.write(f">chr{index}\n")
            for offset in range(0, len(bases), 60):
                genome.write(bases[offset:offset + 60] + "\n")
            gff.write(f"chr{index}\tSYNTHETIC\tmRNA\t{start}\t{start + len(sequence)*3-1}\t.\t+\t.\tID={name};Parent=gene{index}\n")
    alignments = []
    for model, sequence in (("NACHT", nbd), ("WD40", repeat)):
        alignment = inputs / (model + ".afa")
        alignment.write_text(f">synthetic_training_sequence\n{sequence}\n")
        alignments.append((model, alignment))
    return alignments


def _read_table(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_demo(outdir: str | Path) -> dict:
    """Run all eight stages, a positive Exonerate check, and verified resume.

    Requires a new output directory and checks the complete installation before
    creating it. Returns the saved evidence report. Every command uses actual
    executables, and packaged invocation never depends on a source checkout.
    """
    requested = Path(outdir).expanduser()
    if requested.exists() or requested.is_symlink():
        raise ValueError(f"Demo output must be a new directory: {requested}")
    root = requested.resolve()
    if any(char in str(root) for char in "\n\r\x00"):
        raise ValueError("Output paths cannot contain newlines or NUL characters")
    diagnosis = check_installation(include_demo=True)
    if not diagnosis["ok"]:
        raise ValueError("Demo requires a complete installation: " + "; ".join(diagnosis["errors"]))
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / "inputs"
    logs = root / "logs"
    logs.mkdir()
    evidence = {"status": "running", "purpose": _PURPOSE, "funlr_version": __version__,
                "output_dir": str(root), "installation": diagnosis}
    evidence_path = root / "evidence.json"
    write_json(evidence_path, evidence)
    print(_PURPOSE, flush=True)
    print(f"Demo logs and results: {root}", flush=True)
    env = os.environ.copy()
    # Use this exact installation, including with a relative development PYTHONPATH.
    env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent) + os.pathsep + env.get("PYTHONPATH", "")
    env.update(LC_ALL="C", PYTHONHASHSEED="0")
    command_log = logs / "commands.jsonl"

    def execute(command, name):
        invoke(command, logs / name, env, root, command_log)

    try:
        alignments = _write_inputs(inputs)
        for model, alignment in alignments:
            execute([diagnosis["tools"]["hmmbuild"]["path"], "--amino", "-n", model,
                     str(inputs / (model + ".hmm")), str(alignment)], model + "_hmmbuild.log")
        (inputs / "pfam_synthetic.hmm").write_bytes((inputs / "NACHT.hmm").read_bytes() + (inputs / "WD40.hmm").read_bytes())
        demo_config = root / "demo-config.json"
        write_json(demo_config, {"asm": {"asm_enable": 0}, "fusion": {"enabled": 0},
                                "reporting": {"integration": 0},
                                "discovery": {"pfam_nbd_scan": 0}})
        command = [sys.executable, "-m", "funlr", "run", "--config", str(demo_config), "--genome", str(inputs / "genome.fa"),
                   "--proteins", str(inputs / "proteins.faa"), "--gff3", str(inputs / "genes.gff3"),
                   "--nbd-hmms", str(inputs / "NACHT.hmm"), "--pfam", str(inputs / "pfam_synthetic.hmm"),
                   "--outdir", str(root / "run"), "--threads", "1", "--sample", "SYNTHETIC_INSTALLATION_TEST"]
        print("Running stages 0–7 with real bioinformatics tools…", flush=True)
        execute(command, "pipeline.log")
        state_before = json.loads((root / "run/state.json").read_text())
        stages = state_before.get("stages", {})
        if set(stages) != {str(number) for number in range(8)} or any(item.get("status") != "complete" for item in stages.values()):
            raise RuntimeError("Installation test did not complete all eight stages")
        final = {row["protein_id"]: row for row in _read_table(root / "run/final_results/nlr_final_report.tsv")}
        if final.get("synthetic_canonical", {}).get("tier") != "TIER_1A_HIGH_CONFIDENCE":
            raise RuntimeError("Real HMMER scans did not yield the designed canonical architecture")
        if not final.get("synthetic_rescue", {}).get("mp_scaffold"):
            raise RuntimeError("Real miniprot did not map the designed rescue query")

        # miniprot finds this query in the full run, so use isolated copied
        # handoffs to exercise positive Exonerate output without changing it.
        print("Checking positive Exonerate alignment and query identity…", flush=True)
        direct = root / "positive_exonerate"
        (direct / "results").mkdir(parents=True)
        for number in (0, 3, 4, 5):
            shutil.copytree(root / f"run/results/stage{number}", direct / f"results/stage{number}")
        (direct / "results/stage5/queries/exonerate_refine_ids.txt").write_text("synthetic_rescue\n")
        runtime = json.loads((root / "run/work/runtime.json").read_text())
        runtime["config"]["execution"]["output_dir"] = str(direct)
        settings = direct / "runtime.json"
        write_json(settings, runtime)
        execute([sys.executable, "-m", "funlr.stage_worker", "6", str(settings)], "positive_exonerate.log")
        exonerate_rows = _read_table(direct / "results/stage6/tables/exonerate_summary.tsv")
        if not any(row["protein_id"] == "synthetic_rescue" and int(row["n_features"]) > 0 for row in exonerate_rows):
            raise RuntimeError("Real Exonerate output did not preserve the synthetic query ID")
        print("Verifying resume reuses all eight completed stages…", flush=True)
        tool_commands_before = sha256(root / "run/logs/commands.jsonl")
        execute([*command, "--resume"], "resume.log")
        state_after = json.loads((root / "run/state.json").read_text())
        if state_after != state_before or sha256(root / "run/logs/commands.jsonl") != tool_commands_before:
            raise RuntimeError("Resume did not reuse the completed run unchanged")
        evidence.update(status="passed", completed_stages=8, resume_reused_stages=8,
                        full_run_manifest=str(root / "run/manifest.json"),
                        final_report=str(root / "run/final_results/nlr_final_report.tsv"),
                        positive_exonerate_rows=exonerate_rows,
                        evidence_file=str(evidence_path))
        write_json(evidence_path, evidence)
        print(f"Installation test passed. Evidence: {evidence_path}", flush=True)
        return evidence
    except BaseException as exc:
        evidence.update(status="failed", error=str(exc) or type(exc).__name__)
        write_json(evidence_path, evidence)
        raise
