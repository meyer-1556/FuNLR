"""Verified public sequences in a deliberately constructed installation fixture."""
from __future__ import annotations

import csv
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path, PurePosixPath
import sys

from . import __version__
from .doctor import check_installation
from .runner import invoke, sha256, write_json

PURPOSE = "Public-sequence hybrid installation example; NOT reference parity or independent biological validation."
TIERS = {
    "HETE_AAL37299": "TIER_1A_HIGH_CONFIDENCE",
    "HETE_Q8X1P4": "TIER_1B_NEEDS_REVIEW",
    "HETE_A7IQW3": "TIER_1A_HIGH_CONFIDENCE",
    "KELCH_A0A090DCW7": "TIER_4A_REPEAT_ONLY_NO_NBD",
    "TPR_B2B7X0": "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL",
    "ATG1_Q3ZDQ4": "TIER_4A_REPEAT_ONLY_NO_NBD",
}
DEMO_SETTINGS = {"asm": {"asm_enable": 0}, "fusion": {"enabled": 0},
                 "discovery": {"discovery_mode": "BROAD", "pfam_nbd_scan": 0},
                 # The fixture has hand-made gene/mRNA spans, with no CDS.
                 "reporting": {"write_proteins": 0}}
INPUTS = {"genome.fa", "proteins.faa", "annotation.gff3", "emapper.annotations", "db/nbd.hmm", "db/pfam_mini.hmm"}


def verified_assets(resource_root=None):
    """Read and checksum the bundled fixture before creating an output directory."""
    root = resource_root or files("funlr").joinpath("data", "public_sequences")
    manifest = json.loads(root.joinpath("assets.json").read_text())
    if manifest.get("schema_version") != 1 or not INPUTS <= set(manifest.get("sha256", {})):
        raise ValueError("Public example asset manifest is incomplete")
    result = {}
    for name, expected in manifest["sha256"].items():
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or "\\" in name:
            raise ValueError("Invalid path in public example asset manifest")
        data = root.joinpath(*path.parts).read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Public example checksum mismatch: {name}")
        result[name] = data
    return result, manifest


def _table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def verify_outputs(run):
    rows = _table(run / "final_results/nlr_final_report.tsv")
    tiers = {row["protein_id"]: row["tier"] for row in rows}
    if len(rows) != len(TIERS) or tiers != TIERS:
        raise RuntimeError(f"Public example candidate/tier mismatch: {tiers}")
    if "HARD_END" not in next(row["flags"] for row in rows if row["protein_id"] == "HETE_Q8X1P4"):
        raise RuntimeError("Public example lost the designed contig-edge flag")
    strict = _table(run / "final_results/nlr_strict_candidates.tsv")
    if len(strict) != 3 or {row["protein_id"] for row in strict} != {"HETE_AAL37299", "HETE_Q8X1P4", "HETE_A7IQW3"}:
        raise RuntimeError("Public example strict-set mismatch")
    # Inspect the raw alignments, because the retained legacy summary collapses
    # different scaffolds and cannot establish that both inserted loci mapped.
    loci = {}
    for line in (run / "results/stage5/gff/miniprot.gff3").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 9 and fields[2] == "CDS" and "Target=HETE_Q8X1P4 " in fields[8]:
            loci.setdefault(fields[0], set()).add((int(fields[3]), int(fields[4])))
    expected_cds = {(12811, 15093), (15143, 16930)}
    if loci != {"het_e_locus": expected_cds, "het_e_edge": expected_cds}:
        raise RuntimeError(f"Public example did not recover the two expected spliced locus copies: {loci}")
    refinement_counts = {}
    for track, expected in (("priority", {"HETE_Q8X1P4"}),
                            ("comprehensive", {"HETE_Q8X1P4", "HETE_AAL37299", "HETE_A7IQW3"})):
        stage6 = run / "results/stage6"
        if track == "comprehensive":
            stage6 /= "comprehensive"
        refine = (stage6 / "queries/exonerate_refine_ids.work.txt").read_text().split()
        summary = _table(stage6 / "tables/exonerate_summary.tsv")
        if set(refine) != expected or len(refine) != len(expected):
            raise RuntimeError(f"Public example {track} refinement query mismatch: {refine}")
        if {row["protein_id"] for row in summary} != expected or len(summary) != len(expected):
            raise RuntimeError(f"Public example {track} Exonerate mapping mismatch")
        # The duplicated loci tie: either copy may be first in the summary.
        for row in summary:
            if (row["scaffold"] not in {"het_e_locus", "het_e_edge"}
                    or row["strand"] != "+"
                    or not 12811 <= int(row["start"]) < int(row["end"]) <= 16930):
                raise RuntimeError(f"Public example {track} Exonerate locus mismatch: {row}")
        refinement_counts[track] = len(refine)
    return {"candidate_count": len(rows), "strict_candidate_count": len(strict),
            "stage6_refine_count": refinement_counts["priority"],
            "stage6_comprehensive_refine_count": refinement_counts["comprehensive"],
            "spliced_locus_copies": 2, "intron_bp": 49, "tiers": tiers}


def run_public_demo(outdir):
    requested = Path(outdir).expanduser()
    if requested.exists() or requested.is_symlink():
        raise ValueError(f"Demo output must be a new directory: {requested}")
    root = requested.resolve()
    if any(ord(char) < 32 or ord(char) == 127 for char in str(root)):
        raise ValueError("Output paths cannot contain control characters")
    assets, manifest = verified_assets()
    diagnosis = check_installation()
    if not diagnosis["ok"]:
        raise ValueError("Public example requires a complete installation: " + "; ".join(diagnosis["errors"]))
    root.mkdir(parents=True, exist_ok=False)
    logs = root / "logs"
    logs.mkdir()
    evidence_path = root / "evidence.json"
    evidence = {"status": "running", "dataset": "public-sequences", "purpose": PURPOSE,
                "funlr_version": __version__, "installation": diagnosis, "asset_manifest": manifest,
                "evidence_file": str(evidence_path)}
    write_json(evidence_path, evidence)
    print(PURPOSE, flush=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent) + os.pathsep + env.get("PYTHONPATH", "")
    env.update(LC_ALL="C", PYTHONHASHSEED="0")
    try:
        inputs = root / "inputs"
        for name, data in assets.items():
            target = inputs / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        demo_config = root / "demo-config.json"
        write_json(demo_config, DEMO_SETTINGS)
        command = [sys.executable, "-m", "funlr", "run", "--config", str(demo_config), "--genome", str(inputs / "genome.fa"),
                   "--proteins", str(inputs / "proteins.faa"), "--gff3", str(inputs / "annotation.gff3"),
                   "--eggnog", str(inputs / "emapper.annotations"), "--pfam", str(inputs / "db/pfam_mini.hmm"),
                   "--nbd-hmms", str(inputs / "db/nbd.hmm"), "--outdir", str(root / "run"),
                   "--threads", "1", "--sample", "PUBLIC_SEQUENCE_INSTALLATION_TEST"]
        invoke(command, logs / "pipeline.log", env, root, logs / "commands.jsonl")
        state_before = (root / "run/state.json").read_bytes()
        states = json.loads(state_before).get("stages", {})
        if set(states) != {str(n) for n in range(8)} or any(v.get("status") != "complete" for v in states.values()):
            raise RuntimeError("Public example did not complete all eight stages")
        counts = verify_outputs(root / "run")
        commands_before = sha256(root / "run/logs/commands.jsonl")
        invoke([*command, "--resume"], logs / "resume.log", env, root, logs / "commands.jsonl")
        if (root / "run/state.json").read_bytes() != state_before or sha256(root / "run/logs/commands.jsonl") != commands_before:
            raise RuntimeError("Public example resume did not reuse all completed stages unchanged")
        evidence.update(status="passed", completed_stages=8, resume_reused_stages=8, **counts)
        write_json(evidence_path, evidence)
        print(f"Public example passed. Evidence: {evidence_path}", flush=True)
        return evidence
    except BaseException as exc:
        evidence.update(status="failed", error=str(exc) or type(exc).__name__)
        write_json(evidence_path, evidence)
        raise
