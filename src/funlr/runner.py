"""Sequential orchestration; original analysis is retained in bundled stages."""
import hashlib
import csv
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

STAGES_DIR = Path(__file__).resolve().parent / "stages"
STAGES = (
    "stage0_validate", "stage1_discover", "stage2_architecture", "stage3_tier",
    "stage4_export", "stage5_miniprot", "stage6_exonerate", "stage7_finalize",
)
TOOLS = {"samtools": ["--version"], "hmmpress": ["-h"],
         "hmmscan": ["-h"], "seqkit": ["version"], "miniprot": ["--version"],
         "exonerate": ["--version"], "hmmfetch": ["-h"],
         "gffread": ["--version"], "Rscript": ["--version"]}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def r_package_inventory(executable):
    """Record the packages and library paths actually used for domain plots."""
    script = """pkgs <- c('dplyr','ggplot2','RColorBrewer','colorspace','tidyr','scales')
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop(paste('Missing R packages:', paste(missing, collapse=', ')))
cat('LIBRARIES\\n'); cat(.libPaths(), sep='\\n'); cat('\\nPACKAGES\\n')
for (p in pkgs) cat(p, as.character(packageVersion(p)), sep='\\t', fill=TRUE)
"""
    result = subprocess.run([executable, '--vanilla', '-e', script], capture_output=True,
                            text=True, errors='replace', timeout=60)
    if result.returncode:
        raise ValueError('R plot dependency check failed: ' + result.stderr[-2000:])
    return result.stdout.strip()


def tool_inventory(configured=None, scientific=None):
    inventory = {}
    missing = []
    scientific = scientific or {}
    for name, flags in TOOLS.items():
        if scientific is not None:
            if name == 'Rscript' and not (scientific.get('REPORT_PLOTS', 1) and scientific.get('REPORT_PLOT_BACKEND', 'matplotlib') == 'r'):
                continue
            if name == 'gffread' and not (scientific.get('REPORT_INTEGRATION', 1) and scientific.get('REPORT_PROTEINS', 1) and scientific.get('REPORT_TAGGED_GFF3', 1)):
                continue
        executable = shutil.which((configured or {}).get(name, name))
        if not executable:
            missing.append(name)
            continue
        try:
            result = subprocess.run([executable, *flags], stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, errors="replace", timeout=15)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise ValueError(f"Cannot run {name}: {exc}") from exc
        # Exonerate 2.4.0 prints its version and normally exits 1.
        normal_exonerate_version = (name == "exonerate" and result.returncode == 1
                                    and "exonerate" in result.stdout.lower() and "2.4.0" in result.stdout)
        if result.returncode and not normal_exonerate_version:
            raise ValueError(f"Version check failed for {name} (exit {result.returncode}): {result.stdout[:400]}")
        inventory[name] = {"path": str(Path(executable).absolute()), "sha256": sha256(executable),
                           "version_output": result.stdout[:12000].strip()}
        if name == 'Rscript':
            inventory[name]['r_packages'] = r_package_inventory(executable)
    if missing:
        raise ValueError("Missing executables on PATH: " + ", ".join(missing) + ". Install environment.yml and activate it.")
    return inventory


def invoke(command, logfile, env, cwd, command_log):
    entry = {"started": now(), "argv": list(map(str, command)), "cwd": str(cwd), "log": str(logfile)}
    with command_log.open("a") as out:
        out.write(json.dumps(entry) + "\n")
    with logfile.open("w") as log:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            returncode = process.wait()
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            raise
    entry.update(finished=now(), returncode=returncode)
    with command_log.open("a") as out:
        out.write(json.dumps(entry) + "\n")
    if returncode:
        raise RuntimeError(f"{Path(str(command[0])).name} exited {returncode}; see {logfile}")


def file_hashes(root, folder):
    return {str(p.relative_to(root)): sha256(p) for p in sorted(folder.rglob("*")) if p.is_file()}


def hashes_match(root, hashes):
    return bool(hashes) and all((root / name).is_file() and sha256(root / name) == digest
                               for name, digest in hashes.items())


def validate_stage_outputs(root, number):
    """Completion requires the stage's public handoff, including empty schemas."""
    tables = {
        0: ("master_table.tsv", {"protein_id", "scaffold", "start", "end", "protein_length"}),
        1: ("nbd_hits.tsv", {"protein_id", "dom_i_evalue", "ali_len"}),
        2: ("architecture_summary.tsv", {"protein_id", "has_nbd", "has_sensor", "order_invalid", "nbd_confidence"}),
        3: ("tables/tiered_candidates.tsv", {"protein_id", "tier", "flags"}),
        4: ("fasta_export_summary.tsv", {"file", "bytes", "seqs"}),
        7: ("tables/nlr_final_report.tsv", {"protein_id", "tier", "flags"}),
    }
    stage = root / "results" / f"stage{number}"
    if number in tables:
        name, expected = tables[number]
        path = stage / name
        if not path.is_file():
            raise RuntimeError(f"Stage {number} missing required table: {path}")
        with path.open(newline="") as handle:
            columns = next(csv.reader(handle, delimiter="\t"), [])
        if not expected.issubset(columns):
            raise RuntimeError(f"Stage {number} wrote an invalid table header: {path}")
    required = {
        0: ["proteins_clean.faa", "contig_lengths.tsv", "input_qc.json"],
        1: ["nbd_candidate_ids.txt", "union_candidate_ids.txt", "union_candidates.faa"],
        4: ["tiered_candidates.faa", "rescue_priority.faa"],
        5: ["queries/rescue_queries.faa", "queries/exonerate_refine_ids.txt", "tables/miniprot_summary.tsv", "gff/miniprot.gff3"],
        6: ["tables/exonerate_summary.tsv", "gff/exonerate_merged.gff3"],
    }
    for name in required.get(number, []):
        if not (stage / name).is_file():
            raise RuntimeError(f"Stage {number} missing required output: {stage / name}")
    if number == 0:
        qc = json.loads((stage / "input_qc.json").read_text())
        if not isinstance(qc, dict) or qc.get("schema_version") != 1 or not {"genome", "proteins", "gff3"}.issubset(qc):
            raise RuntimeError("Stage 0 wrote an invalid input-QC report")
    if number == 7:
        for name in ("nlr_final_report.tsv", "nlr_strict_candidates.tsv", "tier_summary.tsv", "qc_summary.tsv", "nlr_candidates.bed", "nlr_miniprot_loci.bed", "nlr_exonerate_loci.bed"):
            if not (root / "final_results" / name).is_file():
                raise RuntimeError(f"Missing final report: {name}")
        for parent in (stage, root / "final_results"):
            evidence = parent / "evidence"
            for name in ("candidate_evidence.tsv", "review_and_exclusions.tsv", "nbd_hit_evidence.tsv", "nbd_aligned_segments.faa", "evidence_summary.json"):
                if not (evidence / name).is_file():
                    raise RuntimeError(f"Missing final evidence output: {evidence / name}")


def run(args, inputs, config, warnings, counts):
    import numpy
    import pandas
    import yaml
    from importlib.metadata import version as package_version

    root = Path(args.outdir).expanduser().resolve()
    if any(c in str(root) for c in "\n\r\x00"):
        raise ValueError("Output paths cannot contain newlines or NUL characters")
    if root.is_file():
        raise ValueError(f"Output path is a file: {root}")
    for name in ("results", "work", "logs", "final_results", "state.json", "manifest.json", "resolved_config.json"):
        reserved = root / name
        candidates = [reserved, *reserved.rglob("*")] if reserved.is_dir() and not reserved.is_symlink() else [reserved]
        if any(path.is_symlink() for path in candidates):
            raise ValueError(f"Symlinks in pipeline-owned output paths are unsupported: {reserved}")
    for path in inputs.values():
        if path.is_relative_to(root):
            raise ValueError("Inputs must be outside the output directory")
    tools = tool_inventory(getattr(args, "tools", None), config)
    source_files = sorted(p for p in Path(__file__).parent.rglob("*") if p.is_file() and p.suffix in (".py", ".R", ".json"))
    fingerprint = {
        "version": __version__, "config": config, "threads": args.threads,
        "sample": args.sample, "assembly_version": args.assembly_version,
        "species_id": getattr(args, "species_id", None), "ram_gb_label": getattr(args, "ram_gb", 8),
        "runtime": {"python": platform.python_version(), "pandas": pandas.__version__, "numpy": numpy.__version__, "pyyaml": yaml.__version__, "matplotlib": package_version("matplotlib")},
        "inputs": {key: {"path": str(path), "sha256": sha256(path)} for key, path in inputs.items()},
        "tools": tools, "code": {str(p.relative_to(Path(__file__).parent)): sha256(p) for p in source_files},
    }
    signature = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
    state_path, manifest_path = root / "state.json", root / "manifest.json"
    state = {"signature": signature, "stages": {}}
    if root.exists() and any(root.iterdir()):
        if not args.resume:
            raise ValueError(f"Output directory is not empty: {root}; use --resume for the same run, or a new directory")
        if not state_path.is_file():
            raise ValueError("Cannot resume: output directory has no FuNLR state.json")
        state = json.loads(state_path.read_text())
        if not isinstance(state, dict) or not isinstance(state.get("stages"), dict):
            raise ValueError("Cannot resume: invalid state.json structure")
        if state.get("signature") != signature:
            raise ValueError("Cannot resume: inputs, parameters, runtime, tools, or code changed. Use a new output directory.")
    root.mkdir(parents=True, exist_ok=True)
    lock = root / ".run.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise ValueError(f"Run lock exists: {lock}. If its recorded process has stopped, remove this stale lock before resuming.") from None
    with os.fdopen(descriptor, "w") as handle:
        handle.write(json.dumps({"pid": os.getpid(), "started": now()}) + "\n")
    manifest = {**fingerprint, "signature": signature, "started": now(), "status": "running",
                "warnings": warnings, "input_counts": counts, "upstream_commit": None,
                "analysis_policy": "fungal-nlr-2026.04; see the versioned methods and rescue documentation"}
    try:
        for subdir in ("logs", "results", "work"):
            (root / subdir).mkdir(exist_ok=True)
        write_json(manifest_path, manifest)
        write_json(state_path, state)
        write_json(root / "resolved_config.json", config)
        env = os.environ.copy()
        for key in ("BASH_ENV", "ENV", "CDPATH"):
            env.pop(key, None)
        # The worker must import this exact installation even when the caller
        # used a relative PYTHONPATH and the stage changes working directory.
        env["PYTHONPATH"] = str(Path(__file__).resolve().parent.parent) + os.pathsep + env.get("PYTHONPATH", "")
        env.update({key: str(value) for key, value in config.items()})
        env.update({key: str(value) for key, value in inputs.items()})
        env.update(EGGNOG_FILE=str(inputs.get("EGGNOG_FILE", "")), CUSTOM_HMMS=str(inputs.get("CUSTOM_HMMS", "")),
                   NLR_DIR=str(root), BASE_DIR=str(root), RESULTS_DIR=str(root / "results"),
                   LOG_DIR=str(root / "logs"), SCRIPTS_DIR=str(STAGES_DIR), FUNLR_PYTHON=sys.executable,
                   THREADS=str(args.threads), GENOME_NAME=args.sample, ASSEMBLY_VERSION=args.assembly_version,
                   LC_ALL="C", PYTHONHASHSEED="0")
        command_log = root / "logs" / "commands.jsonl"

        # Work copies keep samtools/hmmpress away from the user's inputs and shared databases.
        genome = root / "work" / "genome.fa"
        if not genome.exists() or sha256(genome) != fingerprint["inputs"]["GENOME_FILE"]["sha256"]:
            shutil.copyfile(inputs["GENOME_FILE"], genome)
            genome.with_suffix(".fa.fai").unlink(missing_ok=True)
        env["GENOME_FILE"] = str(genome)
        for key in ("NBD_HMMS", "PFAM_DB", "CUSTOM_HMMS", "ASM_HMMS", "COMBINED_EFFECTOR_HMMS", "COMBINED_SENSOR_HMMS"):
            if key not in inputs:
                continue
            target = root / "work" / (key.lower() + ".hmm")
            digest = fingerprint["inputs"][key]["sha256"]
            if not target.exists() or sha256(target) != digest:
                shutil.copyfile(inputs[key], target)
                for ext in ("h3f", "h3i", "h3m", "h3p"):
                    Path(str(target) + "." + ext).unlink(missing_ok=True)
            indexes = [Path(str(target) + "." + ext) for ext in ("h3f", "h3i", "h3m", "h3p")]
            recorded = state.get("indexes", {}).get(key, {})
            if not all(p.is_file() and p.stat().st_size > 0 for p in indexes) or not hashes_match(root, recorded):
                invoke([tools["hmmpress"]["path"], "-f", str(target)], root / "logs" / (key.lower() + "_hmmpress.log"), env, root, command_log)
                if not all(p.is_file() and p.stat().st_size > 0 for p in indexes):
                    raise RuntimeError(f"hmmpress did not create all indexes for {key}")
                state.setdefault("indexes", {})[key] = {str(p.relative_to(root)): sha256(p) for p in indexes}
                write_json(state_path, state)
            env[key] = str(target)

        # Stages receive only resolved settings. The worker runs in its own
        # process group, so interruption stops its external tools as well.
        from .config import FunlrConfig
        cfg = FunlrConfig(config)
        cfg.apply_overrides({
            "inputs": {"genome": env["GENOME_FILE"], "proteins": env["PROTEINS_FILE"],
                       "gff3": env["GFF_FILE"], "eggnog": env["EGGNOG_FILE"] or None,
                       "annotations": env.get("ANNOTATIONS_FILE") or None},
            "databases": {"pfam": env["PFAM_DB"], "nbd_hmms": env["NBD_HMMS"],
                          "custom_hmms": env["CUSTOM_HMMS"] or None,
                          "asm_hmms": env.get("ASM_HMMS") or None,
                          "effector_hmms": env.get("COMBINED_EFFECTOR_HMMS") or None,
                          "sensor_hmms": env.get("COMBINED_SENSOR_HMMS") or None},
            "sample": {"sample_id": args.sample, "species_id": getattr(args, "species_id", None),
                       "assembly_version": args.assembly_version},
            "execution": {"output_dir": str(root), "cpu_threads": args.threads,
                          "ram_gb": getattr(args, "ram_gb", 8)},
            "tools": {name: item["path"] for name, item in tools.items()},
        })
        runtime_config = root / "work" / "runtime.json"
        write_json(runtime_config, {"config": cfg.to_dict(), "inventory": tools})

        invalidate = False
        for number, script in enumerate(STAGES):
            previous = state["stages"].get(str(number), {})
            stage_path = root / "results" / f"stage{number}"
            valid = previous.get("status") == "complete" and hashes_match(root, previous.get("outputs", {}))
            if args.resume and not invalidate and valid:
                print(f"[{number}/7] Verified; reusing {script}", flush=True)
                continue
            if not invalidate:
                # Invalidate all downstream metadata AND products, so a failed rerun
                # cannot leave old final reports looking current.
                for downstream in range(number, len(STAGES)):
                    state["stages"].pop(str(downstream), None)
                    old = root / "results" / f"stage{downstream}"
                    if old.is_symlink():
                        raise ValueError(f"Refusing to remove symlinked stage directory: {old}")
                    if old.exists():
                        shutil.rmtree(old)
                final = root / "final_results"
                if final.is_symlink():
                    raise ValueError(f"Refusing to remove symlinked final directory: {final}")
                if final.exists():
                    shutil.rmtree(final)
                invalidate = True
            state["stages"][str(number)] = {"status": "running", "started": now()}
            write_json(state_path, state)
            print(f"[{number}/7] Running {script}", flush=True)
            logfile = root / "logs" / f"stage{number}.log"
            if number == 0:
                # Stage 0 must never trust an old FASTA index during recovery.
                Path(env["GENOME_FILE"] + ".fai").unlink(missing_ok=True)
            invoke([sys.executable, "-m", "funlr.stage_worker", str(number), str(runtime_config)], logfile, env, root, command_log)
            validate_stage_outputs(root, number)
            outputs = file_hashes(root, stage_path)
            if number == 7:
                outputs.update(file_hashes(root, root / "final_results"))
            if not outputs:
                raise RuntimeError(f"Stage {number} returned success but wrote no output files; see {logfile}")
            state["stages"][str(number)].update(status="complete", finished=now(), outputs=outputs)
            write_json(state_path, state)
        manifest.update(status="complete", finished=now(), final_outputs=file_hashes(root, root / "final_results"))
        write_json(manifest_path, manifest)
        print(f"Completed. Final reports: {root / 'final_results'}", flush=True)
    except BaseException as exc:
        manifest.update(status="failed", finished=now(), error=str(exc) or type(exc).__name__)
        write_json(manifest_path, manifest)
        raise
    finally:
        lock.unlink(missing_ok=True)
