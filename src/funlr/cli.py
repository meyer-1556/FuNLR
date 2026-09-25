"""CLI for single runs, sample sheets, and read-only result inspection."""
import argparse
from contextlib import contextmanager
import json
import signal
import sys

from . import __version__
from .runner import STAGES, run
from .validation import validate_inputs


@contextmanager
def termination_handling():
    """Let subprocess cleanup run when a scheduler terminates the CLI."""
    def terminate(signum, frame):
        raise KeyboardInterrupt("Terminated")
    previous_handler = signal.signal(signal.SIGTERM, terminate)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous_handler)


def parser():
    root = argparse.ArgumentParser(prog="funlr", description="Portable fungal NLR discovery and evidence review")
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    commands.add_parser("version", help="Print version")
    commands.add_parser("init-samples", help="Print a TSV sample-sheet header")
    from .batch import add_batch_arguments
    add_batch_arguments(commands.add_parser("batch", help="Validate and run a TSV cohort sequentially"))
    status = commands.add_parser("status", help="Read recorded run or cohort status without changing results")
    status.add_argument("--run-dir", required=True)
    status.add_argument("--verify", action="store_true", help="Also verify recorded output checksums")
    status.add_argument("--json", action="store_true", help="Print structured JSON")
    report = commands.add_parser("report", help="Create a self-contained HTML report from completed results")
    report.add_argument("--run-dir", required=True)
    report.add_argument("--output", required=True, help="New .html file outside pipeline-owned subdirectories")
    init = commands.add_parser("init-config", help="Print a complete starter YAML configuration")
    init.add_argument("--profile", choices=("ILLUMINA", "HIFI"), default="ILLUMINA")
    init.add_argument("--minimal", action="store_true", help="Print only inputs, databases, execution and assembly profile")
    qc = commands.add_parser("input-qc", help="Summarize input masking, ambiguity and annotation structure without external tools")
    qc.add_argument("--config")
    for field in ("genome", "proteins", "gff3"):
        qc.add_argument("--" + field)
    qc.add_argument("--output", help="Write JSON to this file instead of stdout")
    inventory = commands.add_parser("model-inventory", help="Inspect HMM profiles, checksums and known provenance without changing files")
    inventory.add_argument("--directory", required=True)
    inventory.add_argument("--verify-snapshot", action="store_true")
    benchmark = commands.add_parser("benchmark", help="Compare an explicit called-ID set with a labeled reference panel")
    benchmark.add_argument("--reference", required=True, help="Reference JSON manifest")
    benchmark.add_argument("--predictions", required=True, help="One-column TSV containing called IDs")
    benchmark.add_argument("--input", required=True, help="FASTA defining the evaluation units and checksum")
    benchmark.add_argument("--id-column", default="protein_id")
    verify = commands.add_parser("verify-models", help="Verify all five fungal HMM libraries against the reference snapshot")
    verify.add_argument("--directory", required=True)
    fetch = commands.add_parser("fetch-pfam", help="Download and verify Pfam38.2 against the reference checksum")
    fetch.add_argument("--outdir", "--dest", dest="outdir", required=True)
    doctor = commands.add_parser("doctor", help="Check the Python environment and external tools without running an analysis")
    doctor.add_argument("--plot-backend", choices=("matplotlib", "r"))
    doctor.add_argument("--config", help="Use configured tool paths and plot backend; scientific inputs are not checked")
    doctor.add_argument("--self-test", action="store_true", help="Also verify packaged fixture checksums and demo dependencies")
    doctor.add_argument("--models-dir", help="With --self-test, verify the five production HMMs in this directory")
    doctor.add_argument("--demo", action="store_true", help="Also check hmmbuild, used by the installation demo")
    demo = commands.add_parser("demo", help="Run a bundled installation example with real tools")
    demo.add_argument("--outdir", required=True, help="New directory for example inputs and test results")
    demo.add_argument("--dataset", choices=("synthetic", "public-sequences"), default="synthetic",
                      help="Synthetic positive-rescue test (default) or public-sequence hybrid example")
    for command in ("run", "validate"):
        sub = commands.add_parser(command, help="Run stages 0–7" if command == "run" else "Validate file/identifier consistency")
        for flag in ("genome", "proteins", "gff3", "pfam", "eggnog", "annotations"):
            sub.add_argument("--" + flag)
        sub.add_argument("--nbd-hmms", "--nbd_hmms", dest="nbd_hmms")
        sub.add_argument("--custom-hmms", "--custom_hmms", dest="custom_hmms")
        sub.add_argument("--asm-hmms", dest="asm_hmms")
        sub.add_argument("--effector-hmms", dest="effector_hmms")
        sub.add_argument("--sensor-hmms", dest="sensor_hmms")
        sub.add_argument("--profile", "--nlr-profile", "--nlr_profile", dest="profile", choices=("ILLUMINA", "HIFI"))
        sub.add_argument("--discovery-mode", "--discovery_mode", dest="discovery_mode", choices=("STRICT", "BALANCED", "BROAD"))
        sub.add_argument("--config", help="Nested YAML/JSON or legacy flat JSON settings")
        sub.add_argument("-o", "--outdir", "--output-dir", "--output_dir", dest="outdir")
        sub.add_argument("-t", "--threads", "--cpu-threads", "--cpu_threads", dest="threads", type=int)
        sub.add_argument("-r", "--ram-gb", "--ram_gb", dest="ram_gb", type=int,
                         help="Memory allocation label for provenance; does not enforce a limit")
        sub.add_argument("--sample", "--sample-id", "--sample_id", dest="sample")
        sub.add_argument("--species-id", "--species_id", dest="species_id")
        sub.add_argument("--assembly-version", dest="assembly_version")
        if command == "run":
            sub.add_argument("--resume", action="store_true", default=None, help="Verify and reuse stages from the identical run")
            sub.add_argument("--dry-run", "--dry_run", dest="dry_run", action="store_true", default=None,
                             help="Validate and print PLANNED stages; write nothing and require no tools")
        else:
            sub.add_argument("--self-test", action="store_true", help="Installation checks instead of scientific-input validation (alias of doctor --self-test)")
            sub.add_argument("--models-dir", help="With --self-test, verify the five production HMMs in this directory")
    return root


def resolve_args(args):
    from .config import FunlrConfig
    cfg = FunlrConfig.from_yaml(args.config) if args.config else FunlrConfig()
    cfg.apply_overrides({
        "inputs": {key: getattr(args, key) for key in ("genome", "proteins", "gff3", "eggnog", "annotations")},
        "databases": {key: getattr(args, key) for key in ("pfam", "nbd_hmms", "custom_hmms", "asm_hmms", "effector_hmms", "sensor_hmms")},
        "context": {"profile": args.profile},
        "discovery": {"discovery_mode": getattr(args, "discovery_mode", None)},
        "execution": {"output_dir": args.outdir, "cpu_threads": args.threads, "ram_gb": args.ram_gb,
                      "resume": getattr(args, "resume", None), "dry_run": getattr(args, "dry_run", None)},
        "sample": {"sample_id": args.sample, "species_id": args.species_id, "assembly_version": args.assembly_version},
    })
    for section in ("inputs", "databases"):
        for key, value in cfg[section].items():
            setattr(args, key, value)
    for attr, section, key in (("outdir", "execution", "output_dir"), ("threads", "execution", "cpu_threads"),
                               ("ram_gb", "execution", "ram_gb"), ("resume", "execution", "resume"),
                               ("dry_run", "execution", "dry_run"), ("sample", "sample", "sample_id"),
                               ("species_id", "sample", "species_id"), ("assembly_version", "sample", "assembly_version")):
        setattr(args, attr, cfg.get(section, key))
    args.sample = args.sample or "sample"
    args.assembly_version = args.assembly_version or "unspecified"
    args.tools = cfg["tools"]
    args.scientific = cfg.scientific_settings()
    return cfg.scientific_settings()


def installation_report(args):
    """Keep both installation entry points on one diagnostic implementation."""
    from .config import FunlrConfig
    from .doctor import check_installation, check_self_test
    if args.models_dir is not None and not args.self_test:
        raise ValueError("--models-dir requires --self-test")
    if args.command == "validate":
        ignored = ("genome", "proteins", "gff3", "eggnog", "annotations", "pfam",
                   "nbd_hmms", "custom_hmms", "asm_hmms", "effector_hmms", "sensor_hmms",
                   "profile", "discovery_mode", "outdir", "threads", "ram_gb", "sample",
                   "species_id", "assembly_version")
        supplied = [name for name in ignored if getattr(args, name, None) is not None]
        if supplied:
            raise ValueError("--self-test checks installation; use plain validate for input/run options: " + ", ".join(supplied))
    cfg = FunlrConfig.from_yaml(args.config) if args.config else FunlrConfig()
    backend = getattr(args, "plot_backend", None) or cfg.scientific_settings()["REPORT_PLOT_BACKEND"]
    common = {"configured_tools": cfg["tools"], "plot_backend": backend}
    if args.self_test:
        return check_self_test(**common, models_dir=args.models_dir)
    return check_installation(**common, include_demo=args.demo)


def main(argv=None):
    args = parser().parse_args(argv)
    if args.command == "version":
        print(__version__)
        return 0
    if args.command == "init-config":
        from .config import write_example_yaml
        print(write_example_yaml(profile=args.profile, minimal=args.minimal), end="")
        return 0
    if args.command == "init-samples":
        from .batch import sample_template
        print(sample_template(), end="")
        return 0
    try:
        if args.command == "input-qc":
            from .config import FunlrConfig
            from .input_qc import collect_input_qc, write_input_qc
            cfg = FunlrConfig.from_yaml(args.config) if args.config else FunlrConfig()
            cfg.apply_overrides({"inputs": {name: getattr(args, name) for name in ("genome", "proteins", "gff3")}})
            paths = {name: cfg["inputs"].get(name) for name in ("genome", "proteins", "gff3")}
            missing = [name for name, path in paths.items() if not path]
            if missing:
                raise ValueError("Missing input QC files: " + ", ".join(missing))
            result = collect_input_qc(**paths)
            if args.output:
                write_input_qc(result, args.output)
            else:
                print(json.dumps(result, indent=2))
            return 0
        if args.command == "model-inventory":
            from .model_inventory import inventory_models
            print(json.dumps(inventory_models(args.directory, verify_snapshot=args.verify_snapshot), indent=2))
            return 0
        if args.command == "benchmark":
            from .benchmark import compare_reference
            print(json.dumps(compare_reference(args.reference, args.predictions, args.input, id_column=args.id_column), indent=2))
            return 0
        if args.command == "batch":
            from .batch import run_batch
            with termination_handling():
                result = run_batch(args)
            print(json.dumps(result, indent=2))
            return 1 if result.get("status") == "FAILED" else 0
        if args.command == "status":
            from .report import inspect_run, format_status
            result = inspect_run(args.run_dir, verify=args.verify)
            print(json.dumps(result, indent=2) if args.json else format_status(result))
            return 1 if args.verify and result.get("integrity", {}).get("status") != "VERIFIED" else 0
        if args.command == "report":
            from .report import write_html_report
            print(json.dumps(write_html_report(args.run_dir, args.output), indent=2))
            return 0
        if args.command == "verify-models":
            from .databases import verify_models
            print(json.dumps(verify_models(args.directory), indent=2))
            return 0
        if args.command == "fetch-pfam":
            from .databases import fetch_pfam
            print(json.dumps(fetch_pfam(args.outdir), indent=2))
            return 0
        if args.command == "doctor" or (args.command == "validate" and args.self_test):
            report = installation_report(args)
            print(json.dumps(report, indent=2))
            return 0 if report["ok"] else 1
        if args.command == "demo":
            if args.dataset == "public-sequences":
                from .public_demo import run_public_demo as run_demo
            else:
                from .demo import run_demo
            with termination_handling():
                print(json.dumps(run_demo(args.outdir), indent=2))
            return 0
        if args.command == "validate" and args.models_dir is not None:
            raise ValueError("--models-dir requires --self-test")
        config = resolve_args(args)
        inputs, warnings, counts = validate_inputs(args)
        for warning in warnings:
            print("WARNING: " + warning, file=sys.stderr)
        if args.command == "validate":
            print(json.dumps({"valid": True, "counts": counts, "warnings": warnings}, indent=2))
            return 0
        if not args.outdir:
            raise ValueError("Missing --outdir (or execution.output_dir in configuration)")
        if args.dry_run:
            print(json.dumps({"status": "PLANNED", "stages": list(STAGES), "inputs": {k: str(v) for k, v in inputs.items()},
                              "config": config, "outdir": args.outdir, "threads": args.threads}, indent=2))
            return 0
        with termination_handling():
            run(args, inputs, config, warnings, counts)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        message = ("Interrupted; partial demo evidence is in the output directory" if args.command == "demo"
                   else "Interrupted; completed stages can be verified with --resume")
        print(message, file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
