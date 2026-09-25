"""Opt-in installed-package cohort check using the public hybrid fixture.

Two sample IDs deliberately use the same data. This checks isolation, path
resolution, metadata, resume, and reporting; it is not a biological cohort.
Run with the real toolchain active: python tests/real_tools_batch.py --outdir NEW
"""
import argparse
import csv
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml

from funlr import __version__
from funlr.public_demo import DEMO_SETTINGS, verified_assets, verify_outputs
from funlr.report import inspect_run


class ReportTables(HTMLParser):
    """Read displayed cells, so a successful HTML write is not enough to pass."""

    def __init__(self):
        super().__init__()
        self.tables = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.tables.append([])
        elif tag == "tr":
            self.tables[-1].append([])
        elif tag in ("td", "th"):
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ("td", "th"):
            self.tables[-1][-1].append("".join(self.cell))
            self.cell = None


def verify_html_rows(run, html_path):
    parsed = ReportTables()
    parsed.feed(html_path.read_text())
    review = next(table for table in parsed.tables if table[0][0] == "Protein ID")
    keys = ("protein_id", "tier", "rescue_priority", "flags", "nbd_confidence", "fusion_members")
    expected = [[str(row.get(key, "")) for key in keys]
                for row in inspect_run(run)["review_rows"]]
    assert len(expected) > 1, "Fixture must exercise distinct review rows"
    assert len({row[0] for row in expected}) > 1
    assert review[1:] == expected, "Displayed HTML review rows differ from recorded data"
    return len(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root = Path(args.outdir).resolve()
    root.mkdir(parents=True, exist_ok=False)
    assets, asset_manifest = verified_assets()
    for name, content in assets.items():
        path = root / "input files" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    (root / "settings").mkdir()
    settings = {**DEMO_SETTINGS,
                "databases": {"nbd_hmms": "../input files/db/nbd.hmm",
                              "pfam": "../input files/db/pfam_mini.hmm"},
                "inputs": {"eggnog": "../input files/emapper.annotations"}}
    config = root / "settings/common.yaml"
    config.write_text(yaml.safe_dump(settings))
    (root / "sheets").mkdir()
    sheet = root / "sheets/samples.tsv"
    with sheet.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample_ID", "Species_ID", "GENOME_PATH", "PROTEINS_PATH",
                         "GFF3_PATH", "ORGANISM_KINGDOM", "PLOIDY", "MASKED_GENOME_PATH",
                         "NLR_PROFILE", "DISCOVERY_MODE"])
        for sample in ("fixture_A", "fixture_B"):
            writer.writerow([sample, "Podospora_anserina", "" if sample == "fixture_A" else "../input files/genome.fa",
                             "../input files/proteins.faa", "../input files/annotation.gff3",
                             "Funga", "1", "../input files/genome.fa" if sample == "fixture_A" else "unused-masked.fa",
                             "ILLUMINA", "BROAD"])
    env = os.environ.copy()
    env.update(LC_ALL="C", PYTHONHASHSEED="0")
    prefix = [sys.executable, "-m", "funlr"]
    cohort = root / "cohort"
    command = [*prefix, "batch", "--samples", str(sheet), "--config", str(config),
               "--outdir", str(cohort), "--threads", "1"]

    def execute(argv, name):
        result = subprocess.run(argv, cwd=root, env=env, capture_output=True, text=True,
                                timeout=600)
        (root / name).write_text(result.stdout + "\n" + result.stderr)
        if result.returncode:
            raise RuntimeError(f"Command failed ({result.returncode}); see {root / name}")
        return result

    diagnosis = execute([*prefix, "doctor", "--self-test", "--config", str(config)], "self-test.log")
    health = json.loads(diagnosis.stdout)
    assert health["ok"] and health["public_fixture"]["status"] == "VERIFIED"
    assert health["production_models"]["status"] == "NOT_REQUESTED"
    execute([*command, "--dry-run"], "dry-run.log")
    assert not cohort.exists(), "Dry-run created the cohort directory"
    execute(command, "batch.log")
    outputs = {}
    states = {}
    html_rows = {}
    for sample in ("fixture_A", "fixture_B"):
        run = cohort / "samples" / sample
        outputs[sample] = verify_outputs(run)
        states[sample] = (run / "state.json").read_bytes()
        manifest = json.loads((run / "manifest.json").read_text())
        assert manifest["sample"] == sample
        assert manifest["species_id"] == "Podospora_anserina"
        assert manifest["status"] == "complete"
        assert len(json.loads(states[sample])["stages"]) == 8
        execute([*prefix, "report", "--run-dir", str(run), "--output", str(root / f"{sample}.html")],
                f"{sample}-report.log")
        html_rows[sample] = verify_html_rows(run, root / f"{sample}.html")
    execute([*command, "--resume"], "resume.log")
    for sample, state in states.items():
        assert (cohort / "samples" / sample / "state.json").read_bytes() == state
    status_result = execute([*prefix, "status", "--run-dir", str(cohort), "--verify", "--json"],
                            "status.log")
    status = json.loads(status_result.stdout)
    assert status["integrity"]["status"] == "VERIFIED"
    manifest = json.loads((cohort / "batch_manifest.json").read_text())
    assert len(manifest["samples"]) == 2
    assert all(row["metadata"]["PLOIDY"] == "1" for row in manifest["samples"])
    evidence = {"status": "passed", "funlr_version": __version__,
                "purpose": "Two copies of a constructed public fixture; workflow validation, not biological benchmarking.",
                "samples": outputs, "stages_completed_per_sample": 8,
                "resume_preserved_both_stage_states": True, "cohort_integrity": "VERIFIED",
                "dry_run_wrote_nothing": True, "html_reports": 2,
                "html_review_rows_matched": html_rows,
                "installed_self_test": {"ok": health["ok"], "public_fixture": health["public_fixture"]["status"],
                                        "production_models": health["production_models"]["status"]},
                "tsv_compatibility": ["mixed-case headers", "masked genome fallback", "primary genome precedence",
                                      "NLR_PROFILE alias", "per-sample DISCOVERY_MODE"],
                "sample_sheet_sha256": hashlib.sha256(sheet.read_bytes()).hexdigest(),
                "asset_manifest": asset_manifest}
    (root / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
