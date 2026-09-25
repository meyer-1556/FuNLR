"""Pipeline contract checks using synthetic input and explicit executable doubles.

These check orchestration, parsing, and the preserved tier rules, not sensitivity,
specificity, scientific equivalence to an reference run, or real tool compatibility.
Run from the project root: PYTHONPATH=src python -m unittest discover -s tests -v
"""

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ("samtools", "hmmpress", "hmmscan", "seqkit", "miniprot", "exonerate", "hmmfetch", "gffread", "Rscript")
EXPECTED_TIERS = {
    "canonical": "TIER_1A_HIGH_CONFIDENCE",
    "review": "TIER_1B_NEEDS_REVIEW",
    "rescue_edge": "TIER_2A_HIGH_PRIORITY_RESCUE",
    "rescue_internal": "TIER_2B_RESCUE_CANDIDATE",
    "variant": "TIER_3B_ARCHITECTURAL_VARIANT",
    "repeat_only": "TIER_4A_REPEAT_ONLY_NO_NBD",
    "housekeeping": "TIER_4B_LIKELY_STAND_HOUSEKEEPING",
    "priority_only": "TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL",
}


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_ids(path):
    return {line[1:].split()[0] for line in Path(path).read_text().splitlines() if line.startswith(">")}


class PipelineTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="funlr synthetic ")
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.inputs = self.work / "input files with apostrophe's"
        self.inputs.mkdir()
        self.bin = self.work / "fixture bin"
        self.bin.mkdir()
        tool_source = (ROOT / "tests" / "fixture_tool.py").read_text().split("\n", 1)[1]
        for name in TOOLS:
            path = self.bin / name
            path.write_text(f"#!{sys.executable}\n" + tool_source)
            path.chmod(0o755)
        self.log = self.work / "tool-invocations.jsonl"
        self.env = dict(os.environ)
        self.env.update({
            "PATH": str(self.bin) + os.pathsep + str(Path(sys.executable).parent) + os.pathsep + os.environ.get("PATH", ""),
            "PYTHONPATH": str(ROOT / "src"),
            "FUNLR_TEST_DOUBLES": "1",
            "FUNLR_TEST_TOOL_LOG": str(self.log),
        })
        self.out = self.work / "results with apostrophe's"

    def make_inputs(self, scenario="positive", eggnog=True):
        if scenario == "positive":
            proteins = [*EXPECTED_TIERS, "bad_eval", "short_hit"]
        elif scenario == "no_rescue":
            proteins = ["canonical"]
        elif scenario == "zero":
            proteins = ["bad_eval", "short_hit", "no_hit"]
        else:
            raise AssertionError(scenario)
        with (self.inputs / "genome.fa").open("w") as genome, (self.inputs / "proteins.fa").open("w") as fasta, (self.inputs / "genes.gff3").open("w") as gff:
            gff.write("##gff-version 3\n")
            for pid in proteins:
                genome.write(f">chr_{pid}\n" + ("ACGT" * 15000) + "\n")
                fasta.write(f">{pid} synthetic fixture protein\n" + ("M" + "ACDEFGHIKLMNPQRSTVWY" * 42 + "A") + "\n")
                start = 501 if pid in {"review", "rescue_edge"} else 20001
                gff.write(f"chr_{pid}\tSYNTHETIC\tgene\t{start}\t{start + 2399}\t.\t+\t.\tID=gene_{pid}\n")
                gff.write(f"chr_{pid}\tSYNTHETIC\tmRNA\t{start}\t{start + 2399}\t.\t+\t.\tID={pid};Parent=gene_{pid}\n")
                gff.write(f"chr_{pid}\tSYNTHETIC\tCDS\t{start}\t{start + 2399}\t.\t+\t0\tID=cds_{pid};Parent={pid}\n")
        for name in ("nbd.hmm", "pfam.hmm", "custom.hmm"):
            (self.inputs / name).write_text("HMMER3/f [SYNTHETIC TEST PLACEHOLDER — NOT A REAL HMM]\n//\n")
        descriptions = {
            "repeat_only": "WD40 repeat domain",
            "priority_only": "NOD-like hypothetical protein",
            "housekeeping": "CDC6 replication factor",
        }
        if eggnog:
            with (self.inputs / "eggnog.annotations").open("w") as out:
                out.write("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
                for pid in proteins:
                    columns = [pid, "seed", "0", "0", "-", "-", "-", descriptions.get(pid, "-"), *(["-"] * 13)]
                    out.write("\t".join(columns) + "\n")
        (self.inputs / "annotations.txt").write_text("GeneID\tProduct\n" + "".join(f"{pid}\t{descriptions.get(pid, 'hypothetical protein')}\n" for pid in proteins))
        self.original_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.inputs.iterdir()}
        args = [
            "--genome", str(self.inputs / "genome.fa"),
            "--proteins", str(self.inputs / "proteins.fa"),
            "--gff3", str(self.inputs / "genes.gff3"),
            "--annotations", str(self.inputs / "annotations.txt"),
            "--nbd-hmms", str(self.inputs / "nbd.hmm"),
            "--pfam", str(self.inputs / "pfam.hmm"),
            "--outdir", str(self.out), "--threads", "1",
        ]
        if eggnog:
            args.extend(["--eggnog", str(self.inputs / "eggnog.annotations")])
        return args

    def configured_args(self, args):
        # This fixture intentionally exercises every legacy tier via BROAD.
        # Separate tests check BALANCED's narrower candidate union. No ASM model
        # is supplied, so disable that channel explicitly; keep exports/plots on
        # to exercise native tool interfaces through their doubles. This suite
        # explicitly tests the optional R backend; test_plots covers matplotlib.
        settings = {"DISCOVERY_MODE": "BROAD", "ASM_ENABLE": 0, "ENABLE_FUSION": 0, "REPORT_PLOT_BACKEND": "r"}
        remaining = []
        index = 0
        while index < len(args):
            if args[index] == "--config":
                settings.update(json.loads(Path(args[index + 1]).read_text()))
                index += 2
            else:
                remaining.append(args[index])
                index += 1
        config = self.work / "resolved-test-config.json"
        config.write_text(json.dumps(settings, sort_keys=True))
        return [*remaining, "--config", str(config)]

    def run_pipeline(self, args, success=True):
        args = self.configured_args(args)
        result = subprocess.run([sys.executable, "-m", "funlr", "run", *args], cwd=ROOT, env=self.env, text=True, capture_output=True, timeout=120)
        diagnostics = result.stdout + "\n" + result.stderr
        if success and result.returncode:
            logs = sorted((self.out / "logs").glob("stage*.log"), key=lambda p: p.stat().st_mtime)
            if logs:
                diagnostics += "\n" + logs[-1].read_text()[-6000:]
        if success:
            self.assertEqual(result.returncode, 0, diagnostics)
        else:
            self.assertNotEqual(result.returncode, 0, diagnostics)
        return result

    def calls(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()] if self.log.exists() else []

    def assert_inputs_untouched(self):
        self.assertEqual({p.name for p in self.inputs.iterdir()}, set(self.original_hashes), "Indexes must be created under the output directory, not alongside user inputs")
        for name, expected in self.original_hashes.items():
            self.assertEqual(hashlib.sha256((self.inputs / name).read_bytes()).hexdigest(), expected)

    def test_complete_run_exercises_eight_tiers_and_both_rescue_tracks(self):
        args = self.make_inputs()
        args.extend(["--custom-hmms", str(self.inputs / "custom.hmm")])
        self.run_pipeline(args)
        report = read_tsv(self.out / "final_results" / "nlr_final_report.tsv")
        self.assertEqual({row["protein_id"]: row["tier"] for row in report}, EXPECTED_TIERS)
        rows = {row["protein_id"]: row for row in report}
        self.assertEqual(rows["canonical"]["flags"], "PASS")
        self.assertEqual(rows["rescue_edge"]["rescue_priority"], "100")
        self.assertEqual(rows["rescue_internal"]["rescue_priority"], "85")
        self.assertIn("NON_NLR_ANNOT", rows["housekeeping"]["flags"])
        self.assertEqual(rows["review"]["mp_scaffold"], "chr_review")
        self.assertEqual(rows["rescue_internal"]["exo_scaffold"], "chr_rescue_internal")
        strict = read_tsv(self.out / "final_results" / "nlr_strict_candidates.tsv")
        self.assertEqual({row["protein_id"] for row in strict}, {"canonical", "review", "rescue_edge", "rescue_internal", "variant"})
        stage1 = self.out / "results" / "stage1"
        self.assertEqual(set((stage1 / "nbd_candidate_ids.txt").read_text().splitlines()), {"canonical", "review", "rescue_edge", "rescue_internal", "variant"})
        self.assertEqual(fasta_ids(self.out / "results" / "stage4" / "tiered_candidates.faa"), set(EXPECTED_TIERS))
        with (self.out / "results" / "stage3" / "beds" / "nlr_candidates.bed").open() as handle:
            bed = [line.rstrip("\n").split("\t") for line in handle]
        canonical = next(row for row in bed if row[3].startswith("canonical|"))
        self.assertEqual(canonical[1:3], ["20000", "22400"], "GFF 1-based closed to BED 0-based half-open")
        self.assertEqual(set(call["tool"] for call in self.calls()), set(TOOLS))
        plots = self.out / 'final_results' / 'domain_plots'
        self.assertEqual((plots / 'nlr_domains_by_tier.png').read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
        self.assertEqual({p.stem for p in (plots / 'individual').glob('*.png')}, {"canonical", "review", "rescue_edge", "rescue_internal", "variant"})
        self.assert_inputs_untouched()

    def test_no_priority_rescue_candidates_still_runs_comprehensive_mapping(self):
        self.run_pipeline(self.make_inputs("no_rescue", eggnog=False))
        report = read_tsv(self.out / "final_results" / "nlr_final_report.tsv")
        self.assertEqual([(r["protein_id"], r["tier"]) for r in report], [("canonical", EXPECTED_TIERS["canonical"])])
        alignments = [c for c in self.calls() if c["tool"] in {"miniprot", "exonerate"}]
        self.assertEqual([c["tool"] for c in alignments], ["miniprot"])
        self.assertIn("comprehensive", alignments[0]["args"][-1])
        self.assertEqual(fasta_ids(alignments[0]["args"][-1]), {"canonical"})
        self.assert_inputs_untouched()

    def test_zero_candidates_completes_with_header_only_reports(self):
        self.run_pipeline(self.make_inputs("zero", eggnog=False))
        for name in ("nlr_final_report.tsv", "nlr_strict_candidates.tsv"):
            path = self.out / "final_results" / name
            self.assertEqual(read_tsv(path), [])
            self.assertIn("protein_id", path.read_text().splitlines()[0])
        self.assertTrue(all(c["tool"] not in {"miniprot", "exonerate"} for c in self.calls()))
        self.assert_inputs_untouched()

    def test_no_miniprot_loci_sends_every_rescue_query_to_exonerate(self):
        self.env["FUNLR_TEST_MINIPROT_NO_HITS"] = "1"
        self.run_pipeline(self.make_inputs())
        refine = self.out / "results" / "stage5" / "queries" / "exonerate_refine_ids.txt"
        self.assertEqual(set(refine.read_text().splitlines()), {"review", "rescue_edge", "rescue_internal"})
        report = {r["protein_id"]: r for r in read_tsv(self.out / "final_results" / "nlr_final_report.tsv")}
        for pid in ("review", "rescue_edge", "rescue_internal"):
            self.assertEqual(report[pid]["exo_scaffold"], "chr_" + pid)

    def test_hmmscan_failure_is_not_reported_as_success(self):
        self.env["FUNLR_TEST_FAIL_TOOL"] = "hmmscan"
        self.run_pipeline(self.make_inputs(), success=False)
        self.assertFalse((self.out / "final_results" / "nlr_final_report.tsv").exists())
        self.assertTrue(any(c["tool"] == "hmmscan" for c in self.calls()))
        self.assertTrue(all(c["tool"] not in {"miniprot", "exonerate"} for c in self.calls()))
        self.assertEqual(json.loads((self.out / "manifest.json").read_text())["status"], "failed")
        state = json.loads((self.out / "state.json").read_text())
        self.assertNotEqual(state["stages"]["1"]["status"], "complete")
        self.assertFalse((self.out / ".run.lock").exists())

    def test_termination_interrupts_worker_and_resume_recovers(self):
        args = self.make_inputs("no_rescue")
        self.env["FUNLR_TEST_WAIT_TOOL"] = "hmmscan"
        process = subprocess.Popen(
            [sys.executable, "-m", "funlr", "run", *self.configured_args(args)],
            cwd=ROOT, env=self.env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        try:
            deadline = time.monotonic() + 20
            while not any(call["tool"] == "hmmscan" for call in self.calls()):
                if process.poll() is not None or time.monotonic() > deadline:
                    self.fail("The isolated HMMER double did not enter its interruptible wait")
                time.sleep(0.05)
            process.terminate()
            stdout, stderr = process.communicate(timeout=10)
            self.assertEqual(process.returncode, 130, stdout + stderr)
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate(timeout=10)
        self.assertIn("Interrupted", stderr)
        tool_pid = next(call["pid"] for call in self.calls() if call["tool"] == "hmmscan")
        deadline = time.monotonic() + 5
        while True:
            try:
                os.kill(tool_pid, 0)
            except ProcessLookupError:
                break
            if time.monotonic() > deadline:
                self.fail("Terminating the CLI left its external HMMER process alive")
            time.sleep(0.05)
        manifest = json.loads((self.out / "manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["error"], "Terminated")
        self.assertFalse((self.out / ".run.lock").exists())
        self.assertFalse((self.out / "final_results" / "nlr_final_report.tsv").exists())
        self.assertEqual(sum(c["tool"] == "samtools" for c in self.calls()), 1)
        del self.env["FUNLR_TEST_WAIT_TOOL"]
        self.run_pipeline([*args, "--resume"])
        self.assertEqual(sum(c["tool"] == "samtools" for c in self.calls()), 1)
        self.assertEqual(json.loads((self.out / "manifest.json").read_text())["status"], "complete")
        self.assert_inputs_untouched()

    def test_exonerate_failure_is_not_silently_ignored(self):
        self.env["FUNLR_TEST_FAIL_TOOL"] = "exonerate"
        self.run_pipeline(self.make_inputs(), success=False)
        self.assertFalse((self.out / "final_results" / "nlr_final_report.tsv").exists())
        self.assertTrue(any(c["tool"] == "exonerate" for c in self.calls()))

    def test_config_threshold_reaches_candidate_parser(self):
        args = self.make_inputs()
        config = self.work / "stricter threshold.json"
        config.write_text(json.dumps({"MIN_NBD_ALI_LEN": 181}))
        self.run_pipeline([*args, "--config", str(config)])
        self.assertEqual((self.out / "results" / "stage1" / "nbd_candidate_ids.txt").read_text(), "")
        # Raising the strict threshold must reach the parser without disabling
        # the independent relaxed NBD and Pfam discovery channels.
        stage1 = self.out / "results" / "stage1"
        self.assertEqual(set((stage1 / "nbd_candidate_ids_relaxed.txt").read_text().splitlines()), {"canonical", "review", "rescue_edge", "rescue_internal", "variant", "housekeeping"})
        rows = {r["protein_id"]: r for r in read_tsv(self.out / "results" / "stage2" / "architecture_summary.tsv")}
        self.assertEqual(rows["canonical"]["has_nbd_stage1"], "False")
        self.assertEqual(rows["canonical"]["has_nbd_pfam"], "True")

    def test_balanced_union_omits_keyword_only_candidates(self):
        args = self.make_inputs()
        config = self.work / "balanced.json"
        config.write_text(json.dumps({"DISCOVERY_MODE": "BALANCED"}))
        self.run_pipeline([*args, "--config", str(config)])
        ids = {row["protein_id"] for row in read_tsv(self.out / "final_results" / "nlr_final_report.tsv")}
        self.assertEqual(ids, set(EXPECTED_TIERS) - {"repeat_only", "priority_only"})

    def test_input_identifier_mismatch_fails_before_tools_run(self):
        args = self.make_inputs("no_rescue")
        gff = self.inputs / "genes.gff3"
        gff.write_text(gff.read_text().replace("ID=canonical;", "ID=different_id;"))
        self.run_pipeline(args, success=False)
        self.assertEqual(self.calls(), [])

    def test_existing_output_is_not_overwritten(self):
        args = self.make_inputs("no_rescue")
        self.out.mkdir()
        sentinel = self.out / "existing-work.txt"
        sentinel.write_text("Keep my existing result.\n")
        self.run_pipeline(args, success=False)
        self.assertEqual(sentinel.read_text(), "Keep my existing result.\n")
        self.assertEqual(self.calls(), [])

    def test_dry_run_needs_no_external_tools_and_writes_nothing(self):
        args = self.make_inputs()
        self.env["PATH"] = ""
        result = self.run_pipeline([*args, "--dry-run"])
        plan = json.loads(result.stdout)
        self.assertEqual(len(plan["stages"]), 8)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])
        self.assert_inputs_untouched()

    def test_resume_verifies_outputs_and_skips_completed_stages(self):
        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        calls_before = self.calls()
        report = self.out / "final_results" / "nlr_final_report.tsv"
        before = report.read_bytes()
        self.run_pipeline([*args, "--resume"])
        self.assertEqual(self.calls(), calls_before)
        self.assertEqual(report.read_bytes(), before)
        state = json.loads((self.out / "state.json").read_text())
        self.assertEqual({k: v["status"] for k, v in state["stages"].items()}, {str(n): "complete" for n in range(8)})
        self.assertEqual(json.loads((self.out / "manifest.json").read_text())["status"], "complete")

    def test_resume_rejects_changed_configuration(self):
        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        calls_before = self.calls()
        config = self.work / "changed.json"
        config.write_text(json.dumps({"MIN_NBD_ALI_LEN": 151}))
        self.run_pipeline([*args, "--resume", "--config", str(config)], success=False)
        self.assertEqual(self.calls(), calls_before)

    def test_new_evidence_is_recorded_and_repaired_by_resume(self):
        from funlr.report import inspect_run

        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        state = json.loads((self.out / "state.json").read_text())
        self.assertIn("results/stage0/input_qc.json", state["stages"]["0"]["outputs"])
        evidence = self.out / "final_results/evidence/candidate_evidence.tsv"
        self.assertIn("final_results/evidence/candidate_evidence.tsv", state["stages"]["7"]["outputs"])
        original = evidence.read_bytes()
        evidence.write_text("CORRUPTED\n")
        self.assertNotEqual(inspect_run(self.out, verify=True)["integrity"]["status"], "VERIFIED")
        result = self.run_pipeline([*args, "--resume"])
        self.assertEqual(result.stdout.count("Verified; reusing"), 7)
        self.assertEqual(evidence.read_bytes(), original)
        self.assertEqual(inspect_run(self.out, verify=True)["integrity"]["status"], "VERIFIED")

    def test_damaged_output_is_recomputed_without_rerunning_previous_stages(self):
        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        report = self.out / "final_results" / "nlr_final_report.tsv"
        original_report = report.read_bytes()
        architecture = self.out / "results" / "stage2" / "architecture_summary.tsv"
        architecture.write_text("CORRUPTED\n")
        self.run_pipeline([*args, "--resume"])
        self.assertEqual(report.read_bytes(), original_report)
        scans = [c for c in self.calls() if c["tool"] == "hmmscan"]
        self.assertEqual(sum("nbd_whole.domtblout" in " ".join(c["args"]) for c in scans), 1)
        self.assertEqual(sum("pfam.domtblout" in " ".join(c["args"]) for c in scans), 2)

    def test_failed_resume_removes_stale_final_report_and_can_recover(self):
        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        report = self.out / "final_results" / "nlr_final_report.tsv"
        original_report = report.read_bytes()
        (self.out / "results" / "stage1" / "nbd_whole.domtblout").write_text("CORRUPTED\n")
        self.env["FUNLR_TEST_FAIL_TOOL"] = "hmmscan"
        self.run_pipeline([*args, "--resume"], success=False)
        self.assertFalse(report.exists())
        self.assertEqual(json.loads((self.out / "manifest.json").read_text())["status"], "failed")
        state = json.loads((self.out / "state.json").read_text())
        self.assertNotIn("7", state["stages"])
        del self.env["FUNLR_TEST_FAIL_TOOL"]
        self.run_pipeline([*args, "--resume"])
        self.assertEqual(report.read_bytes(), original_report)

    def test_zero_domain_threshold_spellings_keep_original_pfam_fallback(self):
        from funlr.config import read_config

        for index, value in enumerate((0.0, "0e0")):
            with self.subTest(value=value):
                self.out = self.work / f"zero-sentinel-{index}"
                args = self.make_inputs("no_rescue")
                config = self.work / f"zero-{index}.json"
                config.write_text(json.dumps({"MAX_DOM_I_EVAL": value}))
                self.assertEqual(read_config(config)["MAX_DOM_I_EVAL"], "0")
                self.run_pipeline([*args, "--config", str(config)])
                report = read_tsv(self.out / "final_results" / "nlr_final_report.tsv")
                self.assertEqual([(r["protein_id"], r["tier"]) for r in report], [("canonical", EXPECTED_TIERS["canonical"])])

    def test_invalid_regex_is_a_clean_configuration_error(self):
        from funlr.config import read_config

        args = self.make_inputs("no_rescue")
        config = self.work / "invalid-regex.json"
        config.write_text(json.dumps({"NON_NLR_DOMAIN_REGEX": "["}))
        with self.assertRaisesRegex(ValueError, "NON_NLR_DOMAIN_REGEX"):
            read_config(config)
        result = self.run_pipeline([*args, "--config", str(config)], success=False)
        self.assertIn("NON_NLR_DOMAIN_REGEX", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertEqual(self.calls(), [])
        self.assertFalse(self.out.exists())

    def test_resume_rejects_symlinked_top_level_results_and_work_before_tools(self):
        for name in ("results", "work"):
            with self.subTest(name=name):
                self.out = self.work / f"symlink-check-{name}"
                args = self.make_inputs("no_rescue")
                self.out.mkdir()
                target = self.work / f"harmless-target-{name}"
                target.mkdir()
                sentinel = target / "sentinel.txt"
                sentinel.write_text("This temporary target must remain untouched.\n")
                (self.out / name).symlink_to(target, target_is_directory=True)
                result = self.run_pipeline([*args, "--resume"], success=False)
                self.assertIn("Symlinks", result.stderr)
                self.assertEqual(self.calls(), [], "Reject before invoking external tools")
                self.assertEqual(list(target.iterdir()), [sentinel])
                self.assertEqual(sentinel.read_text(), "This temporary target must remain untouched.\n")
                self.assertFalse((self.out / "manifest.json").exists())

    def test_stale_fasta_index_is_regenerated_when_stage_zero_is_damaged(self):
        args = self.make_inputs("no_rescue")
        self.run_pipeline(args)
        report = self.out / "final_results" / "nlr_final_report.tsv"
        original_report = report.read_bytes()
        genome_index = self.out / "work" / "genome.fa.fai"
        original_index = genome_index.read_bytes()
        genome_index.write_text("chr_canonical\t1000\t0\t60\t61\n")
        (self.out / "results" / "stage0" / "master_table.tsv").write_text("CORRUPTED\n")
        self.run_pipeline([*args, "--resume"])
        self.assertEqual(genome_index.read_bytes(), original_index)
        self.assertEqual(report.read_bytes(), original_report)
        self.assertEqual(sum(c["tool"] == "samtools" for c in self.calls()), 2)


if __name__ == "__main__":
    unittest.main()
