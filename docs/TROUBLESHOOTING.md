# Troubleshooting FuNLR

This guide describes the commands and recovery rules in FuNLR 0.5.0a1.
Start with the failing check or stage, preserve its logs, and use a new output
directory when changing the analysis or software. For installation and platform
limits, see [INSTALLATION.md](INSTALLATION.md).

## Choose the check that answers the problem

Run diagnostics in the same environment and compute allocation that will run
FuNLR. An HPC login shell can have a different PATH from its batch jobs.

| Command | What it checks |
| --- | --- |
| `funlr doctor --config settings.yaml` | Supported Python/dependency versions and executable tool probes, using the configured tool paths and plot backend. Scientific inputs are not required. |
| `funlr doctor --self-test --config settings.yaml` | The installation checks, `hmmbuild` for the demos, and checksums of the packaged public-sequence fixture. No searches or plot rendering are performed. |
| `funlr doctor --self-test --models-dir /path/to/models` | Also checks all five production fungal HMM libraries against the fixed reference snapshot. Add `--config settings.yaml` when using configured tool paths. |
| `funlr validate --config settings.yaml` | Required input/database files, FASTA IDs, GFF3 structure and coordinates, and protein-to-annotation mappings. Requires no external tools or output directory. |
| `funlr input-qc --config settings.yaml` | Input syntax and descriptive genome, protein and GFF3 observations. Prints JSON by default; optional `--output` writes a new file. Needs no models/tools and creates no run directory. |
| `funlr run --config settings.yaml --dry-run` | Input validation followed by a `PLANNED` stage/configuration report. Creates no run files and executes no tools. |
| `funlr status --run-dir analysis --verify --json` | Recorded completion and checksums of recorded stage/final outputs. Does not recheck original input bytes or execute tools. |

`funlr validate --self-test` is an alias for `funlr doctor --self-test`. Both
accept `--config` and `--models-dir`; use plain `validate` for scientific-input
options such as `--genome`. The `--models-dir` flag requires `--self-test`.
The diagnostics return JSON and a nonzero exit status if a required check fails.
Missing, nonexecutable or unsuccessfully probed tools fail the check; a filename
alone does not establish a working installation.

When `--models-dir` is omitted, `production_models.status` is `NOT_REQUESTED`.
That means no production-model verification was requested. The fixture check
still runs, and a missing or checksum-mismatched fixture fails on both source
and wheel installations. The small public fixture is included in the installed
package, with its checksum manifest and data notices.

For actual searches, alignments and rendering, run a demo in a **new** directory:

```bash
funlr demo --outdir synthetic-demo
funlr demo --dataset public-sequences --outdir public-demo
funlr status --run-dir public-demo/run --verify --json
```

These are installation exercises. Their small synthetic or hybrid model sets
are not production fungal databases or a biological benchmark. The public demo
combines public sequences with documented synthetic fixture material; see its
[data description](../src/funlr/data/public_sequences/README.md).

## Resume an interrupted run

Keep the same inputs, effective settings, input paths, sample labels, resource
choices, Python runtime, FuNLR code and external tools, then add `--resume`:

```bash
funlr run --config settings.yaml --outdir analysis --threads 8 --resume
```

FuNLR compares the current run fingerprint with the saved one before reuse.
It checks recorded outputs before skipping a completed stage. If a stage's
outputs are missing or changed, that stage and its dependents are recomputed;
their old products and final reports are removed before rebuilding. Preserve
the run directory, including `state.json`, `manifest.json`, `work/` and logs,
rather than retaining only the final TSVs.

| Symptom | Next step |
| --- | --- |
| `Output directory is not empty` | Use `--resume` for the identical run, or select a new directory for a new analysis. |
| `Cannot resume: inputs, parameters, runtime, tools, or code changed` | Compare the saved `manifest.json` and `resolved_config.json` with the intended run. Use a new output directory after an intentional change or upgrade. |
| A stage was recomputed on resume | Inspect `state.json`, the recorded output checksums and `logs/stageN.log`. Modified or missing evidence requires rebuilding its dependents. |
| `Run lock exists` | Check the scheduler and original job before acting. After confirming that job and its worker processes have stopped, remove only the stale `.run.lock`, then resume. A recorded PID is not a liveness check. |
| A batch stopped on one sample | Inspect `batch_state.json` and that sample's logs. Resume the same sheet and settings with `funlr batch ... --resume`; each child run receives the normal integrity checks. |

There is no arbitrary `--force_stage` or standalone `--stage` command in this
release. Do not edit saved state to force reuse or splice stages from different
analyses. A dry run also leaves an existing run's markers and products alone.
See [BATCH.md](BATCH.md) for batch failure/recovery flags and
[REPORTS.md](REPORTS.md) for read-only inspection.

## Inputs and settings

**A genome FASTA, matching protein FASTA and GFF3 are required.** A proteome-only
workflow is not implemented. EggNOG-mapper and Funannotate annotation tables are
optional inputs; supplying or omitting evidence can change results. The fungal
NBD and Pfam libraries are required, and the ASM library is required while its
default evidence channel is enabled. Include the custom domain library to
collect custom domain evidence from the intended model collection. See the full
[input contract](INPUTS_AND_OUTPUTS.md#input-contract).

| Symptom | What to inspect |
| --- | --- |
| Protein/GFF3 mapping failure | Use files from the same annotation. Check FASTA first-token IDs, transcript IDs, CDS `protein_id`/`orig_protein_id` attributes and gene parents. After Stage 0, inspect `results/stage0/id_mapping_report.tsv` and `id_coords.tsv`. |
| Invalid eggNOG table | Preserve its named `#query` header. A nonempty malformed table fails at Stage 0; initial file validation is not every stage parser's schema check. |
| Missing ASM library | Supply `databases.asm_hmms`, or deliberately set `asm.asm_enable: 0` in a new run's configuration. Disabling it changes the evidence collected. |
| A relative path points to the wrong file | YAML paths resolve beside the YAML file. Batch row paths resolve beside the TSV file. Moving those files can change what a relative path means. |
| Changing profile did not change a threshold | Explicit YAML thresholds override profile/discovery presets. Inspect the resolved settings rather than inferring them from the profile label. |

Generate a complete configuration for the intended profile before editing paths:

```bash
funlr init-config --profile HIFI > hifi-settings.yaml
```

A complete ILLUMINA configuration already contains explicit thresholds; changing
only its profile label later preserves those thresholds. For a mixed-profile
batch, start with `funlr init-config --minimal > shared.yaml` and put each row's
profile in `ASSEMBLY_PROFILE` or its `NLR_PROFILE` alias. Conflicting aliases fail
validation. A row's profile and `DISCOVERY_MODE` select existing presets, while
explicit shared YAML thresholds retain precedence. Species and ploidy metadata
do not infer an assembly profile or select different HMMs. See
[batch precedence](BATCH.md#shared-settings-and-precedence) and
[scientific settings](METHODS.md).

## Databases, indexing and disk space

The public code package does not include the five production HMM libraries.
MIT licensing of FuNLR code does not settle redistribution rights for those
models. Their public permissions, immutable distribution locations and component
build provenance remain separate release work; the supplied private snapshot
can be checked with:

```bash
funlr verify-models --directory /path/to/models
funlr fetch-pfam --outdir /path/to/databases/pfam38.2
```

A model verification failure means the files do not match the fixed snapshot
by the recorded checks. Check the exact filename and restore/retransfer the
intended unchanged file. Rebuilt models need their own documented identity and
validation; changing the expected checksum would not establish equivalence.
The Pfam downloader separately checks the decompressed Pfam 38.2 file against
the recorded checksum. See [DATABASES.md](DATABASES.md) and [MODELS.md](MODELS.md).

Unpressed source HMMs are accepted. FuNLR creates private database copies and
pressed indexes under the run's `work/` directory, checks reusable index
evidence, and copies the genome there for indexing. It does not need write
permission beside the source genome or shared HMMs. Allow disk space for these
copies and indexes as well as alignments and final outputs. Pfam alone is about
2.25 GB uncompressed, before its extra working copy and indexes.

For resource failures, inspect the scheduler's termination reason together with
the stage log. `--threads` controls tool allocation within one sample;
`--ram-gb` records an allocation label and does not enforce a memory limit.
Batch samples run sequentially. Follow the site's scratch/compute policies and
retain the complete run when moving results. See [HPC.md](HPC.md).

## Empty results, classifications and review

Successful searches can produce empty candidate or rescue sets. Expected tables
retain headers and empty FASTAs remain empty. A tool failure is recorded as a
failure; it must not be interpreted as zero candidates. The two explicit
fallbacks are failed Pfam NBD extraction and a failed optional whole-proteome
Pfam background scan. They emit explicit warnings and follow the documented
fallback behavior; other checked tool failures stop execution. Inspect logs
before comparing evidence counts across runs.

The classification contains expanded Tier 1–4 subcategories and fusion-review
rows. It is not a three-tier scheme or a promise that every possible tier file
will appear. Inspect `results/stage3/tables/tiered_candidates.tsv`, the export
manifests and `final_results/nlr_final_report.tsv`. Tier-specific exports exist
for categories present in the run. An empty priority-rescue set does not prevent
the independent comprehensive rescue track from running.

Fusion rows are evidence-gated reconstruction hypotheses, marked for review;
they are not validated new genes. Strict candidates remain computational
predictions. An empty strict set does not establish biological absence, and a
predicted NLR architecture does not establish immune function. Use the recorded
domain, alignment and context evidence, with the acceptance and interpretation
rules in [RESCUE.md](RESCUE.md).

For an inspectable snapshot of a completed run:

```bash
funlr report --run-dir analysis --output analysis-report.html
```

The HTML command verifies recorded output checksums before writing. Choose a
new `.html` filename with an existing parent directory, outside protected
pipeline subdirectories. An HTML output cannot overwrite an existing file.
The full TSV remains the source for review rows beyond the report's displayed
limit. See [REPORTS.md](REPORTS.md).

## Information to include in a bug report

Include the FuNLR version, installation method/platform, exact command, effective
configuration, `doctor` output, relevant `logs/stageN.log` and the scheduler's
error or termination message. `logs/commands.jsonl` records executed commands;
`manifest.json` records input, tool and code identities. For batch issues, add
the sample ID and relevant batch state. A small reproducible input helps when
the failure depends on a FASTA/GFF3/table format. Review paths and biological
data before sharing the diagnostic files.
