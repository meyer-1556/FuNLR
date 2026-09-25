# TSV-driven batch execution

`funlr batch` runs the existing eight-stage annotated-genome workflow once per
sample. It validates the complete sheet before starting, executes samples
sequentially, and records a summary without pooling scientific evidence across
samples. It can run within one scheduler job or on a personal computer with the
normal FuNLR runtime installed.

## Start with the included example

From the source repository root:

```bash
funlr batch --samples examples/cohort/samples.tsv \
  --config examples/cohort/settings.yaml --outdir cohort-demo --validate-only
funlr batch --samples examples/cohort/samples.tsv \
  --config examples/cohort/settings.yaml --outdir cohort-demo --threads 2
```

The example repeats the same public-sequence fixture under `demo_a` and `demo_b`.
It verifies workflow organization; it is not an independent two-sample biological
cohort. Its mini databases and disabled evidence channels are installation-test
settings. See [the example](../examples/cohort/README.md).

For a real batch, supply distinct matching genome/protein/GFF3 inputs and the
production databases described in [DATABASES.md](DATABASES.md). An entire sample's
required paths can come from shared YAML, but every resolved sample must satisfy
the normal [input contract](INPUTS_AND_OUTPUTS.md). Proteome-only analysis, genome
assembly and gene prediction are not implemented modes.

## Sample-sheet contract

Print the supported header, then fill the rows with a text editor or export a spreadsheet as TSV:

```bash
funlr init-samples > samples.tsv
```

Use a UTF-8 text file with one header row, **tabs** between columns and one sample
per row. Header names are case-insensitive and surrounding spaces are removed;
duplicate names are rejected after this normalization. The two identifier
columns are required. `SAMPLE_ID` must contain 1–128 letters, digits, dots, underscores or hyphens,
starting with a letter or digit. IDs must be unique even when compared
case-insensitively; use simple labels such as `isolate_01`.
`SPECIES_ID` can be shared by multiple samples.

| Column | Meaning |
| --- | --- |
| `SAMPLE_ID` | Required unique sample identifier; output directory is `samples/SAMPLE_ID/`. |
| `SPECIES_ID` | Required species/group label recorded in the sample run. It does not select a lineage-specific model bundle. |
| `GENOME_PATH` | Genome nucleotide FASTA; maps to `inputs.genome`. |
| `MASKED_GENOME_PATH` | Alternative genome input used when this row's `GENOME_PATH` is blank. If both are set, `GENOME_PATH` wins and a warning is recorded. This does not perform masking. |
| `PROTEINS_PATH` | Matching amino-acid FASTA; maps to `inputs.proteins`. |
| `GFF3_PATH` | Matching gene annotation; maps to `inputs.gff3`. |
| `EGGNOG_PATH` | Optional eggNOG-mapper annotations; maps to `inputs.eggnog`. |
| `ANNOTATION_TABLE_PATH` | Optional Funannotate table with `GeneID`; maps to `inputs.annotations`. |
| `PFAM_PATH` | Pfam HMM library; maps to `databases.pfam`. |
| `NBD_HMMS` | Fungal NBD library; maps to `databases.nbd_hmms`. |
| `CUSTOM_HMMS` | Optional combined domain library; maps to `databases.custom_hmms`. Use it to include custom domain evidence. |
| `ASM_HMMS` | ASM library; maps to `databases.asm_hmms`. Required while ASM scanning is enabled, as it is by default. |
| `EFFECTOR_HMMS`, `SENSOR_HMMS` | Optional separate source-library paths recorded in provenance. They do not replace the scanned `CUSTOM_HMMS` library. |
| `ASSEMBLY_PROFILE` | `ILLUMINA` or `HIFI`; selects assembly presets before explicit scientific overrides. |
| `NLR_PROFILE` | Alias for `ASSEMBLY_PROFILE`. Either may be used; conflicting nonempty values are rejected. |
| `DISCOVERY_MODE` | `STRICT`, `BALANCED` or `BROAD`; selects the existing discovery preset while retaining explicit YAML thresholds. |
| `ASSEMBLY_VERSION` | Assembly/version label recorded in the sample run. |
| `ORGANISM_KINGDOM` | Optional metadata; when supplied, must be `Funga` (case-insensitive). |
| `ORGANISM_KARYOTE` | Optional metadata; when supplied, must be `eukaryote` (case-insensitive). |
| `PLOIDY` | Optional single-line metadata, such as `haploid` or `unknown`; no ploidy normalization or scientific rule changes. |

The organism metadata check validates a supplied declaration; it does not infer organism identity from sequences. These labels do not select models or change thresholds.

Only supported columns are accepted. This is a FuNLR-specific sheet, conceptually
inspired by [EGAP](https://github.com/iPsychonaut/EGAP); an existing EGAP sheet needs
an explicit conversion. There is no repeat-GFF input or inferred input mode.
`FUNLR_PROFILE`, `REPEATS_PATH`, and `PROTEOME_PATH` are unsupported and rejected;
they do not silently become metadata or select unimplemented analyses. A masked
assembly must match the provided annotation and proteins, just like an ordinary
genome input. `CUSTOM_HMMS` names one existing library file, not a comma-separated
list of files to concatenate during execution.

Path cells may be absolute. Relative paths resolve beside the **TSV file**, not
the current working directory. Quote/escape fields according to ordinary TSV
rules if your editor does so; tabs and newlines are not part of a path or label.
An empty optional cell, `None`, `NA` or `null` means **inherit the shared value**.
It does not explicitly disable a database or erase a shared annotation input.
To omit an optional input from some samples, leave it out of shared YAML and set
it only on the rows that use it.

The genome fallback applies within the TSV row: a nonempty `GENOME_PATH` wins,
otherwise a nonempty `MASKED_GENOME_PATH` overrides the shared genome. If both
cells are absent, the shared genome is inherited. When the primary path wins,
the unused masked path is not validated or scanned. The warning appears in the
sample record of the preflight JSON and the saved batch manifest.

## Shared settings and precedence

Keep scientific thresholds and switches in YAML. A short production starting
point is:

```bash
funlr init-config --minimal --profile ILLUMINA > shared.yaml
```

Edit its database paths and any genuinely shared inputs, then place each sample's
own inputs in the TSV. `execution.output_dir` in shared YAML is superseded by the
batch's required `--outdir`; sample destinations are always `ROOT/samples/ID/`.
The complete `init-config` output remains available when you want to explicitly
freeze all defaults, but those written thresholds still override a different
profile selected in a row.

Values resolve in this order:

1. FuNLR defaults and the selected assembly preset.
2. Explicit settings from shared YAML, if supplied.
3. Nonempty supported TSV cells.
4. Explicit batch `--threads` and `--ram-gb` resource overrides.

A TSV `ASSEMBLY_PROFILE` (or `NLR_PROFILE`) and `DISCOVERY_MODE` select existing
presets while retaining explicit scientific thresholds from YAML. Their values
are case-insensitive. With minimal YAML, omitted thresholds can therefore follow
each row's presets. For example, `HIFI` with `STRICT` uses the HIFI assembly
defaults and the strict discovery default, but a YAML `thresholds.eval_pfam` or
`discovery.pfam_enrich_min_log2fc` still takes precedence for that threshold.
The effective settings for every sample are recorded, including
version-dependent defaults. YAML paths resolve relative to the YAML file before
row values are applied; TSV paths resolve relative to the TSV file. Do not move
either file without accounting for its relative paths.

Batch execution controls `resume` and `dry_run` through its own flags. Stored
sample configurations describe the analysis with these operational switches set
to false; they do not silently request recovery or planning on their own.

## Execution, preflight and failures

```bash
funlr batch --samples samples.tsv --config shared.yaml \
  --outdir analysis-batch --threads 8
```

`--config` is optional when the sheet and defaults fully specify valid inputs.
`--outdir` is required. Equivalent sample-sheet options are `--input-tsv` and
`--input_tsv`; use `--samples` for new documentation and scripts.

- `--validate-only` checks every sample's input/configuration contract without
  creating a batch or running tools. It is not a check of every external program;
  run `funlr doctor` in the intended environment as well.
- `--dry-run` validates and reports the batch plan without creating outputs or
  running tools.
- Samples execute in sheet order, one at a time. `--threads` is the per-sample
  tool allocation, not the number of concurrent samples. `--ram-gb` records an
  allocation label; the scheduler/operating system enforces memory limits.
- The default stops after a sample failure and records it. `--keep-going` attempts
  the remaining prevalidated samples, while retaining failed status and a failing
  overall exit status. It does not skip invalid-sheet preflight checks.
  `--stop-on-error` (or `--stop_on_error`) explicitly selects the default;
  combining it with `--keep-going` is an error.

A valid empty candidate set is a completed scientific result. A failed or
unattempted sample has missing counts, not zero candidates. Consult the sample's
manifest and logs before interpreting a blank batch-summary cell.

## Batch output layout

All paths below are relative to the batch `--outdir`:

| Path | Contents |
| --- | --- |
| `samples.tsv` | Exact input-sheet snapshot. |
| `configs/SAMPLE_ID.yaml` | Fully resolved configuration for that sample, with absolute input paths and isolated output directory. |
| `samples/SAMPLE_ID/` | The normal FuNLR run directory: eight stages, final results, manifest, settings, state, logs and working copies. |
| `batch_manifest.json` | Fixed batch identity, source-sheet/configuration provenance and resolved sample records. |
| `batch_state.json` | Recorded progress for each sample. |
| `batch_summary.json`, `batch_summary.tsv` | Per-sample metadata, status, paths and counts from completed scientific results. |

Summary counts distinguish `n_candidates` (original candidates), `n_report_rows`
(including fusion review rows), `n_strict_candidates` and `n_fusion_models`.
Incomplete samples have JSON `null` / empty TSV count cells. These summaries do
not establish equivalent assembly quality, model calibration or biological
independence across rows; use comparable inputs and settings before making
cross-sample interpretations.

## Resume and inspect

Reuse the same sheet, effective settings, inputs, software and resource choices:

```bash
funlr batch --samples samples.tsv --config shared.yaml \
  --outdir analysis-batch --threads 8 --resume
```

Batch resume verifies the fixed sheet and effective sample configurations. Every
child run re-enters normal FuNLR resume, including output-integrity verification;
a prior batch success label alone does not make a sample reusable. Damaged stages
and their dependents are recomputed. Changed scientific settings, inputs or
software require a new output root. Do not edit the generated sheet or configs as
a way to mutate an existing batch.

Use `status` on the batch root or an individual sample. Batch verification checks
the recorded stage/final outputs of all manifest-listed child runs; it does not
revalidate the batch-summary tables or original inputs. HTML reporting accepts
one completed **sample run**:

```bash
funlr status --run-dir analysis-batch --verify --json
funlr status --run-dir analysis-batch/samples/isolate_01 --verify --json
funlr report --run-dir analysis-batch/samples/isolate_01 \
  --output isolate-01-report.html
```

The batch summary supplies the cohort overview. The HTML command verifies the
recorded output checksums before writing the explicitly requested snapshot; it does not
rerun stages or change the batch evidence. Retain both batch and sample manifests
for provenance and recovery. See [REPORTS.md](REPORTS.md) for inspection details
and HTML output rules.
