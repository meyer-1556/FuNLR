# Inspect runs and export an HTML report

These commands read existing run evidence. They do not execute the scientific
pipeline, infer process liveness or change candidate classifications.

## Inspect an individual run or a batch

```bash
funlr status --run-dir /path/to/run
funlr status --run-dir /path/to/run --verify --json
funlr status --run-dir /path/to/batch --verify --json
```

Point at the directory containing `manifest.json` for an individual run, or
`batch_manifest.json` for a batch. The one-command installation demos put the
analysis in an extra `run/` directory: after
`funlr demo --dataset public-sequences --outdir public-demo`, inspect
`public-demo/run/`. A direct `funlr run --outdir analysis` uses `analysis/` itself.

Without `--verify`, status reads metadata and available final tables. It reports
recorded stage completion, candidate counts, warnings, lock-file presence and
consistency issues. **Integrity is `NOT_CHECKED`.** A recorded completion label is
not proof that output bytes are intact. A PID or lock file is also not proof
that an HPC job is still running; consult the scheduler for process status.

`--verify` additionally hashes each unique recorded stage/final output and
compares it with the saved SHA-256. Verification is `VERIFIED` only when recorded
completion and output checks are consistent. Changed, missing, unsafe or
unreadable outputs and inconsistent records are reported. For batches, this
checks the manifest-listed child runs; it does not independently revalidate the
batch summary tables. Original inputs, executable files and manifest authenticity
are not revalidated by this command. Normal pipeline resume performs its own
broader input/configuration/tool checks before reuse.

`--json` emits structured JSON with `schema_version: 1`. The `reporter_version`
identifies the inspecting FuNLR version; the recorded analysis version remains
in provenance. Count fields distinguish original candidates (`n_candidates`),
all final-report rows including fusion hypotheses (`n_report_rows`), strict
candidates and fusion models. Missing results remain unknown rather than being
converted into zero candidates.

An ordinary status query can exit successfully while describing an incomplete
or failed run; it is an observation command. `status --verify` exits nonzero
unless verification succeeds. Use the verified form for a checksum-based script
check, and read the JSON/status detail to understand failures.

## Create an HTML snapshot of one completed sample

```bash
funlr report --run-dir /path/to/run --output /existing/directory/sample-report.html
```

The command first verifies the completed run's recorded output checksums. It
then writes one self-contained HTML file with:

- original, strict and fusion candidate counts;
- tier and flag counts, with multiple flags counted independently per row;
- up to the first 200 review rows in original report order, selected from non-PASS
  flags, fusion hypotheses, `NEEDS_REVIEW` tiers and Tier 2 rescue tiers;
- recorded stage progress/times and warnings;
- available input observations and discovery/rescue evidence summaries;
- recorded elapsed time and output-file sizes;
- the saved settings, input/database hashes, tool versions and inspection record.

Evidence and input-observation summaries are optional for runs created before
those exports were available. The inspector includes them only when their paths
appear in the recorded output checksum mapping. They summarize saved observations
and do not recompute classifications. Output sizes count each recorded path
once, include mirrored final/stage files, and exclude unrecorded working files
and source inputs; they are logical bytes, not filesystem allocation. Elapsed
time is derived from saved timestamps, not measured CPU time or peak memory.
For a resumed run, invocation timing describes the latest invocation.

The report states how many review rows are shown and points to the complete
`final_results/nlr_final_report.tsv` for the remaining evidence. It summarizes
the run's tables; it does not embed the existing domain PNGs or generate new
scientific figures. The HTML contains no scripts, external assets, remote fonts
or network requests, and can be opened locally in a browser.

Choose a new filename ending in `.html` inside an **existing parent directory**.
Existing files are never overwritten. Output symlinks and symlinked parents are
rejected. You may write outside the run, as in the examples, or directly beside
`manifest.json` (for example `run/report.html`). Nested run directories such as
`final_results/` and `results/` are protected so the report does not alter tracked
pipeline products or invalidate their checksums. Missing parent directories are
not created implicitly.

HTML reporting accepts a single completed sample, including
`batch/samples/SAMPLE_ID/`; it does not generate a combined batch HTML page.
Use the batch summary and `status --run-dir batch` for the cohort overview.
A failed validation leaves no new report and returns a nonzero exit status.

The file contains the recorded paths and detailed provenance visible in its
expandable section. Review those contents when choosing what to share. Retain
the original run directory for audit and resume: the HTML is a summary snapshot,
not a replacement for the scientific evidence or independent biological
validation. See [INPUTS_AND_OUTPUTS.md](INPUTS_AND_OUTPUTS.md) and
[BATCH.md](BATCH.md) for the underlying run layout and batch contract.
