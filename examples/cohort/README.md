# Two-row batch workflow fixture

Both rows reference the **same** [public-sequence hybrid fixture](../demo/README.md)
under different sample IDs. This tests sheet parsing, shared configuration,
isolated runs and batch resume. It is not two biological samples, a natural
genome or a comparative benchmark. No extra production databases are needed.

After installing FuNLR and its external tools, run these commands from the source
repository root. Choose a new `cohort-demo` output directory:

```bash
funlr batch --samples examples/cohort/samples.tsv \
  --config examples/cohort/settings.yaml --outdir cohort-demo --validate-only
funlr batch --samples examples/cohort/samples.tsv \
  --config examples/cohort/settings.yaml --outdir cohort-demo --threads 2
funlr batch --samples examples/cohort/samples.tsv \
  --config examples/cohort/settings.yaml --outdir cohort-demo --threads 2 --resume
funlr status --run-dir cohort-demo --verify
funlr status --run-dir cohort-demo/samples/demo_a --verify
funlr report --run-dir cohort-demo/samples/demo_a --output demo-a-report.html
```

TSV paths resolve relative to `samples.tsv`; YAML paths resolve relative to
`settings.yaml`. The first command validates the full sheet without running tools
or creating a batch. Samples execute sequentially into `samples/demo_a/` and
`samples/demo_b/` under the batch output root. Each gets the normal eight-stage
run, `final_results/`, manifest, logs and state. The batch root records the exact
input sheet, resolved sample configurations and batch summary.

The settings match the public fixture: BROAD discovery, no separate Pfam NBD
extraction, no ASM scan, no fusion hypotheses and no protein translation. The
fixture's GFF3 has handwritten spans without CDS features. These settings and
miniature databases must not be copied unchanged into a production analysis.
`PLOIDY=unknown` is only a metadata label; it changes no evidence rules.

With the tested baseline settings, each row yields six report candidates and
three strict candidates. The two executions do not add six independent biological
observations; they rerun identical source data. A successful resume verifies the
same effective configuration and reuses intact stage outputs.

For your own cohort, replace the sample rows with distinct matching genome,
protein and GFF3 files, provide the production HMM libraries, and choose the
appropriate assembly presets. Read [the batch guide](../../docs/BATCH.md) for the
full column and failure-policy contract. Sequence/HMM provenance and reuse terms
remain in [PROVENANCE.json](../demo/PROVENANCE.json) and
[DATA_LICENSES.md](../demo/DATA_LICENSES.md).
