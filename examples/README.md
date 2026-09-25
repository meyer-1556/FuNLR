# Run with your own data

After installing FuNLR and its external tools, run the existing stages in order:

```bash
funlr run \
  --genome /path/to/genome.fa \
  --proteins /path/to/proteins.fa \
  --gff3 /path/to/annotation.gff3 \
  --nbd-hmms /path/to/NLR_NBD_combined.hmm \
  --asm-hmms /path/to/NLR_ASM_combined.hmm \
  --pfam /path/to/Pfam-A.hmm \
  --eggnog /path/to/sample.emapper.annotations \
  --outdir /path/to/new-funlr-run \
  --threads 8
```

Use uncompressed FASTA and annotation-only GFF3. Protein identifiers must map through CDS protein attributes, transcript IDs or parent gene IDs as specified in [the input contract](../docs/INPUTS_AND_OUTPUTS.md). Quote paths containing spaces. Choose a fresh output directory for each run.

Add `--custom-hmms /path/to/NLR_custom_combined.hmm` when available. Omit `--eggnog` if you do not have those annotations; this removes the eggNOG priority evidence and can change the candidate set and tiers.

No configuration file is required to retain the pasted scripts' defaults. `thresholds.json` shows the existing threshold names and the April 2026 ILLUMINA values. Passing `--config examples/thresholds.json` makes those selected defaults explicit; unspecified parameters retain their defaults.

The exact NBD/ASM models and Pfam database are scientific inputs that you must supply. This directory does not substitute synthetic models or toy predictions for them. The executable doubles and placeholder HMMs used by the integration tests are documented separately in `tests/README.md`.
