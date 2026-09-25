# Curated reference comparison

`funlr benchmark` is a read-only comparison harness. It compares an explicitly
selected called-ID set with a labeled reference and checks which reference IDs
occur in the evaluation input. It does not run FuNLR, curate biological truth,
validate gene models or establish an independent biological benchmark.

```sh
funlr benchmark --reference examples/benchmark/reference.json \
  --predictions examples/benchmark/called.tsv \
  --input examples/benchmark/input.faa
```

The included example is entirely synthetic. Its hand-written sequences, labels
and calls demonstrate denominator handling; they are not accuracy evidence.
Expected output includes positive recovery `0.5`, one unavailable positive,
one unresolved call, and null precision and specificity.

## Define one evaluation unit

The reference is JSON with exactly these top-level fields:

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer `1` |
| `benchmark_id` | Nonempty single-line identifier for the reference version |
| `unit` | `protein` or `locus` |
| `prediction_set` | Description of the positive-call set and selection policy |
| `input_sha256` | SHA256 of the exact evaluation FASTA bytes |
| `exhaustive_labels` | Boolean: every input ID has a resolved positive or negative label |
| `synthetic` | Boolean; must be true when any record has synthetic evidence |
| `records` | Nonempty list of the records below |

Each record contains exactly these fields:

| Field | Meaning |
| --- | --- |
| `protein_id` **or** `locus_id` | Use the field matching the declared unit, never a fallback |
| `expected_call` | `positive`, `negative` or `uncertain` |
| `evidence_category` | `experimental`, `curated`, `computational` or `synthetic` |
| `evidence_source` | Single-line citation/curation reference; required for positive/negative labels, otherwise nullable |
| `reference_present` | Boolean, checked against exact first-token FASTA IDs |
| `training_overlap` | `yes`, `no` or `unknown`, supplied by the curator |

The example [reference.json](../examples/benchmark/reference.json) provides a
complete editable template. Unknown fields and duplicate IDs are rejected.
Uncertain is an explicit label; an unlisted record is also unresolved for
scoring purposes. Neither becomes a negative automatically.

The input FASTA defines the evaluation universe, not necessarily an entire
proteome. For protein evaluation, use the same original-protein FASTA used to
produce the calls. For locus evaluation, supply a deliberately constructed
locus-ID FASTA and pass `--id-column locus_id`. Its mapping to genes, isoforms,
fragments and fusions is a curator assertion that this harness cannot validate.
Do not pass genomic contig IDs as though they were locus IDs. Separate protein
and locus analyses avoid counting multiple fragments as independent recovered
loci or pooling gene IDs with protein IDs.

## Export a called set deliberately

The predictions TSV must have exactly one column, named `protein_id` (default)
or `locus_id`. It contains only positive calls under the declared selection
policy. A header-only file is a valid empty call set. Duplicate or out-of-input
IDs fail before any scores are produced. A full final report is refused: its
rows include candidates that are not accepted positive calls.

For the original-protein strict set, first verify the completed run:

```sh
funlr status --run-dir results --verify --json
python examples/benchmark/export_strict_proteins.py \
  results/final_results/nlr_strict_candidates.tsv called.tsv
```

Proceed only if the status reports completed results and verified integrity.
The export script requires the strict-table filename and excludes fusion-model
rows. It creates a new output exclusively and refuses to overwrite one. This
original-protein export therefore does not measure fusion or gene-model rescue
accuracy. Those require a separately adjudicated locus-level reference and
mapping. The comparator records hashes of the reference, input and called set,
but cannot verify their upstream run from a detached ID list; its result states
`run_integrity: NOT_CHECKED` explicitly.

## Read the metrics with their denominators

Positive recovery is recovered labeled positives divided by labeled positives
present in the evaluation FASTA. An absent reference is `UNAVAILABLE`, never a
false negative. A present labeled negative contributes to the negative-panel
call count. Calls with uncertain or unlisted labels appear in
`unknown_called_ids` and never contribute false positives.

Precision and specificity remain null for partial or positive-only panels.
They are available only when `exhaustive_labels=true`, every evaluation FASTA
ID actually has a positive/negative label, and the panel includes negatives.
They describe agreement with those labels in that declared universe, not a
larger proteome or all fungi. Ratios with a zero denominator are null; zero
recovered positives with a nonzero denominator give recovery `0.0`.

Results retain per-record outcomes and positive-recovery strata by evidence
category and training overlap. Training overlap is an assertion, not detected
by the harness. Non-synthetic comparisons still report independence as
`NOT_ESTABLISHED_BY_THIS_HARNESS`. A literature-derived label or sequence can
overlap the HMM training data even when identifiers differ.

## Current evaluation limits

The included synthetic fixture checks metric accounting. It does not measure biological accuracy. A small positive panel establishes recovery only; precision over a complete proteome requires sufficiently complete labels for the declared evaluation universe. Reference provenance, input presence and label quality remain the responsibility of the person supplying a panel.
