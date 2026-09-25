# Synthetic comparison example

These short hand-written sequences, labels and calls test the comparison
harness. They are not real NLR sequences, a biological benchmark or an
independent accuracy estimate. `training_overlap` values are illustrative.

From the repository root:

```sh
funlr benchmark --reference examples/benchmark/reference.json \
  --predictions examples/benchmark/called.tsv \
  --input examples/benchmark/input.faa
```

Two labeled positives are in the input and one is called: recovery is `0.5`.
The absent positive is unavailable. The unresolved called protein is unknown,
not a false positive. Precision and specificity remain null because the labels
are not exhaustive. See [the schema and workflow](../../docs/BENCHMARKS.md).
