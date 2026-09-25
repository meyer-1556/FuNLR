# Public-example verifier regression data

These small tables and raw miniprot alignments were copied from an actual
native-tool run of the bundled public-sequence fixture on 2026-09-09, using
miniprot 0.18-r281 and Exonerate 2.4.0 with the current stages. They exercise
the public-demo output checker without rerunning external tools in unit tests.

The originating run completed stages 0–6 and report-table generation; its later
GFF translation export failed because the handwritten input has no CDS features.
The corrected demo disables that export and is separately checked end to end.
These files alone do not establish successful installation or biological
validation. Input-source attribution and reuse terms are in
`examples/demo/DATA_LICENSES.md` and `PROVENANCE.json`.
