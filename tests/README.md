# Test scope

Install with `python -m pip install -e '.[dev]'`, then run `python -m pytest tests -q` from the repository root.

Unit and integration tests cover input/configuration contracts, threshold and tier boundaries, HMM parsing, deterministic selection, both rescue tracks, fusion-review evidence, literal FASTA exports, output verification, failure recovery, no-write planning and installed-resource checks. Test doubles are explicitly enabled in tests and are never installed in the application.

`test_differential.py` runs the same 12-protein fixture in independent processes with different hash seeds, output directories and equivalent configuration representations. It compares selected scientific files exactly and confirms that a changed classification is detected. This measures reproducible execution of the current rules, not biological accuracy.

Evidence/QC tests cover descriptive denominators and input integrity, accepted raw HMM-hit coordinates/scores, empty exports, excluded/strict membership, incomplete reference panels, unavailable labels, evaluation units and model provenance. The included reference-panel fixture is synthetic.

With actual external tools installed, run both packaged demo commands. `real_tools_smoke.py` also checks positive Exonerate rescue and resume; `real_tools_batch.py` checks two independent sample-sheet runs, status, report contents and verified resume using the public fixture. The pytest demo cases are enabled with `FUNLR_TEST_REAL_TOOLS=1`; they otherwise skip and are run separately in installation CI. `FUNLR_REFERENCE_DIR` optionally enables a reference-data comparison using separately supplied inputs.

These checks do not establish sensitivity, specificity, reconstructed-gene correctness or support for untested platforms. Record source tests, installed-package demos, production dataset comparisons and independent biological benchmarks separately. See [validation](../docs/VALIDATION.md).
