# Validation of 0.5.0a2

## Alpha release checks on 24 September 2026

- **589 tests passed, 3 skipped and 4 subtests passed** on Python 3.12.14. Source files were normalized to mode 0644 to represent browser uploads. The existing pandas warning remains; skip scope is recorded below.
- The wheel was installed in a fresh Python environment outside the checkout. `pip check` and the packaged-fixture self-test passed. Every packaged application file matched the release source.
- Both installed-wheel demonstrations completed all eight stages and verified eight-stage resume with real external tools. The public-sequence fixture produced 6 original candidates and 3 strict candidates; the synthetic fixture exercised positive Exonerate rescue.
- Verified status matched all 252 recorded public-demo outputs. A completed-run HTML report was exported successfully; no browser visual review was performed in this check.
- All 48 runtime Python modules outside the package-version module have identical parsed syntax trees to 0.5.0a1. The only runtime-file edits are the package version/docstring and one comment. Classification and execution policy are retained.

The [alpha check record](validation/ALPHA_RELEASE_CHECKS.json) records versions, hashes, fixture identities and scope. These checks used macOS on Apple Silicon with existing native tools and the documented serial-build Exonerate workaround. They do not establish a clean Conda-lock installation, Linux container execution or hosted-CI success. No new full-genome run or independent biological accuracy evaluation was performed. Production model distribution remains separate.

The earlier records below retain their original dates and software versions. In particular, cohort and saved-dataset replay results are earlier evidence, not repeated claims for this alpha preparation.

# Earlier validation of 0.5.0a1

The checks below were completed on 12 September 2026. Software correctness, installation behavior, reference-output agreement and independent biological accuracy are separate claims.

## Software and installed-package checks

- **589 tests passed, 3 skipped and 4 subtests passed** on Python 3.12.14. The complete suite ran from a source copy with every file normalized to mode 0644, matching browser-upload behavior. The skips are two optional real-tool pytest demos, run separately below, and one optional external reference-data test. One existing pandas future-downcasting warning remains under the pinned pandas version.
- A wheel was installed into a fresh Python environment and used outside the source tree. `pip check` and `twine check` passed. Package fixture checksums and new model metadata were verified. Wheel runtime `.py`, `.R` and `.json` bytes match the current source package.
- The synthetic real-tool demo completed all eight stages, verified reuse of all eight, and exercised positive Exonerate rescue.
- The public-sequence real-tool demo completed all eight stages and verified eight-stage reuse: **6 original candidates, 3 strict candidates, 1 priority and 3 comprehensive refinement queries**. It retains the expected 49-bp intron-bearing locus alignment.
- A two-row installed-wheel cohort completed both runs, verified both saved stage states during resume, passed checksum inspection, produced two HTML reports whose displayed review rows matched the source tables, and preserved no-write planning behavior. Its two samples are copies of the same constructed fixture, not a biological cohort.
- Installed `input-qc`, `model-inventory`, `benchmark` and verified `status` commands passed explicit output checks. The synthetic reference-panel example returns positive recovery 0.5 and null precision/specificity, as intended for its incomplete labeling.

The machine-readable [release record](validation/RELEASE_CHECKS.json) contains counts, runtime identities, tool versions and checksums. Tests cover incomplete/zero evidence, malformed inputs, wrong evaluation units, unavailable references, unlisted predictions, guarded metrics, individual HMM-hit coordinates and accepted-only segment exports. A regression test corrupts the new evidence table and confirms checksum detection and Stage 7 repair through resume.

These native checks used macOS on Apple Silicon, HMMER 3.4, samtools 1.21, SeqKit 2.12.0, miniprot 0.18, gffread 0.12.7 and Exonerate 2.4.0 with the documented **serial-build workaround**. They used existing external executables on PATH and a fresh Python package installation. They are not evidence that the unmodified Apple Silicon Conda lock installs and executes cleanly. See [installation limitations](INSTALLATION.md).

## Scientific decision preservation

An independent static audit compared 35 scientific modules and 181 top-level functions with version 0.4.0a1. It preserved numeric constants, regexes, scientific strings, command arguments and control flow, normalizing only specifically reviewed internal imports, documentation and diagnostic text. All 81 scientific defaults and 12 effective profile/mode/override configurations matched. The original model identities, bundled profile bytes and existing report calculations also matched. [Audit scope and hashes](validation/SCIENTIFIC_CHANGE_AUDIT.json).

Stage 0 adds descriptive QC after input processing; Stage 7 adds evidence exports after the existing final report. These artifacts are now required for stage completion and recorded for resume. Therefore the release preserves classification decisions while intentionally changing output/completion requirements. Static equivalence and synthetic regression tests do not establish biological accuracy.

## Saved full-dataset evidence replay

The new exporter was applied read-only to an existing completed full native run produced by version 0.2.0a2. This release **did not rerun the full native searches or rescue tools** on that dataset.

The replay audited **5,918 individual HMM hits**: **1,397 accepted hits** produced exact sequence segments; 4,521 rejected hits remained in the mapping table without FASTA records. Every exported sequence matched its bounded alignment span in the scanned protein FASTA and its recorded checksum. There were no unassessed hits or export warnings. The evidence reconciled 257 original candidates, one fusion review hypothesis and 55 strict rows. Sizes and modification times of 508 source files remained unchanged.

The [replay record](validation/EVIDENCE_REPLAY.json) binds the exporter implementation and generated outputs. [Dataset identities](validation/REFERENCE_IDENTITIES.json) distinguish saved input/database hashes from replay-read intermediate hashes. This is evidence of correct export/accounting on a real-sized saved analysis, not a new full-pipeline equivalence test or independent fungal benchmark.

## Model and release boundaries

The inventory was exercised against all five exact production libraries: 39 ASM, 41 NBD, 127 combined custom, 19 effector and 108 sensor profile records. The 334 total counts include repeated records in the combined and component libraries. The combined custom file is the exact effector file followed by the sensor file. These are byte-identity findings, not proof of training history, biological specificity or redistribution permission. See [model inventory](MODELS.md).

A new native full-genome run, independently curated fungal positives and hard negatives, HiFi evaluation, clean installation testing on all advertised platforms and measured production resource guidance remain necessary. Docker/Apptainer recipes and hosted CI are included, but were not executed remotely or published during this release preparation. Conda packages, registry images and a production model download bundle are not yet published. The [roadmap](ROADMAP.md) states the acceptance criteria.
