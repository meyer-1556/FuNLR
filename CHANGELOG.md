# Changelog

## 0.5.0a2 — public alpha preparation

- Identify the portable, rule-based pipeline explicitly as an alpha in the README, package metadata and release notes. Distinguish implemented features, installation demonstrations and biological validation.
- Clarify that production model assets and public package/container distribution are separate from the bundled examples.
- Provide browser-only source upload, repository visibility and tagged prerelease instructions.
- Clean remaining internal terminology while retaining scientific credits and third-party notices.
- Retain the 0.5.0a1 scientific decisions, execution safeguards, packaged demos and evidence utilities.

## 0.5.0a1 — evidence review and reference comparisons

- Add descriptive genome, protein and annotation QC through `input-qc` and automatic Stage 0 JSON output. Metrics do not filter or grade inputs.
- Add candidate membership and review/exclusion tables, original per-hit NBD evidence and aligned segment FASTA exports. Preserve raw HMM coordinates and scores instead of treating aggregate intervals as complete domains.
- Add read-only per-profile model inventory with known source metadata, explicit unknowns and optional fixed-snapshot verification.
- Add a checksum-bound reference-panel comparison command with explicit protein/locus units, unavailable-reference tracking and conservative metric definitions for incomplete panels.
- Extend status and HTML review with saved evidence, QC and recorded elapsed-time/output-size summaries.
- Clarify optional-input diagnostics to distinguish missing scientific evidence from unavailable annotation joins and separate model provenance.
- Retain discovery thresholds, tier decisions, rescue rules, default profiles, verified resume, checked subprocesses and headless plotting.
- Consolidate public documentation around features and scientific contracts. Track distribution and validation work in the roadmap.

## 0.4.0a1 — sample sheets and installation diagnostics

- Add case-insensitive sample-sheet columns, assembly-profile aliases, per-sample discovery modes and explicit masked-genome fallback.
- Add installed-resource self-tests and optional production-model identity checks.
- Clear stale pipeline annotation tags when re-tagging existing annotations.

## 0.3.0a1 — portable execution and reporting

- Provide staged Python execution, explicit configuration, TSV cohorts, private indexes, checked commands, input/tool identities and checksum-verified resume.
- Provide dual genome-rescue tracks, guarded fusion review, tagged annotation exports, headless domain plots and offline HTML reports.
- Package synthetic and public-sequence demonstrations, dependency locks and source-based installation/container recipes.
