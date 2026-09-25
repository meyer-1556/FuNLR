# Inputs and outputs

FuNLR starts from an annotated genome. Supply the genome, gene models, proteins, and reference profiles used for the intended analysis. It does not assemble genomes, predict a complete gene set, run eggNOG-mapper, or train production HMMs.

## Input contract

| CLI option / YAML key | Expected file and role |
| --- | --- |
| `--genome` / `inputs.genome` | Uncompressed nucleotide FASTA with unique first-token sequence IDs. A masked genome is accepted; supply the same assembly used for the annotation. |
| `--proteins` / `inputs.proteins` | Uncompressed amino-acid FASTA with unique first-token IDs and nonempty sequences. |
| `--gff3` / `inputs.gff3` | Annotation-only GFF3, with mRNA/transcript features, matching genome sequence names, and coordinates within the genome. Embedded FASTA is unsupported. CDS features are needed for protein translation during integration. |
| `--nbd-hmms` / `databases.nbd_hmms` | Fungal NBD HMM collection for whole-proteome discovery. |
| `--pfam` / `databases.pfam` | Uncompressed HMMER-format Pfam-A HMM database, including NACHT and NB-ARC. A precomputed annotation table is not a substitute. |
| `--asm-hmms` / `databases.asm_hmms` | ASM HMM collection. Required with the default `asm.asm_enable: 1`; explicitly set it to `0` to omit this evidence channel. |
| `--custom-hmms` / `databases.custom_hmms` | Optional combined effector-and-sensor HMM collection scanned in Stage 2. Use it to collect the custom domain evidence represented in the reference model snapshot. |
| `--eggnog` / `inputs.eggnog` | Optional eggNOG-mapper TSV with its original `#query` header. Stage 0 reads named annotation columns, including `PFAMs`; these annotations contribute discovery and false-positive evidence. |
| `--annotations` / `inputs.annotations` | Optional Funannotate annotation TSV with a `GeneID` column matching gene IDs. Used for the annotated copy exported at Stage 7. |
| `--effector-hmms`, `--sensor-hmms` / corresponding `databases` keys | Optional separate source libraries recorded in provenance. They are not automatically combined or scanned in place of `custom_hmms`. |

Every supplied path must point to a nonempty regular file. Omitting an optional evidence source can change results. Preserve database versions, model composition, and checksums when comparing runs; see [DATABASES.md](DATABASES.md).

Stage 0 reconciles protein IDs against GFF3 in this order: a matching CDS `protein_id` or `orig_protein_id`; a matching transcript `ID`; then a matching transcript `Parent` gene ID. For a CDS protein attribute containing `|`, the suffix after the final `|` is used. The first CDS mapping for a transcript is retained, and only its first comma-separated parent is used. This mapping policy does not convert arbitrary GFF dialects. `results/stage0/id_mapping_report.tsv` summarizes the mapping routes; `id_coords.tsv` records each transcript's mapping source. Every supplied protein must have a mapping. Additional annotation rows without a protein sequence are reported.

Protein IDs may contain letters, digits, `_`, `.`, `:`, `|`, `+`, and `-`; they cannot start with `-` or equal `.` or `..`. Whitespace separates the FASTA ID from its description. Use matching files rather than independently renaming proteins or dropping unmatched records.

## Configuration and checks

Generate the profile you intend to run, then edit its paths:

```bash
funlr init-config --profile ILLUMINA > settings.yaml
funlr validate --config settings.yaml
funlr run --config settings.yaml --dry-run
funlr run --config settings.yaml --threads 8
```

Use `--profile HIFI` when generating a HIFI configuration. Profile presets supply defaults; explicit scientific settings in YAML/JSON override them. A complete configuration generated for ILLUMINA contains explicit values, so changing only its profile later does not replace those values with HIFI defaults. See [METHODS.md](METHODS.md) for profile differences and settings retained only for provenance.

Paths inside YAML are relative to that file. Explicit CLI options override their corresponding YAML values. The default profile is ILLUMINA, discovery mode is BALANCED, and ASM, fusion rescue, plots, and annotation integration are enabled. BALANCED excludes description-only additions by default; its candidate union includes relaxed NBD, Pfam NBD, and enabled Pfam-enrichment evidence.

`validate` checks required files, FASTA IDs and sequences, GFF structure, coordinate bounds, and protein mappings without external tools or an output directory. Stage parsers additionally enforce their table schemas; for example, a nonempty eggNOG input without `#query` fails at Stage 0. These checks do not establish biological correctness or suitable model calibration. `run --dry-run` performs validation and prints a plan without creating outputs or running tools.

## Descriptive input observations

```bash
funlr input-qc --config settings.yaml
funlr input-qc --config settings.yaml --output input-observations.json
# Or supply only the three matching biological inputs:
funlr input-qc --genome genome.fna --proteins proteins.faa --gff3 genes.gff3
```

`input-qc` validates the input formats and writes JSON to standard output by
default. Optional `--output` writes a new file instead; its parent directory must
already exist, and existing files, output symlinks and symlinked parents are
rejected. The command needs the genome, proteins and GFF3, but no models or
external tools, and creates no run directory. Stage 0 records the same descriptive measurements in
`results/stage0/input_qc.json` after normal preflight. Schema version 1 records
source paths and SHA-256 values, sequence counts and length distributions,
genome N50, N content, lowercase A/C/G/T content, proteins with internal or
terminal `*`, and GFF3 feature/parent counts.

Lowercase A/C/G/T is consistent with soft masking but does not identify the
masking program or repeats. N bases are ambiguous sequence; they do not establish
hard masking, deletion or gap origin. Protein stop characters and initial M are
observations, not completeness or pseudogene classifications. Transcript-parent
counts describe the annotation structure, not independently inferred gene loci
or ploidy. There is no automatic pass/fail quality grade or scientific threshold
change based on these measurements. Percentages state their denominators;
undefined fractions remain null.

## Output layout

Every path below is relative to the chosen output directory. Stages are numbered 0–7. The final directory is `final_results/`.

| Path | Contents |
| --- | --- |
| `results/stage0/` | Clean protein FASTA, protein/contig lengths, coordinate and ID-mapping tables, descriptive `input_qc.json`, eggNOG light table, `master_table.tsv`, and input/software manifests. `contig_lengths.tsv` has no header. |
| `results/stage1/` | Raw NBD/Pfam-NBD scans; strict and relaxed ID lists; hit tables; Pfam enrichment and description audits; `union_candidate_ids.txt`, `union_candidates.faa`, and union accounting. `nbd_candidate_ids.txt` is an additive alias for the strict list. |
| `results/stage2/` | Raw Pfam/custom/ASM scans, `all_domain_hits.tsv`, the 30-column `architecture_summary.tsv`, `asm_hits_summary.tsv`, enrichment, and parser/filter reports. |
| `results/stage3/tables/` | `tiered_candidates.tsv`, `rescue_priority.tsv`, and `TIER_*.tsv` for tiers present in the run. Adjacent `beds/`, `summaries/`, `qc/`, and `metadata/` retain coordinate exports, counts, audits, and provenance. |
| `results/stage4/` | `TIER_*.ids`/`.faa` for discovered tier tables, `tiered_candidates.ids`/`.faa`, `rescue_priority.ids`/`.faa`, and export manifests/counts. There is no fixed eight-file tier list. |
| `results/stage5/` | Priority miniprot queries, GFF, parsed tables, summaries, QC, metadata, and Exonerate handoffs. The comprehensive track has the same layout under `comprehensive/`. |
| `results/stage6/` | Priority Exonerate queries, raw alignments, merged GFF, summaries, fusion proposals/evidence, and metadata. The comprehensive track is under `comprehensive/`. |
| `results/stage7/` | Final-report working tables, summaries, BED files, plotting/integration artifacts, and status metadata. |
| `final_results/nlr_final_report.tsv` | Candidate report joined with priority and comprehensive rescue evidence, architecture labels, and audited fusion decisions. |
| `final_results/nlr_strict_candidates.tsv` | Combined strict report; `.priority.tsv` and `.comprehensive.tsv` preserve the respective track-specific sets. These remain computational candidates. |
| `final_results/*.bed` | Original and available rescue/fusion coordinates, in zero-based half-open BED format. TSV/GFF coordinates remain one-based closed. |
| `final_results/domain_plots/` | Enabled domain plots, `domain_architecture_summary.tsv`, and Matplotlib/font provenance (or R session information for the optional R backend). Empty plot inputs are reported as skipped rather than producing fabricated plots. |
| `final_results/evidence/` | Candidate/channel accounting, review and exclusion reasons, raw Stage 1 per-hit evidence, exact aligned protein segments and `evidence_summary.json`; described below. |
| `final_results/integration/` | Gene-call tables, Tier 1/2 gene-ID lists, isoform-conflict and export-status audits, and enabled `.withNLR` annotation/GFF/eggNOG copies and translated protein FASTA. Existing gene models are tagged; rescue/fusion models are not inserted into GFF3. |
| `manifest.json`, `resolved_config.json`, `state.json`, `logs/`, `work/` | Input/tool/code identities, effective scientific settings, verified stage outputs, commands/logs, private input/database copies, and indexes. Retain these for audit and resume. |

The priority rescue track follows `rescue_priority.faa`, including a valid empty file. Comprehensive rescue considers Tier 1–3 candidates independently. Thus an empty priority set does not imply that no alignment tool will run. See [RESCUE.md](RESCUE.md) for rescue gates, fusion acceptance, strict sets, and export limitations.

Empty candidate and rescue results are valid after successful searches. Expected tables keep headers and expected empty FASTAs remain empty. Tier-specific files are created for tiers present, so consult export manifests instead of assuming every possible tier exists. Checked tool failures stop the run. The two explicit fallback cases—failed Pfam NBD extraction and a failed optional proteome-background Pfam scan—are logged and documented in [METHODS.md](METHODS.md).

## Evidence exports

Stage 7 writes these files under `final_results/evidence/`, with a working copy
under `results/stage7/evidence/`. They expose saved evidence without changing
candidate unions, tiers, rescue selection or strict sets.

| File | Interpretation |
| --- | --- |
| `candidate_evidence.tsv` | Every final-report ID plus IDs present in saved discovery-channel lists or raw Stage 1 hits. Includes excluded evidence IDs; it is not an inventory of every input protein. |
| `review_and_exclusions.tsv` | A subset of that table with explicit review or exclusion notes. Multiple reasons can apply to one ID. |
| `nbd_hit_evidence.tsv` | Parseable raw fungal-NBD and Pfam-NBD domain-table rows with source file/line, original model/accession, HMM and query coordinates, scores, E-values, original threshold pass flags and segment extraction status. |
| `nbd_aligned_segments.faa` | Exact one-based inclusive query-alignment spans from the cleaned Stage 0 protein FASTA, limited to hits accepted by the original recorded filters with valid coordinates. Per-hit identifiers link to the evidence table. |
| `evidence_summary.json` | Versioned source identities, counts, discovery overlaps, annotation-only IDs, recorded rescue handoffs and final-report/strict-count reconciliation. Missing evidence is labeled unknown or partial. |

Review reasons include `NOT_IN_CANDIDATE_UNION`, `NOT_IN_STRICT_SET`,
`MEMBER_OF_FUSION_HYPOTHESIS`, `FLAGS_REQUIRE_INSPECTION`,
`FUSION_REQUIRES_REVIEW` and `REVIEW_OR_RESCUE_TIER`. They explain recorded
membership and flags. A review row is not automatically a false positive, and
membership in a fusion hypothesis does not validate a reconstructed gene.

The hit table retains all parseable raw hits. FASTA export requires acceptance
by the original recorded strict/relaxed fungal-NBD or Pfam-NBD filters and valid
coordinates. The table's `accepted_by_recorded_filters` and `segment_status`
fields distinguish accepted, rejected, unassessed and invalid mappings. Rows
marked `NOT_ACCEPTED_BY_RECORDED_FILTERS` or `NOT_ASSESSED` have no exported
sequence or segment hash. Missing acceptance evidence is not treated as a pass.

Accepted hits can still represent competing model families: these filters do
not independently establish domain identity. Duplicated or overlapping sequences
are expected when searches report overlapping accepted hits. They are not merged
full NBDs, recovered fusion sequences or a ready-to-use multiple sequence
alignment. Check model identity, pass flags and segment status when selecting
sequences. Invalid coordinates or query-length mismatches retain a mapping row
with an error status and no FASTA record; coordinates are never silently clipped.
Stage 1/2 summary intervals can span multiple hits and should not be substituted
for individual alignment coordinates.

Summary counts distinguish `reported_hits`, `accepted_reported_hits`,
`unassessed_reported_hits` and `exported_segments`. The `scan_status` record says
whether optional Pfam-NBD scanning was requested and completed. Zero optional
hits can reflect a disabled or uncompleted scan and must not be interpreted as
biological absence. Input observations remain in the separate Stage 0 QC JSON.

The HTML and JSON status summaries can expose these recorded observations;
[REPORTS.md](REPORTS.md) explains verification and limits. Independent curated
reference evaluation is described in [BENCHMARKS.md](BENCHMARKS.md).

## Troubleshooting

- **Missing tool or library:** run `funlr doctor --demo` in the same environment and compute allocation. Default production runs use samtools, hmmpress, hmmfetch, hmmscan, SeqKit, miniprot, Exonerate, gffread, and Matplotlib for headless Python plots. Rscript and its packages are needed only for the optional `reporting.plot_backend: r` setting. The installation demo also needs hmmbuild.
- **ID mismatch or empty protein export:** inspect Stage 0 mapping tables, FASTA first-token IDs, transcript/CDS parents, and whether the GFF contains valid CDS features for gffread. The annotation TSV's `GeneID` values must match genes, not independently renamed proteins.
- **Missing ASM database:** provide `--asm-hmms` or explicitly disable `asm.asm_enable`. The portable default does not silently discard this enabled channel.
- **Changed resume signature:** start a new output directory after changing inputs, paths, code, configuration, runtime, or tools. `--resume` verifies an identical run and recomputes damaged stages plus their dependents.
- **Interrupted run:** inspect the failing stage log and scheduler output, then resume with the same command and `--resume`. Normal SIGINT/SIGTERM handling stops worker tools and removes the run lock. After a forced kill, verify the recorded process has stopped before removing a stale `.run.lock`.
- **Space or resource exhaustion:** allow for private genome and HMM copies/indexes plus alignment outputs. `--ram-gb` records an allocation label; the scheduler or operating system enforces memory limits.

For a bug report, provide the FuNLR version, effective configuration, command, doctor output, and relevant stage log. Share a small reproducible input example when possible; redact private paths or data as needed.
