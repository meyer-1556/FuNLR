# Rescue, fusion evidence and final reporting

FuNLR runs native miniprot, Exonerate and gffread through its checked command
runner. The default plots use headless Matplotlib, with a separately configured
R backend available. A nonzero checked tool exit fails the stage even if partial
results already exist. Command logs retain executable paths, arguments and exits.

## Two rescue tracks

Stage 5 writes the priority track under `results/stage5/` and the comprehensive
track under `results/stage5/comprehensive/`. Priority queries come from the
Stage 4 `rescue_priority.faa` export. An existing empty export means no priority
queries. The fallback for older exports combines Tier 2 and Tier 1B FASTAs,
with Tier 3 added when `rescue.include_tier3_in_rescue` is enabled. Optional
`rescue.filter_low_nbd` filtering uses the Stage 3 `nbd_ok` field, or HIGH/MEDIUM
confidence when that field is absent. Comprehensive queries are every Stage 3
Tier 1–3 protein, extracted by exact identifier from the cleaned Stage 0 FASTA.

Each track retains query FASTAs and IDs, `gff/miniprot.gff3`,
`tables/miniprot_features.tsv`, `tables/miniprot_summary.tsv`, hit counts,
missing/extra-hit IDs, and `qc/exonerate_refine_reasons.tsv`.
`queries/exonerate_refine_ids.txt` passes candidates to Stage 6. Candidates
without a miniprot hit or mRNA feature need refinement; mapped candidates need
refinement when they have at least 12 CDS features, an mRNA span of at least
15,000 bp, score at most 20, or mappings on multiple scaffolds. These are fixed refinement-rule thresholds, not measures of model accuracy. Tier, confidence and flag annotations add reasons
after one of these mapping conditions is met.

Miniprot statistics are read from the native `Rank`, `Identity` and `Positive` GFF attributes.
Query association comes from `Target`. The summary convention—first
scaffold with aggregate minimum/maximum coordinates—is retained, including for
multiple mappings. Raw feature rows remain available for inspection; a summary
span alone is not evidence of one contiguous model.

Stage 6 uses the same track layout. It processes up to
`rescue.exon_maxn` original refinement queries per track (zero means
all), then the fusion hypotheses. Protein-to-genome alignment retains the
configured Exonerate score, percent and maximum-intron settings. Outputs include
`queries/exonerate_all.ids`, `queries/exonerate_refine_ids.work.txt`, query FASTAs, `gff/per_query/`,
`gff/exonerate_merged.gff3`, `tables/exonerate_summary.tsv`, and
`tables/fusion_manifest.tsv`/`fusion_evidence.tsv`. Merging uses only explicit
query/output pairs from the current invocation and annotates their query IDs.
Both stages write headed tables when there are no queries or hits, and skip
the unnecessary external alignment command. Each track also has a manifest and
`metadata/run_info.json`.

`EXON_MAX_INTRON_RETRY` is recorded but inactive: zero for ILLUMINA and 10,000 for HIFI. No second alignment pass is implemented. Changing this recorded setting does not enable a retry. `rescue.exon_max_intron` controls
the maximum intron length of the actual alignment command.

## Fusion hypotheses and the acceptance guard

Fusion construction uses a connected-component heuristic, confidence requirement, exclusion flags, member limit and gap settings. Defaults allow up
to three members, a 50 kb ILLUMINA or 20 kb HIFI base gap, extension up to 100 kb
for hard-end evidence or 50 kb for repeat evidence, and MEDIUM or higher NBD
confidence. The configured overlap limit also applies. A confident NBD member
without a sensor is paired with a sensor member without a confident NBD.

Members are ordered in transcript direction: ascending genomic coordinates
on `+`, descending on `-`. Fusion detection runs in both tracks even if there are no original
queries requiring refinement. Track-specific raw IDs are retained in manifests.

A successful Exonerate command or an alignment to just one concatenated member
does **not** justify replacing the original candidates. `fusion_evidence.tsv`
requires one valid vulgar alignment, on the expected scaffold and strand, whose
matched alignment operations overlap **every query member** and **every original
member locus**. Matched blocks are checked individually; bounding spans and
separate alignments cannot fill an unaligned member. Malformed operations and
inconsistent consumed lengths fail this check. Rejection reasons are retained.

This is an explicit necessary evidence guard, not a biologically validated
fusion classifier. It introduces no arbitrary percentage-coverage cutoff.
Accepted rows remain review hypotheses, labelled
`FUSION_SPANNING_ALIGNMENT;REVIEW_REQUIRED`; alignment does not establish a
functional NLR, a correct gene model, or biological fusion. Query span fraction
and the evidence rule are recorded so the decision is inspectable.

The final report collapses accepted hypotheses with the same member set across
tracks. A stable `FUSION_` identifier derives from a hash of sorted member IDs;
track names and raw aliases remain in `fusion_source` and `fusion_raw_id`.
The best accepted alignment supplies report coordinates. Unsupported hypotheses
never mark their original members as replaced. Ordinary `fused_into` values are normalized to empty strings; missing values do not remove unrelated strict candidates.

## Final tables, tracks and annotation exports

`final_results/` contains `nlr_final_report.tsv`,
`nlr_strict_candidates.tsv`, `nlr_strict_candidates.priority.tsv`,
`nlr_strict_candidates.comprehensive.tsv`, `tier_summary.tsv`, `qc_summary.tsv`,
and six BED tracks. The corresponding working tables, summaries and BEDs remain
under `results/stage7/`.

The combined strict set includes Tier 1A/1B, 2A/2B, 3A/3B, and Tier 2D with HIGH
NBD confidence, plus accepted fusion review hypotheses. Original members of
accepted hypotheses are omitted from strict sets but remain in the full report.
Priority strict candidates additionally require a priority miniprot or
Exonerate mapping; comprehensive strict candidates require a comprehensive
miniprot mapping. Each track includes its own accepted fusion hypotheses. `NON_NLR` flags do not
independently remove a row that otherwise satisfies its tier rule.

The report provides a dominant-span, conceptual
effector–NBD–sensor architecture summary using domain-name patterns.
This summary is not a reannotation of the recovered alignments.

With `reporting.integration` enabled, `final_results/integration/` contains NLR
calls, gene lists, counts, optional annotation/eggNOG joins, and tagged GFF3 and
translated protein exports. GFF3 tagging adds NLR metadata to existing mRNA
features. It does not insert rescued Exonerate or fusion CDS models. gffread
translates the original annotated structures; resulting headers carry NLR
metadata. The first report isoform supplies the gene-level call, with disagreements exposed in `isoform_conflicts.tsv`.
`export_status.tsv` records written and skipped exports.

`reporting.write_tagged_gff3`, `write_tagged_eggnog` and `write_proteins` control
their respective exports. Protein export requires tagged GFF3; it is skipped
if GFF3 export is disabled. The Funannotate annotation-table join needs optional
`inputs.annotations` with a `GeneID` column. Disabling integration skips all
these exports and removes the need for gffread.

## Plots

`reporting.plots` enables the Stage 7 renderer. Its six summary PNGs,
individual Tier 1–3 protein diagrams and `domain_architecture_summary.tsv` are
written to `domain_plots/` in both Stage 7 and `final_results/`. The renderer uses the integrated report's domain taxonomy. The
length-versus-priority plot uses Stage 3 protein lengths rather than treating the
last domain coordinate as protein length. Fusion diagrams are not synthesized.

The default `reporting.plot_backend: matplotlib` uses a headless Agg renderer
and a private plotting cache. Its
`plot_provenance.json` records the actual Matplotlib version, DPI, backend, cache
and font identities. The Python and R backends can produce different pixels; the scientific tables retain the same interpretation.

The optional `reporting.plot_backend: r` uses the R renderer and writes `R_sessionInfo.txt`. It requires a separately supplied R runtime
and packages, which are checked rather than installed during a run. R PNG export
uses the base device, Quartz on macOS and Cairo on Linux, with a white background.
The default installation does not require R; see [INSTALLATION.md](INSTALLATION.md)
for the optional environment.

Both paths check expected PNG files and signatures. Plot destinations are
refreshed to prevent stale images from earlier attempts. With no eligible domain
hits, plotting records a skip in `metadata/plot_status.json`. The separate
`funlr report` command creates an HTML overview of completed results; it does not
rerun this plotting stage or modify candidate evidence. Its interface and
verification scope are documented in [REPORTS.md](REPORTS.md).

## Validation limits

Tests reject alignments confined to one fusion member and exercise acceptance
when one valid alignment spans every required query member and source locus.
Other checks cover reverse-strand ordering, query gaps hidden by bounding spans,
wrong loci, cross-track duplicates, strict-candidate preservation, empty tracks,
native miniprot attributes and failed graphics output. These check software
behavior and evidence accounting; biological sensitivity and specificity require
separate evaluation. See [VALIDATION.md](VALIDATION.md).
