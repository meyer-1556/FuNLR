# Discovery, architecture and classification methods

FuNLR combines domain searches, annotation evidence, protein architecture and
assembly context to classify computational fungal NLR candidates. The default
configuration uses the ILLUMINA assembly profile and BALANCED discovery mode.
Preserve the effective settings and model identities when comparing analyses.
Profile names describe configurable assumptions; they do not establish assembly
quality or biological accuracy.

This document covers Stages 0–4. See [RESCUE.md](RESCUE.md) for Stages 5–7,
rescue tracks, fusion evidence, strict reports and annotation export limitations.

## Stages 0–4

| Stage | Scientific behavior |
| --- | --- |
| 0 — input tables | Protein FASTA IDs are first tokens. Transcript coordinates are mapped through matching CDS protein attributes, transcript IDs, then parent gene IDs. The first CDS mapping is retained. The coordinate table defines the master-table universe. eggNOG columns are read from `#query`, including `PFAMs`; annotation text is assembled from the same fields. |
| 1 — candidate discovery | Whole-proteome fungal NBD discovery has strict and relaxed parses. Optional Pfam NACHT/NB-ARC discovery, enrichment of eggNOG PFAM annotations among strict NBD candidates, and audited description keywords feed the configured candidate union. The union is intersected with supplied protein IDs by default. |
| 2 — domain architecture | Candidate proteins are scanned against Pfam and the supplied combined custom library. N-terminal ASM windows form a separate evidence channel. Filtered hits feed the sensor/effector groups, conservative NBD calls, confidence, interval/orientation logic, false-positive flags, and 30-column architecture summary. |
| 3 — tiers and flags | Profile-aware tiering uses the decision order: raw NBD evidence, LRR/orientation exclusions, low-confidence housekeeping and fragment handling, canonical sensor categories, rescue context, and remaining architectural variants. Flags and score terms use the rules below. |
| 4 — protein exports | Present `TIER_*.tsv` files are discovered dynamically, and their proteins plus all tiered candidates and the priority rescue subset are exported. Optional rescue FASTA confidence filtering remains off by default. |

The Stage 0 reconciliation policy is specific: the suffix after the final `|` in a CDS `protein_id`/`orig_protein_id` is considered; the first comma-separated CDS parent is used. A matching CDS ID takes precedence over a matching transcript ID, then a matching gene ID. Preflight requires every supplied protein to map and validates GFF coordinates. It reports annotation rows without a protein sequence. It does not silently discard unmatched proteins.

### Discovery defaults

| Setting | Default and effect |
| --- | --- |
| `EVAL_NBD`, `MIN_NBD_ALI_LEN` | Strict domain i-Evalue ≤ `1e-5`, aligned length ≥ 180 aa. |
| `EVAL_NBD_RELAXED`, `MIN_NBD_ALI_LEN_RELAXED` | Relaxed domain i-Evalue ≤ `1e-3`, aligned length ≥ 120 aa. |
| `PFAM_NBD_SCAN` | `1`: fetch and scan Pfam NACHT and NB-ARC; parse with relaxed thresholds. |
| `DISCOVERY_MODE` | `BALANCED`; description-only matches do not expand the union by default. |
| `USE_PFAM_ENRICHMENT_FOR_UNION` | `1`; eligible enriched PFAM annotations can add proteins. |
| `PFAM_ENRICH_MIN_FG_COUNT`, `PFAM_ENRICH_MIN_LOG2FC`, `PFAM_ENRICH_TOP_N` | At least 2 strict-foreground proteins, log2 fold change ≥ 1.0, at most 50 enriched PFAM names. Foreground counts use deduplicated hit protein IDs; background counts use master-table rows. Repeated PFAM labels are counted once per row. |
| `INCLUDE_DESC_KEYWORDS_IN_UNION` | `0`; keyword IDs are still written for audit. |
| `DESC_KEYWORD_REQUIRE_STRUCTURAL_HINT`, `DESC_KEYWORD_STRUCTURAL_HINT_MODE` | `1`, `ANY`, relevant when description expansion is enabled. Default hint sources are relaxed NBD, Pfam NBD, and Pfam enrichment. |
| `UNION_INTERSECT_WITH_PROTEOME_IDS` | `1`; raw IDs, filtered IDs, missing IDs, and counts are retained. |

For each PFAM name, the foreground and background rates are
`(count + 1) / (number_of_rows + 2)`; enrichment is `log2(foreground_rate / background_rate)`.
Foreground rows are deduplicated by protein ID, while the background denominator
is the master-table row count. Extra annotation rows can therefore affect the
background. These are descriptive enrichment heuristics, not significance tests.

`STRICT` raises the enrichment log2 threshold to 1.5 and leaves description expansion off. It still uses relaxed NBD/Pfam discovery; its name does not mean “strict NBD IDs only.” `BROAD` lowers the enrichment threshold to 0.5, enables description expansion, and removes its structural-hint requirement. Explicit configuration values override these presets.

The keyword expression is broad, including terms such as `HET` without added word boundaries. Domain grouping is name based. Changing these expressions to improve specificity would be a separate scientific revision.

### Architecture and confidence

With ILLUMINA defaults, general Pfam/custom hits need domain i-Evalue ≤ `1e-3` and alignment length ≥ 20 aa. `MAX_DOM_I_EVAL=0` means use `EVAL_PFAM`; a nonzero value overrides that Stage 2 threshold. These filters use alignment coordinates and HMMER's reported independent domain E-values.

`NBD_MODE=strict` selects the Stage 1 strict membership list. Any selected Stage 1 member receives `HIGH` confidence. Otherwise, a conservative NACHT/NB-ARC hit receives `HIGH` at i-Evalue ≤ `1e-6` and alignment-length/HMM-length coverage ≥ 0.45, `MEDIUM` at ≤ `1e-3` and coverage ≥ 0.25, and `LOW` otherwise. This coverage is an alignment-length/HMM-length proxy, not a measure of complete-protein recovery. STAND-like or NBD-like names alone are tracked separately and do not establish the conservative NBD call.

`NBD_MODE=relaxed` selects relaxed Stage 1 membership and therefore also makes that selected membership `HIGH`. The interval input remains the strict `nbd_hits.tsv`; changing membership mode does not switch the interval table. `legacy` uses the additive strict-list alias. These details matter when comparing confidence or orientation across configurations.

Recognized SSFR sensor groups remain WD40, ANK, TPR, HEAT, and KELCH. Non-SSFR sensor groups remain C2H2_ZF, ZZ, SPRY, RING, and KINASE. LRR-containing architectures follow an explicit exclusion policy in these rules; this is not a universal claim about fungal NLR biology. A sensor is considered C-terminal when at least one sensor starts at or beyond NBD end + 10 aa; other upstream or overlapping sensor hits do not invalidate that valid C-terminal hit. Effector detection distinguishes any recognized effector from one ending at or before NBD start − 10 aa. An overlap of at least 20 aa between a flagged false-positive domain and the chosen NBD interval triggers the overlap flag. Stage 1 intervals take precedence over domain-derived NBD intervals. These summary intervals can be the bounding range across multiple hits. They are not guaranteed to represent one complete NBD. The per-hit evidence export preserves individual alignment coordinates; see [INPUTS_AND_OUTPUTS.md](INPUTS_AND_OUTPUTS.md#evidence-exports).

ASM scanning is enabled by default and actually scans only the first 200 aa. Accepted ASM hits require i-Evalue ≤ `1e-3`, score ≥ 20 bits, and alignment length ≥ 12 aa. The summary records the lowest-i-Evalue hit per protein; tied hits retain the first input hit. ASM supports flags and rescue prioritization and is not a sensor. The combined custom effector/sensor library is scanned separately; supplying separate provenance libraries does not build that combined database automatically.

Tiering uses `nbd_ok` for HIGH/MEDIUM confidence and uses an expanded tier vocabulary. This includes non-SSFR canonical candidates, low-priority fragments, integrated decoys, LRR/orientation exclusions, and the HIFI suspected-pseudogene category. The presence of a named tier in code is not evidence that a particular run can reach it: with the Stage 2 SSFR-only repeat-region flag, Tier 2D's “repeat-like region without a recognized sensor” condition is unreachable. Rescue scores remain heuristics clipped to 0–100, not calibrated probabilities.

## HIFI defaults

Generate the intended configuration with `funlr init-config --profile HIFI`. The HIFI preset supplies the following differences; explicit YAML/JSON values take precedence. A complete ILLUMINA configuration has explicit thresholds and will retain them if only its profile field is changed later.

| Setting | ILLUMINA | HIFI |
| --- | --- | --- |
| `EVAL_PFAM` | `1e-3` | `5e-4` |
| `MIN_DOM_ALI_LEN` | 20 aa | 25 aa |
| `HARD_END_BP` / alias `CONTIG_END_BP` | 10,000 bp | 1,000 bp |
| `SUBTELO_BP` | 0 (off) | 50,000 bp |
| `REPEAT_PROX_BP` | 10,000 bp | 50,000 bp |
| `NBD_CONF_EV_HIGH`, `NBD_CONF_COV_HIGH` | `1e-6`, 0.45 | `1e-8`, 0.50 |
| `NBD_CONF_EV_MED`, `NBD_CONF_COV_MED` | `1e-3`, 0.25 | `1e-5`, 0.30 |
| `EXON_MAX_INTRON_RETRY` (recorded only; inactive) | 0 | 10,000 bp |

The Python configuration applies profile presets before explicit overrides. Selected Stage 1 membership receives HIGH confidence under both profiles; fallback domain evidence uses the profile-specific confidence thresholds.

Tier and rescue-score behavior depends on the profile: for example, isolated NBD-only candidates can be Tier 2B under ILLUMINA and Tier 4G under HIFI. Profiles encode assembly-context assumptions. They are not inferred from filenames, and their thresholds have not been evaluated for performance across broader inputs.

## Empty evidence and deterministic output

- **Empty evidence and uninitialized state:** zero-hit NBD/Pfam/ASM parses and zero-candidate tables retain defined schemas. Optional empty eggNOG evidence produces empty description IDs. Per-protein effector flags reset before each candidate, so a no-hit protein cannot inherit the previous protein's effector state.
- **Header and boolean handling:** nonempty eggNOG files without `#query` fail clearly. Missing domain cells and missing/false boolean values are handled as absent evidence, rather than accidentally treating a missing value as the literal domain `nan` or the text `False` as true. Explicit textual domain tokens follow the tier helper's domain policy.
- **Explicit LOW confidence:** Stage 3 honors a supplied LOW value. It only infers HIGH/MEDIUM from older evidence columns when the confidence column is absent; a Pfam NBD flag alone does not promote an explicit LOW call.
- **Structural-hint intersection:** `ALL` includes every requested source, including valid empty lists. An empty requested source makes the intersection empty instead of silently weakening the requirement.
- **Stable selection ties:** Pfam enrichment sorts by decreasing log2 fold change, decreasing foreground count, then PFAM name. Priority rescue sorts by decreasing score, then protein ID. This makes top-N and refinement-cap selection reproducible under tied scores.
- **Checked exports:** SeqKit extraction must succeed and return exactly the requested first-token IDs, without missing, unexpected, or duplicate records. An expected empty ID list produces an empty FASTA; a missing upstream table or failed export is an error.
- **Enabled ASM is explicit:** the portable default requires its database. Users must disable ASM explicitly to omit it, instead of silently losing an enabled evidence channel when a profile path is unavailable.
- **Private files and recorded execution:** samtools indexes the private genome copy; HMMER indexes private database copies. User inputs remain untouched. Effective settings, file/tool/code identities, commands, and stage output hashes support verified resume. These execution changes do not add biological evidence.

Two tool-failure fallbacks are explicit: if both Pfam NACHT/NB-ARC models cannot be extracted, Stage 1 logs that Pfam NBD evidence is empty; if an optional whole-proteome Pfam background scan fails, Stage 2 logs its fallback to the candidate background. Other checked scan/export failures stop the run. Inspect the stage logs and `run_info.txt` when comparing results.

Rescue scoring uses substring matching: `VERY_SHORT_PROTEIN` also matches the `SHORT_PROTEIN` score term, giving both +10 and −40 before other terms. These terms are heuristic, not calibrated evidence weights.

## Recorded settings without implemented scientific effects

The configuration records the controls below, but the effective Stage 0–4 code does not implement the named effects. Non-default values are rejected. Their presence in resolved configuration does not mean the described analysis was performed.

| Setting and retained default | Actual behavior |
| --- | --- |
| `ASM_WINDOW_NBD_FLANK_AA=100` | No NBD-flank windows are scanned; only the N-terminal window is implemented. |
| `ASM_ALLOW_MULTIHIT=1` | The ASM summary always uses the single best i-Evalue hit per protein. Filtered hit rows remain in the domain-hit table. |
| `ASM_LOW_COMPLEXITY_MASK=1` | No extra low-complexity masking step is performed by the ASM scanner. |
| `ASM_COUNTS_AS_SENSOR=0` | ASM does not contribute to `has_sensor`; changing that role is unsupported. |
| `ENABLE_REPEATLIKE_SENSOR=1` | No additional repeat-like sensor channel is activated. |
| `REPEATLIKE_GROUPS_CSV=WD40,ANK,TPR,HEAT,KELCH,COILED_COIL,LOW_COMPLEXITY,DUF,RICH_REGION` | `has_any_repeat_region` is calculated from the existing SSFR groups only. The additional named groups are not added. |
| `TIER_REPEATLIKE_UNKNOWN_SENSOR=TIER_2D_PROBABLE_UNKNOWN_SENSOR_REPEATLIKE` | The tier label is fixed; this setting does not rename it. |
| `HIC_IDENTITY_THRESHOLD=85` | No Stage 3 HIC identity detector is implemented. |
| `REPEAT_COLLAPSE_SCORE_BONUS=20` | No additional repeat-collapse score term is applied. |

## Validation and interpretation

Targeted tests exercise empty evidence, ID reconciliation, thresholds, orientation,
ASM filtering, confidence, tier rules, literal exports and failure/resume behavior.
Tool doubles establish command contracts; they do not perform sequence searches
or gene prediction. The [validation record](VALIDATION.md) separates fixture
checks, native-tool execution and comparisons of selected scientific artifacts.
A passing software test or agreement with a reference run does not establish
sensitivity, specificity, immune function or generality across fungal lineages.
Use an independently curated evaluation set for those questions; see
[BENCHMARKS.md](BENCHMARKS.md) and [ROADMAP.md](ROADMAP.md).
