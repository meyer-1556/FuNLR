# Public-sequence hybrid installation fixture

This small fixture runs the existing FuNLR pipeline on public fungal protein
sequences and authentic Pfam profiles, with a real intron-bearing locus embedded
in artificial contigs. It is an installation and regression test. It is not a
natural genome, a validated gene annotation, or a production-model equivalence dataset.

## Run the bundled fixture

With FuNLR and its external tools installed, run from any directory:

```bash
funlr demo --dataset public-sequences --outdir public-sequence-demo
```

Choose a new output directory. The default `funlr demo --outdir synthetic-demo`
continues to run the deterministic synthetic installation test, which also checks
a positive Exonerate alignment. Neither demo requires a database download.

To inspect or run these source files directly, use the supplied configuration:

```bash
funlr validate --config /path/to/FuNLR/examples/demo/config.yaml
funlr run --config /path/to/FuNLR/examples/demo/config.yaml
```

Configuration paths resolve relative to `config.yaml`, including its `output`
directory. To place results elsewhere, add `--outdir /path/to/new/results` to the
run command. Use `--resume` only to verify and reuse the same completed run. HMM
and genome indexing occurs on working copies inside the output directory.

## What is biological and what is constructed

| Component | Content and limitation |
| --- | --- |
| `genome.fa` | The 5,413-bp *Podospora anserina* locus [AF323583.1](https://www.ncbi.nlm.nih.gov/nuccore/AF323583.1) appears twice, surrounded by random DNA. A third contig is entirely random DNA. The native 49-bp intron is preserved. |
| `proteins.faa` | Seven FASTA records match public NCBI/UniProt sequences exactly. They contain six distinct amino-acid sequences: AAL37299.1 and Q8X1P4 identify the same protein. |
| `annotation.gff3` | Handwritten gene/mRNA spans. Two spans correspond to the relocated native locus; the other five are nominal intervals on random DNA. No CDS or exon features are supplied. |
| `emapper.annotations` | Handwritten eggNOG-mapper-style rows with manually adapted protein names. eggNOG-mapper was not run; numeric zeros and dashes are placeholders. |
| `db/nbd.hmm` | Pfam NACHT **PF05729.19** and NB-ARC **PF00931.29** only. These are not the production fungal `Sordariales_ND-domains.hmm` library. |
| `db/pfam_mini.hmm` | A selected 16-family Pfam subset for this fixture, including the same two NBD profiles. It is not a complete Pfam database. |

`het_e_locus` is 29,413 bp (12,000 random bases on each side of the native
insert); `het_e_edge` is 20,413 bp (12,000 upstream and 3,000 downstream);
`controls` is 40,000 random bases. The generator used random seed 7. All three
contigs were reproduced exactly during the provenance audit without executing
the supplied generator.

Every protein sequence and all 16 HMM records were checked against official
source downloads on **2026-09-09**. The packaged versions and SHA-256 hashes,
source URLs, transformations, and verification boundaries are recorded in
[PROVENANCE.json](PROVENANCE.json). Source endpoints may change; the recorded
accession versions and hashes identify this specific fixture.

## Observed behavior in the merged implementation

A real-tool run of the pipeline produced six candidates:

| Demo ID | Tier |
| --- | --- |
| `HETE_AAL37299` | `TIER_1A_HIGH_CONFIDENCE` |
| `HETE_Q8X1P4` | `TIER_1B_NEEDS_REVIEW` |
| `HETE_A7IQW3` | `TIER_1A_HIGH_CONFIDENCE` |
| `KELCH_A0A090DCW7` | `TIER_4A_REPEAT_ONLY_NO_NBD` |
| `TPR_B2B7X0` | `TIER_4C_NO_NBD_NO_REPEAT_LOW_SIGNAL` |
| `ATG1_Q3ZDQ4` | `TIER_4A_REPEAT_ONLY_NO_NBD` |

The three HET-E records enter the strict set. `HETE_Q8X1P4` receives the
`HARD_END` flag because its assigned span is 3,483 bases from the contig end.
The TPR protein has no retained Stage 2 domains: all its `TPR_1` and `TPR_2`
matches exceed the default domain i-E-value cutoff of `1e-3` (the best is
`0.0015`). The updated grouping function recognizes both names as TPR; filtering
occurs before grouping. The updated discovery rules include ATG1 through
its kinase domain; it carries `NON_NLR_ANNOT` and does not enter the strict set.
ORC2 is absent. This small selected example does not measure specificity.

Miniprot maps `HETE_Q8X1P4` across the native intron in both duplicated loci.
The original per-protein summary still aggregates multiple loci using the first
scaffold and minimum/maximum coordinates, so inspect raw GFF alignments when
interpreting these mappings. Feature rows without `Target` are excluded from
query summaries. Rank, identity and positive-match statistics are read from
their native miniprot attributes.

Mapping to multiple scaffolds now triggers refinement. Stage 6 processes one
priority query (`HETE_Q8X1P4`) and three comprehensive queries (all HET-E records),
and checks native Exonerate mappings for each. Either duplicated locus may be
first in the alignment summary. The default synthetic demo also retains its
separate **positive Exonerate identity check**. High-confidence tiering of
`HETE_A7IQW3`, whose assigned genomic interval
is random DNA, illustrates that this workflow classifies the supplied protein
architecture and annotations; it does not validate every supplied gene model.

The fixture configuration disables ASM scanning, fusion hypotheses and the
separate Pfam NBD extraction scan, and uses BROAD discovery with its bundled
miniature databases. Plotting and annotation tagging remain enabled. Protein
translation is disabled because the handwritten GFF3 has no CDS features;
adding invented CDS models would misrepresent this example. The same settings
apply to the one-command demo and manual configuration. See
`docs/RESCUE.md` in a source
checkout for the updated rescue and reporting rules.

## Sources and reuse

The fixture includes source provenance and portable configuration. Read [DATA_LICENSES.md](DATA_LICENSES.md) for UniProt
attribution, Pfam CC0 terms, NCBI molecular-data policy, source accessions, and
modifications. The software license does not replace these data terms.
