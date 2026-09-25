# Data sources, attribution, and reuse terms

This notice covers the public-sequence hybrid fixture in this directory.
The FuNLR software license does not relicense the third-party database content.
All source identities below were checked on 2026-09-09. Exact bundled bytes,
sequence hashes, profile accessions and versions, source-response hashes, and
modifications are recorded in [PROVENANCE.json](PROVENANCE.json).

## Pfam / InterPro profiles — CC0 1.0

The 16 HMM records in `db/pfam_mini.hmm`, and the two-record subset in
`db/nbd.hmm`, are unmodified official Pfam HMM records distributed through
InterPro. Every record matched an official download byte-for-byte, excluding
blank lines between records in the concatenated libraries.

Pfam/InterPro downloadable data are supplied under the
[CC0 1.0 Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/),
as stated in the [InterPro data license](https://www.ebi.ac.uk/interpro/about/license/).
Credit: Pfam, maintained by EMBL-EBI and the Xfam consortium, distributed through
InterPro. For a research analysis, cite the applicable Pfam and InterPro database
publications using their [citation guidance](https://www.ebi.ac.uk/interpro/about/).
No Pfam seed alignment was retrained for this fixture.

These profiles are PF00004.36, PF00023.37, PF00069.32, PF00270.36, PF00271.38,
PF00400.39, PF00515.35, PF00735.25, PF00931.29, PF01344.32, PF01580.25,
PF02463.26, PF05729.19, PF07646.22, PF07719.24, and PF12796.14. The NBD subset
contains only PF05729.19 and PF00931.29.

## UniProt protein records and names — CC BY 4.0

Credit: **UniProt Consortium**, UniProtKB. UniProt applies
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
to copyrightable parts of its databases; see the
[UniProt license and disclaimer](https://www.uniprot.org/help/license/).
The following records contribute sequences and manually adapted names:

| Demo ID | Source record | Sequence version |
| --- | --- | --- |
| `HETE_Q8X1P4` | [Q8X1P4](https://www.uniprot.org/uniprotkb/Q8X1P4/entry) | 1 |
| `HETE_A7IQW3` | [A7IQW3](https://www.uniprot.org/uniprotkb/A7IQW3/entry) | 1 |
| `KELCH_A0A090DCW7` | [A0A090DCW7](https://www.uniprot.org/uniprotkb/A0A090DCW7/entry) | 1 |
| `TPR_B2B7X0` | [B2B7X0](https://www.uniprot.org/uniprotkb/B2B7X0/entry) | 1 |
| `ORC2_B2B6T2` | [B2B6T2](https://www.uniprot.org/uniprotkb/B2B6T2/entry) | 1 |
| `ATG1_Q3ZDQ4` | [Q3ZDQ4](https://www.uniprot.org/uniprotkb/Q3ZDQ4/entry) | 1 |

Changes: amino-acid sequences are unchanged, FASTA headers were replaced with
demo IDs, protein names were manually adapted, and the records were combined
with artificial contigs, coordinates, and placeholder annotation fields.
Only Q3ZDQ4 is a reviewed UniProtKB/Swiss-Prot record; the other five are
unreviewed UniProtKB/TrEMBL records. A public sequence match does not establish
experimental validation of each protein or annotation. Q8X1P4 cross-references
AAL37299.1 and has the same sequence; the two demo rows are duplicates, not
independent allele controls.

## NCBI / GenBank nucleotide locus and translation

Credit: Espagne, Balhadere, Penin, Barreau, Turcq, and colleagues for the
*Podospora anserina* het-e sequence submission and associated work;
NCBI/GenBank for distribution. The source records are
[AF323583.1](https://www.ncbi.nlm.nih.gov/nuccore/AF323583.1) and
[AAL37299.1](https://www.ncbi.nlm.nih.gov/protein/AAL37299.1).
Associated publication: Espagne et al. (2002), *HET-E and HET-D belong to a new
subfamily of WD40 proteins involved in vegetative incompatibility specificity
in the fungus Podospora anserina*, Genetics 161:71–81
([PubMed 12019224](https://pubmed.ncbi.nlm.nih.gov/12019224/)).

The [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/#data)
places no NCBI restriction on use or distribution of these molecular database
records. NCBI also explains that contributors may retain rights; this notice
does not assign a new software license to those data.

Changes: the native nucleotide sequence is unchanged within each of two
insertions, flanked by generated random DNA. The native intron is retained;
its CDS annotation is `join(811..3093,3143..4930)`. The AAL37299.1 amino-acid
sequence is unchanged and its header was renamed `HETE_AAL37299`.

## Constructed demonstration materials

Random flanks, nominal coordinates and the minimal eggNOG-style table are
constructed demonstration materials, not experimental measurements or outputs
from eggNOG-mapper. Contributed construction materials retain BSD-3-Clause terms;
see the project `NOTICE` and `src/funlr/THIRD_PARTY_LICENSE.txt` for the retained
copyright notice and conditions. Sequence, annotation and HMM bytes are unchanged
from the reference fixture; configuration and documentation describe their
portable use and scientific limitations.
