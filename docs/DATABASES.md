# Databases and model identity

FuNLR requires versioned HMM libraries for the evidence channels selected in its
configuration. The fixed reference snapshot records profile counts, file sizes
and SHA-256 values in [MODEL_SNAPSHOT.json](MODEL_SNAPSHOT.json). These identities
make a particular collection reproducible; they do not establish that every
profile is NLR-specific or that the collection can be rebuilt from original
sequence sources.

| File | Profiles | Use |
| --- | ---: | --- |
| `NLR_NBD_combined.hmm` | 41 | Strict and relaxed NBD discovery |
| `NLR_custom_combined.hmm` | 127 | Custom domain architectures |
| `NLR_effectors_combined.hmm` | 19 | Component of the custom library; optional separate provenance input |
| `NLR_sensors_combined.hmm` | 108 | Component of the custom library; optional separate provenance input |
| `NLR_ASM_combined.hmm` | 39 | Amyloid signaling motif scan |
| `Pfam-A.hmm`, release 38.2 | 30,134 | Pfam NBD and architecture evidence |

Counts refer to HMM profiles, including repeated model parameters under distinct
names. Sensor-like repeat domains occur in other protein classes too. Retain
false-positive and context evidence when interpreting candidates.

## Acquire and check Pfam

```bash
funlr fetch-pfam --outdir /path/to/databases/pfam38.2
```

This downloads the fixed EBI Pfam 38.2 HMM archive and verifies the decompressed
file against the reference SHA-256. A matching existing file is reused; a
different file is not overwritten. The file is approximately 2.25 GB
uncompressed. Allow extra space for download, verification, the pipeline's
private copy and its indexes. `funlr run` does not download reference data.

The source is [EBI Pfam 38.2](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/).
The uncompressed HMM SHA-256 is
`4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`.
The compressed archive's upstream MD5 and release record are retained in
[PFAM38_2_VERIFICATION.json](PFAM38_2_VERIFICATION.json).

## Check and inspect the fungal libraries

```bash
funlr verify-models --directory /path/to/combined_hmms
```

Verification checks the fixed file identities. The model inventory provides
profile-level identity and traceability separately; see [MODELS.md](MODELS.md)
and [MODEL_INVENTORY.json](MODEL_INVENTORY.json). A model that differs from this
snapshot is not necessarily invalid, but it represents a different analysis
input and needs its own provenance and evaluation.

Set `databases.nbd_hmms`, `custom_hmms` and `asm_hmms` to the relevant files.
Separate `effector_hmms` and `sensor_hmms` inputs record component identities;
they do not automatically form a replacement for the scanned `custom_hmms`
library. ASM is enabled by default and therefore requires its library unless
explicitly disabled. Omitting custom-domain or annotation evidence can change
results. See [INPUTS_AND_OUTPUTS.md](INPUTS_AND_OUTPUTS.md).

## Composition of the custom library

The fixed `NLR_custom_combined.hmm` is reproduced byte for byte by concatenating
`NLR_effectors_combined.hmm` followed by `NLR_sensors_combined.hmm`. All 127 model
blocks match those components; no additional NBD or ASM blocks are present.
The combined file SHA-256 is
`d0ebec72a6092276f4cc16229869d5f6787de5f3949c54661d125998d68954e7`.

The collection includes HeLo/HET/Ses/Goodbye-related profiles and WD40, ankyrin,
TPR, HEAT, kelch and other domains. Some names contain `__dup` suffixes. The
snapshot retains those names and all model parameters. The verified composition
is recorded in [MODEL_COMPOSITION.json](MODEL_COMPOSITION.json).

To make a new copy from the exact components:

```bash
python scripts/reconstruct_custom_hmm.py --models-dir /path/to/combined_hmms --output /path/to/new/NLR_custom_combined.hmm
```

The helper checks both component hashes and the final hash and refuses an
existing output. This reconstruction establishes the custom file's composition;
it does not recover the original component-building command or rebuild profiles
from publication sequences.

## Source records and unresolved provenance

The retained build manifest, logs and curation records support the ancestry
summarized in [MODEL_SOURCE_AUDIT.json](MODEL_SOURCE_AUDIT.json). Component files
match the manifest's recorded identities. This connects the records to the
snapshot, without independently reproducing every model from publication data.

| Library | Reported source material |
| --- | --- |
| NBD | Bonometti Sordariales NBD profiles; six Dyrka FASTA-derived models built with MAFFT and `hmmbuild`; models built from Wojciechowski supplementary sequence table `pcbi.1010787.s002` with taxonomic splitting. A reassignment record moves an NBD-like model from the effector collection; `NACHT_sigma` is present in the NBD file. |
| Effectors | Wojciechowski `pcbi.1010787.s008.hmm` and the Bonometti N-terminal extract selected using S15 (`pgen.1011739.s018`), followed by category reassignment and name handling. |
| Sensors | Bonometti non-NB S12 selection (`pgen.1011739.s015`) from the C-terminal/N-terminal library. An optional S013 heuristic appears in the build code but is not recorded in the final manifest. |
| ASM | Wojciechowski `pcbi.1010787.s014.hmm`, `s020.hmm` and `s015.hmm`. |

The source records cite Bonometti et al. (2025), Wojciechowski et al. (2022) and
Dyrka et al. (2014). The Wojciechowski-derived NBD profiles were built from a
sequence table; they are not an unchanged published NBD HMM download. The build
signature records minimum group sizes of 50 sequences per phylum and 25 per
class, with at most 5,000 sequences per taxonomic group.

The sensor build log reports 78 extracted profiles but 108 final profiles.
Inspection finds 28 groups identical after removing only their `NAME` lines,
accounting for 29 extra sensor records. The effector collection adds one such
pair, `SesA` and `SesA__dup2`. Thus the custom library contains 30 records whose
model content duplicates another retained record under a different name.

Build code reuses an intermediate `uniq` directory and concatenates all model
files in it. Retained intermediate files are a possible explanation, not proof
of the original cause. Those intermediate directories are unavailable. Removing
records or rebuilding components may affect search statistics and downstream
evidence and requires a separately evaluated model release. A traceable rebuild
needs immutable source inputs, clean intermediate directories and a profile-level
source map.

## Distribution and permissions

Production fungal HMMs are supplied separately. The Python wheel contains model
metadata and small installation fixtures, not these production databases. MIT
licensing of FuNLR's original code does not settle model redistribution rights;
retained third-party code and demo data keep their own notices. Public model
permissions, citations, immutable hosting and complete build provenance remain
release work described in [ROADMAP.md](ROADMAP.md).

The small synthetic and public-sequence demo libraries exercise installation and
execution. They are not substitutes for research databases or a biological
benchmark. Neither a successful demo nor an exact model hash validates a
candidate's function.
