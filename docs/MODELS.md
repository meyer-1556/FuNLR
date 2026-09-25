# Model identities and source evidence

FuNLR uses a fixed reference model snapshot. Production fungal HMMs are not
bundled with the package: complete redistribution permissions and a reproducible
rebuild from archived training sources have not been established. The bundled
demo models exercise installation and must not replace the production libraries.

The snapshot has five libraries, with these profile-record counts:

| Library | Records | Pipeline role | Extra records identical except NAME |
| --- | ---: | --- | ---: |
| NLR_NBD_combined.hmm | 41 | Nucleotide-binding domain discovery | 0 |
| NLR_custom_combined.hmm | 127 | Combined effector/sensor architecture | 30 |
| NLR_effectors_combined.hmm | 19 | Effector component | 1 |
| NLR_sensors_combined.hmm | 108 | Sensor component | 29 |
| NLR_ASM_combined.hmm | 39 | Amyloid signaling motif annotation | 0 |

These are profile records, not counts of independently validated NLR-specific
families. The custom library is the exact byte concatenation of the effector
library followed by the sensor library. Consequently, the table also counts
component records again inside the combined library.

## Inspect a model directory

```sh
funlr model-inventory --directory /path/to/models
funlr model-inventory --directory /path/to/models --verify-snapshot
```

Both commands print JSON to standard output and leave the directory unchanged.
The first inventories direct `*.hmm` files; subdirectories, compressed files and
pressed HMMER sidecars are outside its scope. Library symlinks are rejected.
The second additionally invokes the existing `verify-models` check: all five
fixed libraries must match their recorded sizes, profile counts and SHA256
checksums, and the custom concatenation must match. Additional direct HMM files
are inventoried but are not part of that five-library verification.

Inventory schema version 1 records each profile in file order, including its
`name`, nullable `accession`, model `length`, nullable `nseq`, byte count and
SHA256. `nseq` reports the HMM header, not an independently recovered training
sequence list. A second hash omits only the NAME line; duplicate groups are
reported without deleting or merging any models. Other header and parameter
lines remain in that comparison. This is a header/record-structure check, not
a substitute for HMMER parameter validation.

Known roles, citations and source-audit descriptions are attached only when
both filename and whole-library SHA256 match the packaged metadata. A familiar
filename with changed bytes has `snapshot_identity: UNRECOGNIZED`, with null
role, citations and source audit. A null value means unknown, not absent or
disallowed. Profile-level publication identity, construction provenance and
redistribution permission remain null where the evidence does not establish
them. A source citation is not a license grant or a specificity claim.

The checked reference inventory is [MODEL_INVENTORY.json](MODEL_INVENTORY.json).
It contains metadata and hashes, not HMM parameter matrices. The versioned
library metadata ships in `funlr/data/model_sources.json`; snapshot identities
ship in `funlr/data/model_snapshot.json`. The supporting audit is
[MODEL_SOURCE_AUDIT.json](MODEL_SOURCE_AUDIT.json).

## What the source audit establishes

The NBD library combines reported Sordariales NBD profiles, six sequence-derived
models built with MAFFT and hmmbuild, and additional locally rebuilt models from
a supplemental sequence table. The effector library combines published and
curated terminal-domain resources; the sensor library uses curated non-NB
terminal-domain selections. The ASM library combines three supplemental HMM
resources. These origins are described in
[Bonometti et al.](https://doi.org/10.1371/journal.pgen.1011739),
[Wojciechowski et al.](https://doi.org/10.1371/journal.pcbi.1010787), and
[Dyrka et al.](https://doi.org/10.1093/gbe/evu251), with the precise reported
component relationships and evidence limits retained in the source audit.

Component bytes agree with the supplied manifest. The complete historical
rebuild is still unknown: locally rebuilt profiles are not automatically
identical to published models, the sensor extraction and final record counts
differ, and some duplicate records differ only in name. The original custom
assembly command is also unknown even though its final byte composition is
verified. Current releases preserve those bytes and defaults.

A future model release should archive original assets and checksums, selection
lists, alignments, tool versions, build parameters, immutable per-profile source
identities and redistribution terms. It should explicitly map old profiles to
new profiles and benchmark changed calls. Removing duplicate models or replacing
HMMs changes a scientific database release and requires that comparison.

General domain annotation uses the verified [Pfam38.2 release](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/).
See [DATABASES.md](DATABASES.md) for download and configuration instructions.
