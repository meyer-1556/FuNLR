# Maintenance roadmap

The alpha provides annotated-genome NLR candidate discovery, evidence tiers, targeted genome rescue and inspectable outputs. This roadmap covers maintenance, distribution and validation of those existing functions. It does not promise new scientific capabilities.

## Implemented in 0.5

- Descriptive sequence and annotation QC with explicit denominators.
- Candidate membership, review flags and individual accepted NBD HMM-hit exports.
- Profile inventories, library checksums and source metadata that retain unresolved provenance.
- An explicit reference-panel comparison command with checksum-bound inputs and defined evaluation units.
- Recorded QC, evidence, elapsed time and output sizes in status and HTML review.

## Release and maintenance priorities

| Work | Completion criterion |
| --- | --- |
| Production model distribution | Confirm redistribution terms for every library, publish immutable assets with checksums and source identifiers, and document unresolved construction history. Do not silently replace the byte-verified libraries. |
| Installation coverage | Exercise clean Linux x86_64 and Apple Silicon installations, Linux container builds and an Apptainer deployment on an HPC system. Record exact versions and workarounds. |
| Published packages | Publish tested Conda packages and versioned container images; check cold downloads and operation on offline compute nodes. Recipes alone are not distributed software. |
| Validation records | Extend regression cases for documented inputs and edge conditions. Keep software checks, agreement with saved reference outputs and independent biological evaluation separate. Report evaluation inputs, reference availability and limitations. |
| Resource guidance | Measure wall time, CPU allocation, storage and peak memory where available on representative completed runs. Do not infer production requirements from a tiny demo. |
| Output contracts | Document stable column definitions, coordinate conventions, missing-value states and schema versions. Preserve raw evidence and compatibility notes for changed formats. |
| User documentation | Improve installation troubleshooting, configuration examples and failure messages using reported problems. Keep current functionality separate from unimplemented settings. |
| Release archiving | Archive tagged source and distributions, publish checksums and add the archived release identifier to citation metadata. |

Scientific-policy changes require an explicit version, documented rationale and evaluation before becoming defaults. Counts alone do not establish biological gains or losses, and an alignment does not independently validate a proposed gene model.

## Design references

- [EGAP](https://github.com/iPsychonaut/EGAP): sample-sheet organization and explicit orchestration.
- [PlantLRR-PRR](https://github.com/PHYTOPatCAU/RLP_identification): staged evidence inspection and workflow documentation; its plant cell-surface receptor rules are not the fungal classification policy.
- [RefPlantNLR](https://doi.org/10.1371/journal.pbio.3001124): explicit reference collections and evaluation units.
- [DefenseFinder](https://github.com/mdmparis/defense-finder) and [PADLOC](https://github.com/padlocbio/padloc): versioned model resources and inspectable evidence.
- [NLR_Sordariales](https://github.com/bonospora/NLR_Sordariales), [FixingHetDE](https://github.com/SLAment/FixingHetDE) and [MolEvoNLRs](https://github.com/SLAment/MolEvoNLRs): relevant fungal NLR analysis and gene-model work.

These references inform documentation and engineering choices. FuNLR does not claim to reproduce their biology or measured performance.
