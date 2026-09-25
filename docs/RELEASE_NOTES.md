# FuNLR 0.5.0a2 — portable alpha

FuNLR is a **rule-based, Python-orchestrated alpha implementation of a fungal NLR candidate-discovery and targeted gene-model-review workflow**, adapted from an HPC-oriented analysis pipeline. It runs as a local application within a workstation or compute-node allocation. It does not require a particular scheduler or hosted service.

## Available functionality

- Annotated-genome input: matching genome FASTA, protein FASTA and GFF3, with optional annotation tables.
- Profile-HMM discovery, domain architectures, evidence tiers, review flags and priority/comprehensive genome-rescue tracks using miniprot and Exonerate.
- Candidate and strict-set tables, FASTA/BED exports, annotation tagging, domain plots and offline HTML reports.
- TSV-driven sequential batch runs, input QC, model inventories and checksum-bound reference-panel comparisons.
- Candidate/exclusion evidence and original per-hit NBD coordinates, scores and aligned sequence segments.
- Recorded settings, commands, input/tool identities, output checksums and verified resume.
- Packaged synthetic and public-sequence demonstrations for installation and regression checks.

The current analysis uses explicit domain, sequence and genomic-context rules.

## Installation and data

Download and extract the source archive, then follow the [versioned installation guide](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/docs/INSTALLATION.md). The Python package needs external bioinformatics executables; a pip install alone does not supply them. Dependency locks and an installer are included, with platform-specific validation limits documented.

The bundled demonstrations run without the production fungal-model bundle. **Production fungal HMM libraries must be supplied separately**; a public download bundle and its redistribution documentation are still pending. Pfam acquisition and database configuration are described in the [database guide](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/docs/DATABASES.md). Public Conda-channel packages and published container images are not yet available.

## What this alpha establishes

Software tests, installation demonstrations, reference-output checks and independent biological accuracy are distinct forms of evidence. The [validation record](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/docs/VALIDATION.md) describes the checks performed, their environments and their limits. The demonstrations are constructed regression fixtures, not independent fungal benchmarks. The reference-panel comparison command provides evaluation infrastructure; it does not itself establish model accuracy or independence from training data.

Candidate tiers and rescue scores are heuristic evidence categories. They do not prove biological function or the correctness of a reconstructed gene model. Fusion hypotheses require review. Independent biological accuracy has not been established. The [validation record](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/docs/VALIDATION.md) and installation guide describe the current platform and testing limits.

This is an **alpha prerelease**: interfaces, output schemas and scientific policies may change as validation develops. Retain inputs, configuration, software/database versions and intermediate evidence when using it for research. Report reproducible installation or analysis problems through [GitHub Issues](https://github.com/meyer-1556/FuNLR/issues), omitting sensitive paths and data.

Original FuNLR code is MIT licensed; incorporated components and bundled data retain their own notices. See [LICENSE](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/LICENSE), [NOTICE](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/NOTICE) and [CITATION.cff](https://github.com/meyer-1556/FuNLR/blob/v0.5.0a2/CITATION.cff). Cite the external tools and profile resources used in each analysis as well.
