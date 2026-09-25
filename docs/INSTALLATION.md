# Installation and supported platforms

Download the FuNLR source ZIP, extract it, and run these commands in the directory containing install.py:

```bash
python3 install.py --prefix "$PWD/funlr-env" --demo
export PATH="$PWD/funlr-env/bin:$PATH"
funlr doctor --demo
```

The bootstrap needs Python 3.9 or newer, internet access, and write access to the destination and its parent. Micromamba also uses its standard `~/.conda/environments.txt` registry, so the home-directory registry must be writable. The installed application has its own Python 3.12.14. Root/admin privileges, git, an existing Conda installation, and changes to shell startup files are not required. Prefixes must be new and outside `src/`; an existing prefix is never overwritten. Spaces in normal paths are supported. Control characters are rejected.

The installer selects a platform-specific `@EXPLICIT` lock, downloads micromamba 2.9.0 over HTTPS, verifies its SHA-256, installs the locked packages, builds and installs FuNLR without fetching additional Python dependencies, runs doctor, and optionally runs the demo. It ignores inherited pip/Conda/mamba settings that could redirect the installation. It records `funlr-installation.json` inside the prefix. Package downloads are cached in `.funlr-bootstrap` beside the destination. The registry and cache are the only intended locations outside the chosen environment and optional sibling demo output. Temporary source staging under that cache is removed on success.

```bash
# Inspect the exact plan without downloads, execution, or writes.
python3 install.py --prefix /shared/software/funlr-0.5.0a1 --dry-run
# Use an already installed, exact micromamba 2.9.0 binary.
python3 install.py --prefix /shared/software/funlr-0.5.0a1 --micromamba /path/to/micromamba
```

A failed install retains partial files for diagnosis and reports a nonzero exit. Choose a new prefix for a retry, or inspect and remove only the failed prefix yourself. If a cached bootstrap checksum fails, remove the named corrupt archive and retry. Do not disable certificate or checksum verification. Normal network proxies/CA settings may be needed on your site.

## Platform contract and validation status

| Platform | Installation path | Status of this archive |
| --- | --- | --- |
| Linux x86_64, glibc >=2.28 | Exact Conda lock; Linux Docker/Apptainer image | Exact dependency lock supplied; fresh installation and container checks require recorded results in [VALIDATION.md](VALIDATION.md) |
| macOS Apple Silicon (arm64), macOS >=11 | Exact Conda lock | Exact dependency lock supplied. A complete new-prefix installation remains unverified; native runtime checks and limitations are reported separately in [VALIDATION.md](VALIDATION.md). |
| Windows | Linux x86_64 environment through WSL2, or Linux container | No native Windows package; not tested here |
| Intel macOS | Linux x86_64 container/VM | No native dependency lock in this archive; not tested here |
| Linux ARM, Power, other CPU architectures | Requires a separately solved and tested dependency set | Not claimed supported |

The application has no fixed cluster paths, required module names or scheduler API. It does not mean every CPU architecture, Linux distribution, filesystem, and HPC policy has been certified. The first release targets Linux x86_64; expand this table only after testing. See [Micromamba's guide](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for the environment manager's prefix behavior.

## What is installed

Runtime tools are HMMER 3.4 (hmmscan/hmmpress/hmmfetch and hmmbuild for the demo), samtools 1.21 on Linux (the retained 1.24 pin in the Apple Silicon lock), SeqKit 2.12.0, miniprot 0.18, Exonerate 2.4.0, and gffread 0.12.7. Python dependencies are pandas 2.2.3, NumPy 1.26.4, PyYAML 6.0.3, and Matplotlib 3.10.8. The default plotting backend is headless Matplotlib; R and its packages are not required by the default installation. Full package builds, transitive dependencies, URLs, and checksums are in [environments/](../environments/README.md). Direct-version environment.yml remains for inspection and future lock regeneration; installing it through a fresh solve can choose different builds.

The current locks contain 103 packages / approximately 149.3 MiB of compressed downloads on Linux, and 97 packages / approximately 112.5 MiB on Apple Silicon, plus a small bootstrap. Installed environments and package caches require more space. These sizes exclude Pfam, fungal models, biological inputs, and result files; production resource requirements must be measured. Neither demo needs production databases. The installed package also supports `funlr demo --dataset public-sequences --outdir /path/to/new-public-demo`; that checked public-sequence hybrid fixture does not replace the production fungal profiles.

## Optional R rendering

An optional R renderer is available for the domain figures.
Provide a compatible R 4.4 installation and its dplyr, ggplot2, RColorBrewer,
colorspace, tidyr, and scales packages, then change the existing reporting
section in your configuration:

```yaml
reporting:
  plot_backend: r
```

Check that environment with `funlr doctor --plot-backend r`. The `Rscript`
executable can be selected in the configuration's `tools.Rscript` setting.
The run records the selected backend and checks the requested plot outputs.
Switching backends changes the run fingerprint, so use a new output directory.

## Other installation paths

For an existing compatible dependency environment:

```bash
python -m pip install .
funlr doctor --demo
funlr demo --outdir /path/to/new-demo
```

`pip` does not install the external bioinformatics tools. The wheel is usable with Python 3.11–3.12 and the declared Python dependencies; doctor checks the tools on PATH. Avoid mixing incompatible site modules and an environment. A binary can exist yet fail to start because a shared library is missing; doctor reports that error.

The local Conda build recipe is under packaging/conda/. It is a recipe for maintainers, not an already published `conda install funlr` channel package. Container recipes and build instructions are in [CONTAINERS.md](CONTAINERS.md). No image has been pushed to Quay or another registry.

## Offline compute nodes

Install on an allowed networked node into a shared prefix with the same platform, then leave the environment at that absolute path. Download biological inputs/databases separately, record their releases and SHA-256 values, and copy them to the cluster before submitting a run. FuNLR discovery and the demo need no network once the software is installed. For transfer between machines, prefer a tested Apptainer image; do not assume copying a Conda directory relocates its embedded paths correctly. A complete offline bootstrap bundle is not included in this archive.

## Exonerate runtime limitation

The stock Exonerate 2.4.0 binary on a tested Apple Silicon host intermittently aborted on full-reference queries. FuNLR rejects its nonzero exit rather than accepting partial output. An optional [verified-source serial build](../packaging/exonerate-serial/README.md) is supplied for maintainers; it requires a compiler, GLib headers and build tools. Its small compatibility patch fixes compilation of the upstream serial mode without changing alignment rules. It is not automatically installed by the locked installer. A successful run using that separate executable does not certify the lock's unmodified Exonerate binary.

Use the validation record to distinguish this tested workaround from the still-pending Linux/container and lock-only installation checks. This remains an alpha until each advertised production installation route passes the real reference, not just the small demos.

## Combined installation self-test

```bash
funlr doctor --self-test
funlr doctor --self-test --config settings.yaml --models-dir /path/to/combined_hmms
```

The first command checks native tool version commands and executable hashes, the
pinned Python dependency versions, HMMER's demo-build tool and all bundled public
fixture checksums. It works from an installed wheel as well as a source checkout.
Missing or damaged bundled resources fail the check. No searches, plots or output
directories are created. The native demonstration remains the execution check.

`--config` uses the same executable paths and plotting backend as the pipeline.
`--models-dir` additionally verifies the five fixed reference libraries; without it,
production models are labeled `NOT_REQUESTED`. Pfam and the genome/annotation
inputs are not validated by this installation check. Run plain `funlr validate
--config settings.yaml` for the scientific input contract. `funlr validate
--self-test` is a compatibility alias for the combined installation check.
