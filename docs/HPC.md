# Running FuNLR on an HPC system

FuNLR runs one annotated genome on one compute node. The Python application does not submit Slurm jobs, load environment modules, or download reference data during `run`. Your scheduler launches the installed command. `--threads` controls supported native tools within that node; it does not distribute stages across nodes. Exonerate queries are processed serially in the priority and comprehensive tracks.

## Prepare the environment and configuration

Install on a node that permits network access and environment creation, using a prefix visible from the compute nodes. The compute node must have a compatible OS and architecture; an installed Conda prefix is not assumed relocatable. Check the environment from an allocation with `funlr doctor --demo`, which also checks hmmbuild for the installation demonstration. Default production output uses gffread for integration and headless Matplotlib for plots. R and its packages are required only when `reporting.plot_backend: r` is selected. See [INSTALLATION.md](INSTALLATION.md).

Use your site's supported Apptainer setup if a container is preferable. Build or obtain and test the image before submitting production work, and bind the input/database/output locations as required by that site. The repository's container instructions are in [CONTAINERS.md](CONTAINERS.md); cluster-specific container policies still apply. Scheduler and container combinations have not been universally validated.

Generate the intended profile first:

```bash
funlr init-config --profile ILLUMINA > settings.yaml
# Edit sample, input, database, and output paths in settings.yaml.
funlr validate --config settings.yaml
funlr run --config settings.yaml --dry-run
```

Choose HIFI at configuration generation for the corresponding long-contig assumptions. Explicit values in a full configuration override profile presets, so changing just its profile later does not reset all thresholds. The defaults are BALANCED discovery, enabled ASM, both rescue tracks, enabled fusion, plots, and annotation integration. Provide the ASM and combined custom profiles used for the intended analysis, or explicitly disable/omit channels and record that change. See [INPUTS_AND_OUTPUTS.md](INPUTS_AND_OUTPUTS.md) and [METHODS.md](METHODS.md).

## Submit one run

The example [run.slurm](../examples/hpc/run.slurm) requests one node, one task, eight CPUs, 16 GB RAM, and four hours. These illustrate submission syntax; they are not measured production requirements. Edit the account, partition, memory, time, and other site-specific settings before submission:

```bash
export FUNLR_PREFIX=/shared/software/funlr-env
export FUNLR_CONFIG=/shared/project/sample/settings.yaml
sbatch examples/hpc/run.slurm
```

The script prepends the installed environment's `bin/` directory to PATH and passes `SLURM_CPUS_PER_TASK` to `funlr run --threads`. Set `execution.output_dir` in the configuration to a new writable run directory. Allocate at least one CPU and keep `--threads` within the job allocation. If recording memory with `--ram-gb` or `execution.ram_gb`, use the allocation value; it is provenance, not a memory cap.

PBS, LSF, SGE, or another scheduler can launch the same ordinary command within a site-specific single-node job script:

```bash
funlr run --config /shared/project/sample/settings.yaml --threads 8
```

Those scheduler-specific submission paths have not each been tested. For several genomes, use `funlr batch` inside one allocation to process samples sequentially, or submit independent jobs with distinct configurations and output directories. See [BATCH.md](BATCH.md). FuNLR does not itself submit jobs or distribute samples across scheduler nodes.

## Storage, resources, and restart

Keep input genomes, GFF3, proteins, annotation tables, and versioned HMM databases at stable readable paths. FuNLR makes private genome/database copies and indexes under the output directory; source files and shared databases are not modified. Budget storage for those copies, pressed databases, raw scans, two alignment tracks, plots, and integration exports. Scratch output can reduce contention if appropriate for your site. Retain the entire run directory, including `work/`, `logs/`, `manifest.json`, `resolved_config.json`, and `state.json`, when archiving or preparing to resume.

Benchmark the actual genome/proteome, database versions, and enabled stages to choose memory and wall time. HMM database size and the number of candidates/refinement queries affect resource use; the presence of the comprehensive track means a run without priority rescue candidates can still invoke aligners. See [VALIDATION.md](VALIDATION.md) for validation and benchmarking evidence. Neither a passing tiny demo nor replaying archived HMMER outputs establishes a production allocation.

After an interrupted run, use the same environment, paths, scientific configuration, sample metadata, and thread/allocation labels:

```bash
funlr run --config /shared/project/sample/settings.yaml --threads 8 --resume
```

Resume verifies input contents and paths, code, tool executables, runtime/library identities, settings, and saved stage artifacts. It skips verified stages and recomputes damaged stages and their dependents. A different node is usable only when these identities remain compatible. Moving the environment, changing tools or settings, or relocating the output directory can require a fresh run; do not edit manifests to force reuse.

Normal SIGINT/SIGTERM handling stops the stage's external processes and removes the lock. Hard kills, node loss, or filesystem failures can leave `.run.lock`; first confirm that its recorded process and job have stopped, then remove the stale lock before resuming. Inspect the scheduler output and `logs/` to distinguish a successful empty result from a failed search.

Only installation and reference acquisition need network access. A configured production run uses local files and does not contact a remote service. Model availability and redistribution rights are separate from installation; consult [DATABASES.md](DATABASES.md) before distributing a database or container containing it.
