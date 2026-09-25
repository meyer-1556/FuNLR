# Exact dependency locks

`conda-linux-64.lock` and `conda-osx-arm64.lock` use Conda's `@EXPLICIT` format: each line identifies a concrete package build URL and its MD5 checksum. `locks.json` additionally records SHA-256, package size, version/build, dependencies, solver, and source specification. The installer checks that the selected lock and manifest agree before using them. The pinned bootstrap is independently SHA-256 verified.

These files resolve all runtime dependencies, not FuNLR itself or any biological database. `install.py` installs the source version after creating the environment. Locks avoid fresh solving of floating transitive dependencies at installation time. They do not guarantee that upstream artifacts remain hosted forever, that every ABI/hardware combination works, or that two runs using different profile databases are scientifically equivalent. Archive validated release environments/images and reference data alongside the source release.

The active locks were generated using micromamba 2.9.0 from
`environment.yml` and `environments/environment-osx-arm64.yml` on
2026-09-10 UTC. Matplotlib 3.10.8 replaces the R dependency stack in the default
installation:

| Active lock | Packages | Compressed dependency downloads |
| --- | ---: | ---: |
| `conda-linux-64.lock` | 103 | Approximately 149.3 MiB |
| `conda-osx-arm64.lock` | 97 | Approximately 112.5 MiB |

These totals exclude the bootstrap, installed files, caches, FuNLR source and
biological databases. Both locks resolved successfully. Resolution does not
execute the resulting binaries or validate a fresh installation. The Linux
solve used virtual Linux/glibc overrides because it ran on macOS. The advertised
Linux floor remains glibc >=2.28 pending native validation even though the new
solution's virtual dependencies allow glibc 2.17. Apple Silicon requires macOS
>=11. See [the validation record](../docs/VALIDATION.md) for actual execution
evidence and the local installation restriction.

## Optional R rendering

The default lock includes Matplotlib. The optional `reporting.plot_backend: r`
renderer requires a separately supplied compatible R environment; see
[installation](../docs/INSTALLATION.md). Validate and record that custom runtime.

## Regenerating active locks

A maintainer updating dependencies should solve environment.yml on each matching target, retain the solver JSON, regenerate the explicit URLs/checksums and locks.json from the same solution, and review changed packages. For example, the read-only solving step is:

```bash
micromamba --no-rc create --dry-run --platform linux-64 \
  --prefix /tmp/funlr-lock-plan --file environment.yml --json > linux-solve.json
```

Review the changed package builds. Keep direct versions aligned with
pyproject.toml, doctor, the Conda recipe, and documentation. Run fresh
installation, doctor, demos, regression tests, scientific regression comparisons, and both
container checks before declaring the new distribution validated. The current
release record identifies checks still pending. Checksums identify content;
license and citation obligations remain with the corresponding dependencies.
