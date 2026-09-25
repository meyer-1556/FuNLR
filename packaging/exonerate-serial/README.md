# Optional serial Exonerate build

This is a maintainer workaround for a native Exonerate 2.4.0 memory error found
while testing FuNLR on Apple Silicon. It builds a separate executable; it does
not replace the Conda lock's executable or alter FuNLR's checked exit handling.

The stock binary intermittently aborted or segfaulted on identical valid
protein-to-genome commands. An AddressSanitizer build of the original source
reproduced a heap-use-after-free in `Alphabet_destroy`, called by
`Match_destroy_all` during final cleanup. `FastaDB_close` had already freed the
alphabet. A completed-output footer therefore does not establish success.

The upstream `--disable-pthreads` configuration selects serial execution, but
the release contains two compilation defects in that configuration. The
companion patch guards an unused callback with `USE_PTHREADS` and calls the
`job_func` function parameter in the serial submission branch. It does not
change alignment scoring, thresholds, filtering, or models. Assertions remain
disabled, which is the upstream release default. Enabling assertions causes an
independent failure in the upstream protein self-comparison initialization and
does not diagnose the reported cleanup error.

`--cores 1` still uses the worker queue. Zero cores leave that queue without a
worker. Exhaustive alignment changes the algorithm and still uses a worker
queue when threading is compiled in. These are not equivalent runtime fixes.

## Build

Requirements are Python 3.11 or later, a working C compiler, `make`, `patch`,
`pkg-config` (or `pkgconf`), and GLib development headers and libraries. On
macOS the compiler normally comes from Xcode Command Line Tools. Use paths
without whitespace because the upstream build system embeds paths in compiler
commands. Both output directories must be new and separate.

From the repository root, after activating an environment that has these build
requirements:

```sh
python packaging/exonerate-serial/build.py \
  --glib-prefix "$CONDA_PREFIX" \
  --work-dir "$PWD/../exonerate-build" \
  --install-prefix "$PWD/../exonerate-serial"
```

If necessary, specify the real build tools with `--pkg-config /path/to/pkgconf`
and `--cc /path/to/clang`. The script downloads the original EBI source, verifies
its SHA256, applies the local patch, builds sequentially, and installs under the
specified prefix. Sequential make avoids an upstream generated-header race.
It retains the source archive, build tree, logs, and `build-provenance.json`.
The provenance records commands, return codes, source and patch hashes, compiler
flags, and the resulting binary hash. A build failure stops the script.

For a FuNLR run, set the configuration's `tools.exonerate` to the resulting
absolute `exonerate-serial/bin/exonerate` path. For installed demonstration and
diagnostic commands, place that directory first on `PATH`:

```sh
export PATH="$PWD/../exonerate-serial/bin:$PATH"
funlr doctor --demo
funlr demo --outdir serial-synthetic-demo
funlr demo --dataset public-sequences --outdir serial-public-demo
```

## Validation and limits

On macOS arm64, each of two separately built serial release binaries passed
three repeats of both a no-hit query and a positive query that exposed failures
with the stock binary (six successful commands per binary). A separate serial
AddressSanitizer build also passed all six repeats without a reported sanitizer
error. Repeated outputs were identical per query, and the positive scientific
output matched a successful stock run, excluding command and hostname metadata.
See `PROVENANCE.json` for the hashes and build settings. These are bounded runtime
checks; the exact mechanism of the original cleanup error remains unproven, and
the checks do not establish universal race freedom or Exonerate correctness.

This workaround requires a compiler and additional build tools and is not
automatically installed by `install.py`, the Conda recipe, or the containers.
Successful runs using this executable do **not** validate the unmodified locked
Exonerate package. Validate the lock-only installation, both demos, and a real
reference run on each advertised platform before publishing that installation
route. Linux builds and container behavior are separate checks.

The downloaded Exonerate source, this patch to that source, and the resulting
executable are licensed under GNU GPL version 3. The unchanged upstream license is
included as `COPYING`. The original `build.py` helper is MIT-licensed under
FuNLR's root `LICENSE`. Any redistribution of
the patched executable must meet that license's source-distribution requirements;
this directory is a build recipe and patch, not a published binary package.

Source: [EBI Exonerate distribution](https://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/).
The original archive checksum is also recorded in
[Bioconda's Exonerate recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/exonerate/meta.yaml).
