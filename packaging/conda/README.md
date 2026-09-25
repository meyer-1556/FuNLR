# Local Conda package build

This recipe builds the current checkout and includes all Python and external-tool runtime dependencies. It does not refer to an invented release archive or upload anything to a channel. Run from the repository root on Linux x86_64, with `conda-build` installed in a separate builder environment:

```bash
conda create --name funlr-builder --override-channels --channel conda-forge conda-build
conda run --name funlr-builder conda build packaging/conda \
  --override-channels --channel conda-forge --channel bioconda \
  --output-folder /absolute/path/outside-the-checkout/funlr-conda \
  --no-anaconda-upload
```

Replace the output path with a writable directory outside this checkout. The recipe's relative source path resolves to the checkout root. Its build script is invoked explicitly through a shell, so a source ZIP or browser upload works without executable mode bits. Conda's package test imports FuNLR, checks the tools, and runs the synthetic and public-sequence real-tool demos.

To install the resulting local channel, use its absolute path as the first channel:

```bash
conda create --name funlr-local --override-channels \
  --channel file:///absolute/path/outside-the-checkout/funlr-conda \
  --channel conda-forge --channel bioconda funlr=0.5.0a1
conda run --name funlr-local funlr doctor --demo
```

The recipe pins direct application dependencies. Conda resolves their transitive dependencies at package installation; use the checked explicit environment locks and `install.py` when the exact tested environment is required. Local package building and installation must be verified before channel submission. See [container and CI checks](../../docs/CONTAINERS.md).

This local recipe currently targets Linux. The separately solved Apple Silicon environment lock uses samtools 1.24 instead of 1.21 to resolve its R dependency set; this recipe has not been adapted or tested for that platform. The optional serial Exonerate workaround is also separate from this recipe.
