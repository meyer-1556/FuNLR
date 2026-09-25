# Linux amd64 only. This is the verified platform manifest for micromamba 2.9.0.
# Build with `docker build --platform linux/amd64 ...`, including on Apple Silicon.
FROM mambaorg/micromamba:2.9.0@sha256:5681ae3caa12844c41d31e1d70f636fbaba31b01bdfd1b5b521742d7e490614d

LABEL org.opencontainers.image.source="https://github.com/meyer-1556/FuNLR" \
      org.opencontainers.image.title="FuNLR" \
      org.opencontainers.image.version="0.5.0a1" \
      org.opencontainers.image.licenses="MIT AND BSD-3-Clause AND CC-BY-4.0 AND CC0-1.0"

# Exact package URLs and hashes avoid a new dependency solve at build time.
COPY --chown=$MAMBA_USER:$MAMBA_USER environments/conda-linux-64.lock /tmp/funlr.lock
RUN micromamba install --yes --name base --file /tmp/funlr.lock && \
    micromamba clean --all --yes

WORKDIR /home/mambauser/funlr
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md LICENSE NOTICE ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
RUN micromamba run --name base python -m pip install --no-deps --no-build-isolation . && \
    micromamba run --name base funlr doctor --demo

# Explicit PATH also supports `apptainer exec` after conversion of this image.
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "funlr"]
CMD ["--help"]
