# Containers, Conda packaging, and installation checks

The Docker and Apptainer recipes use the checked Linux x86_64 environment lock and the same application as the Python installer. They currently target Linux amd64; Apple Silicon users can build/run that Docker image through Docker's amd64 support, while the native installer uses the separate macOS arm64 lock. No FuNLR image has been published to a registry.

Both recipes pin `mambaorg/micromamba:2.9.0` to the Linux amd64 OCI manifest digest `sha256:5681ae3caa12844c41d31e1d70f636fbaba31b01bdfd1b5b521742d7e490614d`. The tag's index and this platform manifest were fetched from the official Docker Hub registry on 2026-09-09, and their SHA-256 values were checked against the manifest bytes. Package builds come from `environments/conda-linux-64.lock`, which lists explicit URLs and checksums. This avoids resolving a new set of dependencies during a container build. [Micromamba lockfile documentation](https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html#using-a-lockfile)

## Docker

From the repository root, with Docker available:

```bash
docker build --platform linux/amd64 --tag funlr:0.5.0a1 .
docker run --rm --platform linux/amd64 funlr:0.5.0a1 doctor --demo
mkdir -p "$PWD/docker-output"
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD/docker-output,dst=/output" \
  funlr:0.5.0a1 demo --outdir /output/demo
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD/docker-output,dst=/output" \
  funlr:0.5.0a1 demo --dataset public-sequences --outdir /output/public-demo
```

The explicit platform selects the Linux amd64 dependencies, including when the host is an ARM Mac. Running as the invoking user's UID/GID keeps bind-mounted result files writable by that user on Linux. The test writes `docker-output/demo/evidence.json`; successful evidence records eight completed stages, eight reused stages, and positive Exonerate output. [Docker build platform option](https://docs.docker.com/reference/cli/docker/buildx/build/#platform), [Docker bind mounts](https://docs.docker.com/engine/storage/bind-mounts/)

For real data, bind the input/database directories read-only and an output directory writable, then use the normal `funlr run` arguments after the image name. All paths passed to FuNLR must be the paths visible inside the container. For example, `--mount type=bind,src=/absolute/host/databases,dst=/databases,readonly` makes the HMM files available beneath `/databases`. The image includes the tools and small demo-only profiles; the user supplies the full scientific HMM databases for production runs.

## Apptainer

On a Linux x86_64 machine with Apptainer build support, start at the repository root:

```bash
apptainer build --fakeroot funlr-0.5.0a1.sif Apptainer.def
apptainer exec --cleanenv funlr-0.5.0a1.sif funlr doctor --demo
mkdir -p "$PWD/apptainer-output"
apptainer run --cleanenv \
  --bind "$PWD/apptainer-output:/output" \
  funlr-0.5.0a1.sif demo --outdir /output/demo
apptainer run --cleanenv \
  --bind "$PWD/apptainer-output:/output" \
  funlr-0.5.0a1.sif demo --dataset public-sequences --outdir /output/public-demo
```

The definition builds directly from the pinned base and explicit lock; it does not require a prepublished FuNLR image or a Docker daemon. Its `%test` checks the Python installation and external executables. The full demo is then run against the resulting SIF. If an HPC site does not permit local image builds, build the SIF on a permitted Linux builder and copy it to the cluster. Apptainer's build permissions and namespace support are host requirements. [Apptainer definition files](https://apptainer.org/docs/user/latest/definition_files.html), [Apptainer installation requirements](https://apptainer.org/docs/admin/latest/installation.html)

For a scientific run, add binds such as `--bind /host/data:/data:ro,/host/databases:/databases:ro,/host/results:/output` and pass `run` plus the input flags after the SIF filename. The same private working-copy and resume rules apply inside either container.

## Conda package

`packaging/conda/meta.yaml` is a local-source recipe for this checkout, with the MIT project license, retained BSD-3-Clause code notice, CC-BY-4.0/CC0-1.0 demo-data notices, Python entry point, external-tool dependencies, and both real-tool demonstrations as its package tests. The [recipe instructions](../packaging/conda/README.md) describe building and installing a local channel without publishing it.

The recipe uses `source.path` because a public release archive does not yet exist. After release, replace that path with the actual source archive URL and its verified checksum before submitting to a public channel. A Conda package's direct requirements do not fix every transitive package build; the explicit environment locks and installer remain the route for recreating the selected environment exactly. [Conda recipe metadata](https://docs.conda.io/projects/conda-build/en/latest/resources/define-metadata.html)

## CI evidence and its limits

`.github/workflows/tests.yml` runs ordinary regression tests on Python 3.11 and 3.12 and builds a wheel and source distribution. Its installation matrix exercises `install.py` on Linux x86_64 (`ubuntu-24.04`) and native Apple Silicon (`macos-14`), using each platform's lock. Both jobs are configured to reinstall the built wheel outside the checkout and run `doctor`, the synthetic installation demo, and the packaged public-sequence example. A separate job builds the Docker image and runs both cases with a writable output bind. The synthetic case requires positive Exonerate output and verified resume. The public case requires eight completed stages, six candidates, three strict candidates, and checked refinement alignments. The workflow saves distribution files and test evidence as GitHub Actions artifacts. It does not push images, upload to package channels, or require registry credentials. [GitHub runner architectures](https://docs.github.com/en/actions/reference/runners/github-hosted-runners)

`.github/workflows/apptainer.yml` automatically builds, inspects, and runs both demonstrations in the SIF on pushes or pull requests that change its definition, environment locks, application source, package metadata, license notices, or the workflow itself. Uploading the code batch to the new branch therefore starts this check without modifying `main`. It installs Apptainer 1.5.3 from the official release, verifies the package's SHA-256, and saves installation evidence plus the SIF hash. The SIF itself remains temporary. All CI jobs have a 30-minute timeout. Manual dispatch is also available once the workflow exists on the repository's default branch. [GitHub workflow triggers](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore), [GitHub manual workflows](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)

Recipe and workflow configuration is not evidence that an image or published package works. The current [validation record](VALIDATION.md) distinguishes completed runtime checks from pending fresh-prefix, hosted CI, container and Conda-build checks. A passing fixture establishes execution and behavior on that fixture; it does not establish biological accuracy or production-scale suitability. Record successful CI URLs and image/package hashes before describing a distribution route as validated.
