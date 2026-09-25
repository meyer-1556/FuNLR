# Contributing to FuNLR

For a bug, include `funlr version`, `funlr doctor --demo`, installation method, configuration, the failing stage log and expected versus observed behavior. Share only data you can disclose.

Install development dependencies with `python -m pip install -e '.[dev]'`, then run `python -m pytest tests -q`. With real tools installed, run both `funlr demo --outdir /path/to/new-synthetic-demo` and `funlr demo --dataset public-sequences --outdir /path/to/new-public-demo`. See [validation](docs/VALIDATION.md) for the limits of these checks.

Keep execution fixes separate from scientific rule changes. Changes to thresholds, domain vocabularies, tiering, rescue interpretation or model libraries require documented rationale, before/after regression comparisons and independent biological evaluation. Do not replace reference outputs solely to make changed behavior pass. Record missing and uncertain evidence explicitly. Incomplete positive reference panels must not be used to estimate specificity or precision.

Add tests for meaningful failure modes. Keep external-tool doubles in tests; production commands must execute and check real tools. Avoid new runtime dependencies unless their benefit warrants the installation burden. Document API/output schema changes and check the installed wheel outside the source directory.

Original code is MIT licensed. Preserve third-party notices, scientific attribution and dataset terms. Adding code does not authorize redistribution of external model libraries. For substantial changes, describe the proposal in an issue and refer to the [roadmap](docs/ROADMAP.md). Public documentation should explain functions and evidence without internal development shorthand.
