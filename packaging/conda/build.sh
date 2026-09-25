#!/bin/sh
set -eu

# meta.yaml invokes this through bash, so browser uploads do not need to
# preserve its executable bit. Conda supplies the build interpreter.
"${PYTHON}" -m pip install --no-deps --no-build-isolation --ignore-installed .
