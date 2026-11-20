#!/usr/bin/env bash
# tools/linux_dev.sh <preset> — Linux developer build inside a conda/micromamba env
# (BLADEGEN_ENV, default ~/mamba/envs/bladegen) providing OCCT, Eigen, Catch2, cmake, ninja, g++.
set -euo pipefail
PRESET=${1:-linux-release}
ENV=${BLADEGEN_ENV:-$HOME/mamba/envs/bladegen}
export PATH="$ENV/bin:$PATH"
export CXX=${CXX:-$(ls "$ENV"/bin/*-g++ 2>/dev/null | head -1 || command -v g++)}
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD=${BLADEGEN_BUILD_ROOT:-$SRC/build}
cd "$SRC"
cmake --preset "$PRESET" -B "$BUILD/$PRESET" -DCMAKE_PREFIX_PATH="$ENV" -DPython_EXECUTABLE="$ENV/bin/python"
cmake --build "$BUILD/$PRESET" --parallel
ctest --test-dir "$BUILD/$PRESET" --output-on-failure
