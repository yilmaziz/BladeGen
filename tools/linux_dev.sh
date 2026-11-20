#!/usr/bin/env bash
# tools/linux_dev.sh <preset> [--no-test] — Linux / macOS developer build inside a conda/micromamba env
# (BLADEGEN_ENV, default ~/mamba/envs/bladegen) providing OCCT, Eigen, Catch2, cmake, ninja, g++.
set -euo pipefail
PRESET=${1:-linux-release}
RUN_TESTS=1
[ "${2:-}" = "--no-test" ] && RUN_TESTS=0
ENV=${BLADEGEN_ENV:-$HOME/mamba/envs/bladegen}
export PATH="$ENV/bin:$PATH"
if [ -z "${CXX:-}" ]; then
    CXX=$(ls "$ENV"/bin/*-g++ "$ENV"/bin/*-clang++ 2>/dev/null | head -1 || true)
    [ -n "$CXX" ] || CXX=$(command -v g++ || command -v clang++)
fi
export CXX
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD=${BLADEGEN_BUILD_ROOT:-$SRC/build}
cd "$SRC"
cmake --preset "$PRESET" -B "$BUILD/$PRESET" -DCMAKE_PREFIX_PATH="$ENV" -DPython_EXECUTABLE="$ENV/bin/python"
cmake --build "$BUILD/$PRESET" --parallel
if [ "$RUN_TESTS" = 1 ]; then ctest --test-dir "$BUILD/$PRESET" --output-on-failure; fi
