#!/usr/bin/env bash
# tools/designer.sh [--rebuild] [design.json]
# Opens the BladeGen Designer (C++ / Qt, MDI) on Linux / macOS (needs a display).
#  * uses the conda / micromamba environment BLADEGEN_ENV (default: $CONDA_PREFIX, else ~/mamba/envs/bladegen,
#    else ~/miniconda3/envs/bladegen; create it once with  conda env create -f environment.yml)
#  * builds bladegen_designer (preset linux-release / macos-release) when it is missing, or always with --rebuild
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV=${BLADEGEN_ENV:-${CONDA_PREFIX:-}}
for c in "$ENV" "$HOME/mamba/envs/bladegen" "$HOME/miniconda3/envs/bladegen" "$HOME/micromamba/envs/bladegen"; do
    if [ -n "$c" ] && [ -d "$c/lib/cmake/Qt6" ]; then ENV=$c; break; fi
done
[ -d "$ENV/lib/cmake/Qt6" ] || { echo "Qt 6 not found; create / update the environment: conda env update -f environment.yml"; exit 1; }
PRESET=linux-release
[ "$(uname)" = Darwin ] && PRESET=macos-release
REBUILD=0; DESIGN=()
while [ $# -gt 0 ]; do
    case "$1" in --rebuild) REBUILD=1 ;; *) DESIGN=("$1") ;; esac; shift
done
BUILD=${BLADEGEN_BUILD_ROOT:-$SRC/build}/$PRESET
EXE=$BUILD/bin/bladegen_designer
if [ "$REBUILD" = 1 ] || [ ! -x "$EXE" ]; then
    echo "building the BladeGen Designer (first time: a few minutes) ..."
    BLADEGEN_ENV=$ENV bash "$SRC/tools/linux_dev.sh" "$PRESET" --no-test
fi
export PATH="$ENV/bin:$PATH" LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$EXE" "${DESIGN[@]}"
