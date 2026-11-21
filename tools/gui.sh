#!/usr/bin/env bash
# tools/gui.sh [--rebuild] [design.json]
# Opens the BladeGen section designer window (Linux / macOS; needs a display).
#  * uses the conda / micromamba environment BLADEGEN_ENV (default: $CONDA_PREFIX, else ~/mamba/envs/bladegen,
#    else ~/miniconda3/envs/bladegen; create it once with  conda env create -f environment.yml)
#  * builds the Python module (preset 'python') when it is missing or older than this GUI, or always with --rebuild
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV=${BLADEGEN_ENV:-${CONDA_PREFIX:-}}
for c in "$ENV" "$HOME/mamba/envs/bladegen" "$HOME/miniconda3/envs/bladegen" "$HOME/micromamba/envs/bladegen"; do
    if [ -n "$c" ] && [ -x "$c/bin/python" ] && "$c/bin/python" -c "import tkinter, matplotlib" 2>/dev/null; then ENV=$c; break; fi
done
[ -x "$ENV/bin/python" ] || { echo "conda environment not found; create it with: conda env create -f environment.yml"; exit 1; }
REBUILD=0; DESIGN=()
while [ $# -gt 0 ]; do
    case "$1" in --rebuild) REBUILD=1 ;; *) DESIGN=("$1") ;; esac; shift
done
BUILD=${BLADEGEN_BUILD_ROOT:-$SRC/build}/python
export PATH="$ENV/bin:$PATH" PYTHONPATH="$BUILD/python${PYTHONPATH:+:$PYTHONPATH}"
if [ "$REBUILD" = 1 ] || ! "$ENV/bin/python" -c "import bladegen, bladegen.gui; bladegen.Naca65Section" >/dev/null 2>&1; then
    echo "building the BladeGen Python module (first time: a few minutes) ..."
    BLADEGEN_ENV=$ENV bash "$SRC/tools/linux_dev.sh" python --no-test
fi
exec "$ENV/bin/python" -m bladegen.gui "${DESIGN[@]}"
