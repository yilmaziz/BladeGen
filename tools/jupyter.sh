#!/usr/bin/env bash
# tools/jupyter.sh [--rebuild] [--port N]
# Opens the BladeGen notebooks in JupyterLab (Linux / macOS).
#  * uses the conda / micromamba environment BLADEGEN_ENV (default: $CONDA_PREFIX, else ~/mamba/envs/bladegen,
#    else ~/miniconda3/envs/bladegen; create it once with  conda env create -f environment.yml)
#  * builds the Python module (preset 'python') when it is missing, or always with --rebuild
#  * starts JupyterLab in python/notebooks with the module on PYTHONPATH
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV=${BLADEGEN_ENV:-${CONDA_PREFIX:-}}
for c in "$ENV" "$HOME/mamba/envs/bladegen" "$HOME/miniconda3/envs/bladegen" "$HOME/micromamba/envs/bladegen"; do
    if [ -n "$c" ] && [ -x "$c/bin/python" ] && "$c/bin/python" -c "import jupyterlab" 2>/dev/null; then ENV=$c; break; fi
done
[ -x "$ENV/bin/python" ] || { echo "conda environment not found; create it with: conda env create -f environment.yml"; exit 1; }
PORT=8888; REBUILD=0
while [ $# -gt 0 ]; do
    case "$1" in --rebuild) REBUILD=1 ;; --port) PORT=$2; shift ;; esac; shift
done
BUILD=${BLADEGEN_BUILD_ROOT:-$SRC/build}/python
if [ "$REBUILD" = 1 ] || ! ls "$BUILD"/python/bladegen/_core*.so >/dev/null 2>&1; then
    echo "building the BladeGen Python module (first time: a few minutes) ..."
    BLADEGEN_ENV=$ENV bash "$SRC/tools/linux_dev.sh" python --no-test
fi
export PATH="$ENV/bin:$PATH" PYTHONPATH="$BUILD/python${PYTHONPATH:+:$PYTHONPATH}"
echo "opening JupyterLab on http://localhost:$PORT  (Ctrl+C stops it)"
exec "$ENV/bin/python" -m jupyter lab --notebook-dir="$SRC/python/notebooks" --port="$PORT"
