"""Execute every notebook in python/notebooks headless (ctest: python_notebooks).

Each notebook runs top to bottom with BLADEGEN_NB_STATIC=1 (widgets are skipped and
every plot is drawn once with its default values).  The run fails on any cell error
or on a notebook that produces no figure.  The .ipynb files are not modified.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
NB_DIR = HERE / "notebooks"


KERNEL = "bladegen-run"


def register_kernel() -> None:
    """A kernelspec for THIS interpreter (a global 'python3' spec may be another Python)."""
    root = Path(tempfile.mkdtemp(prefix="bladegen-kernels-"))
    d = root / "kernels" / KERNEL
    d.mkdir(parents=True)
    (d / "kernel.json").write_text(json.dumps({
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "BladeGen (run_notebooks)", "language": "python"}))
    os.environ["JUPYTER_PATH"] = str(root) + os.pathsep + os.environ.get("JUPYTER_PATH", "")


def run(path: Path) -> tuple[bool, str]:
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(nb, timeout=1200, kernel_name=KERNEL,
                            resources={"metadata": {"path": str(path.parent)}})
    t0 = time.perf_counter()
    try:
        client.execute()
    except Exception as exc:  # noqa: BLE001 - report every failure
        msg = str(exc)
        return False, f"{type(exc).__name__}: {msg[:600]} ... {msg[-2500:]}"
    figures = sum(1 for c in nb.cells if c.cell_type == "code"
                  for o in c.get("outputs", []) if "image/png" in o.get("data", {})
                  or "application/vnd.plotly.v1+json" in o.get("data", {}))
    if figures == 0:
        return False, "no figure produced"
    return True, f"{figures} figures, {time.perf_counter() - t0:.1f} s"


def main() -> int:
    os.environ.setdefault("BLADEGEN_NB_STATIC", "1")
    os.environ.setdefault("MPLBACKEND", "Agg")
    register_kernel()
    books = sorted(NB_DIR.glob("*.ipynb"))
    if len(sys.argv) > 1:
        books = [b for b in books if any(k in b.name for k in sys.argv[1:])]
    if not books:
        print("no notebooks found")
        return 0
    failed = 0
    for b in books:
        ok, msg = run(b)
        failed += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {b.name}: {msg}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
