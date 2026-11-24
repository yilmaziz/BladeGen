"""
turbolenza_lint.py
==================
Compiles the TURBOlenza sources (src/turbolenza, tests/turbolenza, apps/designer/turbolenza and the
designer files they touch) with clang++ -fsyntax-only and the warning flags of BladeGen's GCC/Clang CI
(cmake/CompilerFlags.cmake), so the Linux and macOS -Werror builds can be checked on Windows.

    tools\\turbolenza_lint.cmd [file ...]

Needs a clang++ (conda env "tl-lint": conda create -n tl-lint -c conda-forge clangxx) and the
bladegen env for the headers; run from an MSVC developer environment (the .cmd sets it up), since
clang on Windows uses MSVC's standard library.
"""

from __future__ import annotations

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV = Path(os.environ.get("BLADEGEN_ENV", Path.home() / "miniconda3" / "envs" / "bladegen"))
CLANG = Path(os.environ.get("TL_CLANG", Path.home() / "miniconda3" / "envs" / "tl-lint" / "Library" / "bin" / "clang++.exe"))

# cmake/CompilerFlags.cmake (GCC/Clang), plus macOS clang's -Wunused-const-variable
WARNINGS = ["-Wall", "-Wextra", "-Wpedantic", "-Wshadow", "-Wconversion", "-Wsign-conversion", "-Wnon-virtual-dtor",
            "-Wold-style-cast", "-Wcast-align", "-Wunused", "-Woverloaded-virtual", "-Wdouble-promotion", "-Wformat=2",
            "-Wunused-const-variable", "-Wshadow-all", "-Wno-unknown-warning-option"]  # -Wshadow-all ~ GCC -Wshadow
INC = ENV / "Library" / "include"
SYSTEM = [INC / "eigen3", INC, INC / "opencascade", INC / "qt6", INC / "qt6" / "QtCore", INC / "qt6" / "QtGui",
          INC / "qt6" / "QtWidgets", INC / "qt6" / "QtConcurrent"]
DEFINES = ["TURBOLENZA_HAVE_OPENMP", "TURBOLENZA_HAVE_PARDISO", 'BLADEGEN_SOURCE_DIR="."', 'BLADEGEN_TEST_OUTPUT="."',
           'BLADEGEN_VERSION="lint"', "NOMINMAX", "WIN32_LEAN_AND_MEAN", "_USE_MATH_DEFINES", "BLADEGEN_HAVE_OCCT_VIEWER",
           "_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH"]  # clang 19 with a newer MSVC STL


def files(args: list[str]) -> list[Path]:
    if args:
        return [Path(a).resolve() for a in args]
    out = sorted((ROOT / "src" / "turbolenza").glob("*.cpp")) + sorted((ROOT / "tests" / "turbolenza").glob("*.cpp"))
    out += sorted((ROOT / "apps" / "designer" / "turbolenza").glob("*.cpp"))
    out += [ROOT / "apps" / "designer" / name for name in ("PlotWidget.cpp", "DesignerWindow.cpp", "main.cpp")]
    return out


def check(path: Path) -> tuple[Path, int, str]:
    cmd = [str(CLANG), "-fsyntax-only", "-std=c++2b", "-fopenmp", "-fms-extensions", *WARNINGS,
           f"-I{ROOT / 'src'}", f"-I{ROOT / 'tests' / 'turbolenza'}", f"-I{ROOT / 'apps' / 'designer'}",
           *[f"-isystem{p}" for p in SYSTEM], *[f"-D{d}" for d in DEFINES], str(path)]
    run = subprocess.run(cmd, capture_output=True, text=True)
    return path, run.returncode, run.stderr


def main() -> int:
    todo = files(sys.argv[1:])
    warnings = errors = 0
    with ThreadPoolExecutor(max_workers=8) as pool:
        for path, code, err in pool.map(check, todo):
            own = [line for line in err.splitlines() if ("warning:" in line or "error:" in line)
                   and str(ROOT).lower().replace("\\", "/") in line.lower().replace("\\", "/")]
            if own or code:
                print(f"== {path.relative_to(ROOT)} (exit {code})")
                for line in own or err.splitlines()[:20]:
                    print("  " + line)
            warnings += sum("warning:" in line for line in own)
            errors += sum("error:" in line for line in own) + (1 if code and not own else 0)
    print(f"{len(todo)} files: {warnings} warnings, {errors} errors in BladeGen files")
    return 1 if warnings or errors else 0


if __name__ == "__main__":
    sys.exit(main())
