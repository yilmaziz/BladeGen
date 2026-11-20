"""BladeGen: parametric C2 turbomachinery blades (C++ engine, OpenCASCADE).

Everything geometric is computed by the C++ library (module ``bladegen._core``);
this package only re-exports it and, on Windows, makes the OpenCASCADE DLLs
findable before the extension is imported.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _add_dll_dirs() -> None:
    if os.name != "nt":
        return
    dirs = []
    try:
        from . import _paths  # written by CMake
        dirs.append(_paths.OCCT_BIN_DIR)
    except Exception:  # pragma: no cover - only in a broken install
        pass
    dirs.append(os.path.join(sys.prefix, "Library", "bin"))
    dirs += [d for d in os.environ.get("BLADEGEN_DLL_DIRS", "").split(os.pathsep) if d]
    for d in dirs:
        if d and Path(d).is_dir():
            try:
                os.add_dll_directory(d)
            except OSError:
                pass


_add_dll_dirs()

from ._core import *  # noqa: E402,F401,F403
from ._core import __version__  # noqa: E402,F401

try:
    from . import _paths as _p

    SOURCE_DIR = Path(_p.SOURCE_DIR)
except Exception:  # pragma: no cover
    SOURCE_DIR = Path(__file__).resolve().parent
