# TURBOlenza v2 in BladeGen: the port of the v1 solver

TURBOlenza v2 is the blade-to-blade compressible potential-flow solver of Aungier (2003), Ch. 5,
running inside BladeGen as its own window. It is a C++ port of the verified Python package
`turbolenza.b2b` (TURBOlenza v1, `C:\Users\ASUS\Documents\turbolenza`, commit 67599a6).
BladeGen stays the main program: TURBOlenza takes its blade section from a BladeGen design layer.

## What is ported, and what is not

| v1 module | v2 file (`src/turbolenza/`) | notes |
|---|---|---|
| geometry.py: `ParametricBladeCurve`, `SurfaceCurves` | `Curve.hpp/.cpp` | done (contract) |
| geometry.py: cascade generation (`make_cascade` minus the blade, `place_cascade`, throat, overlap check, `validate_cascade` frame checks) | `Cascade.hpp/.cpp` | the NACA 65 construction of v1 is **not** ported: the blade comes from BladeGen |
| grid.py (C2 path, sections 1-9) | `Grid.hpp/.cpp`, data in `GridData.hpp` | the C1 / PCHIP path is not ported (v2 is C2 only) |
| config.py, metrics.py, coefficients.py, discretization.py, boundary_conditions.py, matrix_assembly.py, linear_solver.py, flow.py, eos.py, kutta.py, convergence.py, solver.py, postprocess.py | one `.hpp/.cpp` each, same names in CamelCase | production paths only (below) |
| controller.py | `Case.hpp/.cpp` | the run layer the window calls |
| monitor.py, plots.py | the Qt window (`apps/designer/turbolenza/`) | |

Production paths only. The v1 study options are not ported: `transition_form = "flux"`,
`kutta_eval = "fixed_m"`, `kutta_method = "aungier_13_57*"`. v2 keeps the v1 defaults:
expanded transition rows, the Kutta condition at the last blade station, and the "current"
correction (Aungier 2000 Eqs. 13-57 and 13-59 in our frame).

## Rules

1. **Exact port.** v2 uses the same equations, discretization, boundary relations, node
   numbering, iteration order, relaxations, tolerances and refusals as v1. Each function cites
   two things in its header: the v1 function it ports and Aungier's equation number, copied from
   the v1 docstring. Do not "improve" a formula. A v1 behaviour that looks wrong goes into the
   report, unchanged in the code.
2. **Refuse, never repair** (v1 rule). An inadmissible grid, a singular matrix, a sonic node or a
   non-finite state throws a typed exception with v1's message. Nothing is clipped, floored or
   regularized.
3. **Tests that can fail.** Each module gets two kinds of test:
   - a check against v1 reference data (`tests/turbolenza/data/*.json`, written by
     `tools/turbolenza_reference/dump_v1.py`, git-ignored; such a test SKIPs when the data is
     missing);
   - at least one check that does not need v1: a closed form, the module's defining equation,
     or degenerate input that must be refused.

   Tolerances follow "measure before tolerance": measure the C++-vs-v1 difference, then set
   the bound at about 10x the measured value. Write the measured value in a comment next to
   the bound. Round-off differences against NumPy/SciPy are expected: SuperLU against Eigen
   SparseLU, and scipy's de Boor against our basis recurrence.
4. **Threads do not change the answer.** OpenMP parallelizes only per-node maps
   (`#pragma omp parallel for` over i, with no reductions in floating point). Sums, integrals,
   medians and mass-flow integrals stay serial, in v1's order. The same case gives bit-identical
   results on 1 and N threads with SparseLU. PARDISO is a different factorization: its results
   agree with SparseLU to round-off, not bit for bit, and the report says so. MSVC has OpenMP
   2.0 only: signed loop indices, no `collapse`.
5. **Zero warnings.** MSVC `/W4 /WX` (`BLADEGEN_WARNINGS_AS_ERRORS=ON`), and clean code for GCC
   `-Wall -Wextra -Wconversion -Wsign-conversion -Wshadow -Wold-style-cast` (CI builds Linux and
   macOS).
6. **Layout.** Namespace `BladeGen::Turbolenza`. Fields are `Field` (row-major
   `(n_m, n_eta)`, the v1 C order). Angles are in radians inside the code; a `_deg` name is
   degrees. The style follows BladeGen's `src/` (see `naca/Naca65.cpp`).

## Build and test (Windows, the `bladegen` conda env)

```
tools\turbolenza_build.cmd <build-dir> test_turbolenza ["catch2 filter"]
```

This configures `build\<build-dir>` once (MSVC Release, `/WX`), builds the target, and runs it
when the target is `test_turbolenza`. Use your own build dir so parallel work does not collide.
The designer target is `bladegen_designer`.

Generate the reference data again with:

```
C:\Users\ASUS\miniconda3\envs\newblade\python.exe tools\turbolenza_reference\dump_v1.py C:\Users\ASUS\Documents\turbolenza
```

## Aungier's grid rules (looked up 2026-11-24, `aungier-reference`)

- The grid lines at the far inlet and exit stations lie along the flow angle. Between the far
  station and the LE (TE), the side-line angle beta varies linearly in m (Aungier 2003 Sec. 5.3,
  p. 86, Eq. 5-33; 2000 Sec. 13.1, p. 237). v1 already does this, and v2 keeps it.
- The first and last grid points on the blade surfaces lie inside the blade passage, with one
  node on each surface (2003 pp. 89-90, Fig. 5-5; 2000 pp. 247-249, Figs. 13-3 and 13-4). No
  distance is prescribed. v2 enforces this on every grid the solver gets (`require_ss_ps`).
- Constant dm and deta (2003 p. 87). Neither book gives an extension length.
