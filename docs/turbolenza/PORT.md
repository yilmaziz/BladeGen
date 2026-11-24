# TURBOlenza v2 in BladeGen: the port of the v1 solver

TURBOlenza v2 is the blade-to-blade compressible potential-flow solver of Aungier (2003), Ch. 5,
running inside BladeGen as its own window. It is a C++ port of the verified Python package
`turbolenza.b2b` (TURBOlenza v1, github.com/yilmaziz/TURBOlenza, commit 67599a6).
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

## Not a port: the boundary layer

v1 has no boundary layer. `BoundaryLayer.hpp/.cpp` (Aungier 2003 Sec. 5.6 = Aungier 2000 Ch. 15)
and `Viscous.hpp/.cpp` (its two-way coupling by transpiration) are new. Rules 2 to 6 apply
unchanged. In place of rule 1, every equation cites Aungier's equation number as read from the
page images, and every choice Aungier does not make is marked "ours" in the header. With
`CaseSettings::viscous` off, nothing of it runs: the solve, its records and its files are v1's,
bit for bit. The only changes to the ported modules are the transpiration values in the boundary
relations and the periodic continuation (`Transpiration`, `BoundaryConditions::station_jump`) and
psi_m on the walls in the flow recovery. Without a `Transpiration` both are v1's.

## Not a port: the stream sheet and rotating rows

v1 solves the planar, stationary cascade (r = b = 1, omega = 0). v2 adds Aungier's full Eqs. 5-3,
5-9, 5-15 to 5-26 and 5-51 (2003, printed pp. 82-89, read from the page images) on an axisymmetric
stream sheet r(m), b(m), and rotating rows, plus Aungier 2000 Eqs. 13-57 and 13-59 with their r
(p. 245). Rules 2 to 6 apply; every equation cites its number, every choice of ours is marked.

- `StreamSheet.hpp` holds r(m), b(m) as tables (cubic spline). The grid stays planar in
  (m, y = r_LE theta); the sheet enters through sigma = r / r_LE: S = sigma (y1 - y0),
  tan(beta) = sigma dy/dm, rho b in A, B, C, the source Q = 2 S omega sin(phi), mu = mdot' / (rho b),
  H' = I + (r_i omega)^2 / 2 per station, the Kutta step times r_te / r_out.
- The side lines integrate r d theta / dm = tan(beta) numerically (Gauss-Legendre).
- Without a sheet (and omega = 0) nothing of this runs: the solve is v1's, bit for bit. A cylinder
  whose radius is r_LE is the planar cascade bit for bit (tested).
- The sign convention: TURBOlenza's +theta is Aungier's, omega = -2 pi rpm / 60. Positive rpm moves
  the blades toward -theta, so with axial inflow beta1 > 0 (a compressor rotor). `Flow.hpp`
  velocity_triangle states it.
- `BladeGenSource.hpp` stream_surface_blade places a BladeGen layer on its stream surface. b(m) is
  the normal spacing of BladeGen's streamline family (our model).
- The boundary layer is refused on a sheet that is not a cylinder (Eq. 5-106 carries b(x); the
  march has b = 1) and on rotating rows.

## Build and test (Windows, the `bladegen` conda env)

```
tools\turbolenza_build.cmd <build-dir> test_turbolenza ["catch2 filter"]
```

This configures `build\<build-dir>` once (MSVC Release, `/WX`), builds the target, and runs it
when the target is `test_turbolenza`. Use your own build dir so parallel work does not collide.
The designer target is `bladegen_designer`.

Generate the reference data again with:

```
python tools\turbolenza_reference\dump_v1.py <path-to-turbolenza-v1>
```

## Aungier's grid rules (looked up 2026-11-24, `aungier-reference`)

- The grid lines at the far inlet and exit stations lie along the flow angle. Between the far
  station and the LE (TE), the side-line angle beta varies linearly in m (Aungier 2003 Sec. 5.3,
  p. 86, Eq. 5-33; 2000 Sec. 13.1, p. 237). v1 already does this, and v2 keeps it.
- The first and last grid points on the blade surfaces lie inside the blade passage, with one
  node on each surface (2003 pp. 89-90, Fig. 5-5; 2000 pp. 247-249, Figs. 13-3 and 13-4). No
  distance is prescribed. v2 enforces this on every grid the solver gets (`require_ss_ps`).
- Constant dm and deta (2003 p. 87). Neither book gives an extension length.

## Not a port: the general curvilinear operator

v1 and Aungier solve the expanded equation (Eq. 5-21) on Aungier's sheared H-grid only.
`GeneralMesh`, `GeneralOperator`, `GeneralSolver`, `GeneralAdapter` and `MeshFormat` are new: the same
physics, in conservative (flux) form, on any structured multi-block mesh of one stream sheet in the
conformal plane (m' = ∫dm/r, θ), with blade walls, periodic and interface connections, inlet and outlet
flow-angle conditions and the Kutta condition by trailing-edge tags. The derivation and the scheme are in
`GeneralOperator.hpp`, the solve and the rotation sign in `GeneralSolver.hpp`, the file format in
`docs/turbolenza/mesh_format.md`. Rules 2 to 6 apply unchanged; the density loop reuses Eos, Convergence,
Kutta (`next_exit_angle`) and LinearSolver. `solve_general` is a separate entry point: nothing of the
Aungier solve changes. Choices that are ours are marked in the headers (the run rule of the nodal differences
along a boundary, the cell-vertex control volumes and their closure face at H-grid transition stations,
the general form of the exit-angle correction).

At an H-grid transition station (the last inlet and the first exit station) the two periodic images of a
node meet the suction and the pressure wall, a blade thickness h_f apart; the control volume is closed by
one face across that gap (GeneralOperator.hpp).  Without it (0ede48c) the discrete metric identities failed
there by h_f / 2 and the scheme was first order (sharp edge) or about half order (round edge); with it they
hold to round-off and the manufactured H-grid test is second order with a sharp edge and first order with a
round one.  Aungier's own transition rows are of the same low order (independent review of 0ede48c): on a
manufactured problem both converge to the exact solution with different error constants.

Open items of the general operator (independent review of 407b811, not settled):

1. On the NACA 65 cascade at fixed angle the inlet-corner psi(0,0) of the two solvers differs by a steady
   0.0195 (M1 0.3: Aungier -0.0421, -0.0444, -0.0450; general -0.0631, -0.0644, -0.0644 at very_coarse,
   coarse, medium), each changing by less than 1e-3 over the last level.  The manufactured tests point at
   Aungier's transition rows (order 0.2 to 0.5 with a round edge, nearly stalling for one edge position,
   against order 1 for the general operator); not proven for the cascade.  Needs fine / very_fine levels
   or a leading-edge topology without transition cells (an O- or C-grid).
2. beta_out depends by about 0.9 deg on how the trailing-edge node speeds are differenced (the run rule;
   M1 0.3 medium: 24.91 with it, 24.02 without): the Kutta condition reads speeds that do not converge next
   to a round edge on this layout.  Both solvers share this limitation of Aungier's grid.
3. The synthetic Kutta cascade test bounds beta_out only to about 1.5 deg.
