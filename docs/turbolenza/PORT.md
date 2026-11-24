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

### The Fig. 5-15 study (2026-11-24, `feat/turbolenza-bl-study`)

**Question.** The medium-grid one-way run matches Fig. 5-15 on the pressure surface: transition at
x/c 0.180 against 0.173. On the suction surface it does not: transition at 0.213 against 0.130,
separation at 0.939 against 0.746. Why?

**Answer.** We found no defect in the boundary-layer code. Two things explain the gap:

- **The abscissa (our reading of the figure).** Fig. 5-15 plots against Aungier's blade-station fraction
  f = (m − m_first)/(m_last − m_first) on 24 stations, not against the chordwise x/c of each wall node.
  The evidence is in `validation/aungier2003_fig5_15/README.md`, "Which abscissa". Our transition at x/c
  0.213 is at f = 0.146, one station after f = 0.122, where Re_θ = 248 just misses 250.
- **The inviscid input.** Fed with Aungier's own edge Mach numbers (Fig. 5-6, digitized here, on his
  24 stations), `march_surface` lands on the figure's stations:

  | | Fig. 5-15 | ours on Fig. 5-6's Mach |
  |---|---|---|
  | SS transition | station 3 (0.130) | station 3, on the threshold (see below) |
  | SS separation | station 17 (0.739) | station 17 |
  | PS transition | station 4 (0.174) | station 4 |

  The SS transition cell is not firm. Re_θ at station 3 is 250.0 to 250.35, depending on the digitization.
  An independent digitization of Fig. 5-6 gives station 4, and with ±0.004 noise station 3 comes out in
  only 55 % of draws. The SS separation at station 17 (70 % under noise) and the PS transition at station 4
  (100 %) are robust.

  Our suction-surface Mach distribution differs from Fig. 5-6 on all five grid levels (the levels agree with
  each other within 0.01 for 0.2 ≤ f ≤ 0.8). Near the LE ours accelerates more slowly (f = 0.043:
  0.837 against 0.870). In 0.3 ≤ f ≤ 0.6 ours is 0.03 to 0.05 lower; in 0.78 ≤ f ≤ 0.95 it is 0.01 to
  0.07 higher (medium, f = 0.5: 0.755 against 0.803; f = 0.87: 0.590 against 0.537). So our suction
  surface diffuses less before f ≈ 0.8 and then abruptly at the end. Our separation (f 0.90 to 0.95)
  follows that abrupt rear diffusion. On our input it stays within 0.902 to 0.951 for every boundary-layer
  ingredient in the table below.
- **The H levels still differ.** The stations match, but on Aungier's input our H does not reproduce the
  figure's levels (the independent check found this):
  - Our turbulent H before separation is systematically low. On SS stations 12 to 16 it is 0.05 to 0.18
    below the figure (station 16: 2.42 to 2.46 against 2.60). On the PS, for f 0.5 to 0.9, it is 0.03 to
    0.05 low.
  - Our laminar H at station 1 misses too: SS 2.99 against 2.88, PS 3.00 against 3.21 (see "Discretization
    of the laminar march" below).

  So the agreement is in the stations, not in the H levels.

**How firm this is.** The Fig. 5-6 Mach numbers are our own rough digitization: ±0.004, with the f = 0
values extrapolated. One of its points, f 0.7958, lies on the experiment's dashed curve, not on the
potential-flow line (0.611 against 0.599). Re_θ at suction-surface station 3 is 250.3 on our
digitization, right at the threshold. With ±0.004 of random noise on the Mach numbers, the transition stays
at station 3 in 64 % of 300 draws (station 4 in the rest); the independent check found 55 %. The
separation stays at station 17 in 71 % (18: 11 %, 19: 18 %). Aungier's grid (24 stations,
uniform in m) is inferred from the vertices. Our arc length per station is applied to his grid.

**Ingredients, one at a time.** Each cell gives the transition station f, the Re_θ = 250 crossing f in
brackets, and the separation f. "Ours" is the medium grid with 42 stations. "Fig. 5-6" is Aungier's Mach
on 24 stations. SS = suction surface, PS = pressure surface.

| variant | SS, ours | SS, Fig. 5-6 | PS, ours | PS, Fig. 5-6 |
|---|---|---|---|---|
| as implemented (Re_c 1.034e6) | 0.146 (0.124) / 0.927 | 0.130 (0.130) / 0.739 | 0.171 | 0.174 |
| Re_c 5.5e5 (at M1 0.65) | 0.220 / 0.927 | 0.261 / 0.739 | 0.293 | 0.261 |
| Re_c × 1.28 (P0 = 1 atm) | 0.122 / 0.927 | 0.130 / 0.826 | 0.146 | 0.130 |
| c_f 1 + Λ/6 (Eq. 5-123 as printed) | 0.220 / 0.951 | 0.261 / 0.826 | 0.317 | 0.304 |
| F cubic coefficient 1.61875 (derived) | 0.146 / 0.927 | 0.130 / 0.739 | 0.171 | 0.174 |
| K without b0 | 0.146 / 0.927 | 0.174 / 0.826 | 0.171 | 0.174 |
| Re_θ on the inlet ρ1 W1 / μ1 | 0.171 / 0.927 | 0.217 / 0.826 | 0.146 | 0.130 |
| θ at station 0 = flat plate over the LE arc | 0.098 / 0.927 | 0.087 / 0.739 | 0.171 | 0.130 |
| laminar march in (ρ W² θ)² | 0.122 / 0.927 | 0.130 / 0.739 | 0.171 | fails |
| H_k separation 2.2 / 2.6 | 0.902 / 0.927 (sep) | 0.739 / 0.826 (sep) | – | – |
| E × 0.9 / × 1.1 | 0.927 / 0.927 (sep) | 0.739 / 0.826 (sep) | – | – |

What the table shows:

- **Reynolds number.** Aungier gives no Re, chord or P0 for this case. On his Mach numbers, both transition
  stations match for Re_c = 1.0e6 to 1.24e6. Our reading (fixed P0, T0 from Re 5.5e5 at M 0.30, so Re_c
  1.034e6 at M 0.65) lies at the low edge of that range. Re_c = 5.5e5 at M 0.65 contradicts the figure on
  both surfaces. This supports our reading; it does not confirm it.
- **c_f form.** The 2003 print (1 + Λ/6) contradicts the figure on both surfaces. The default
  (2 + Λ/6, Centrifugal 2000 Eq. 15-20) matches it.
- **F coefficient.** The cubic coefficient of F has no effect at these Λ.
- **Grid levels** (one way, all inside layout; SS transition f / x/c, separation f, PS transition f):

  | level | SS transition f / x/c | SS separation f | PS transition f |
  |---|---|---|---|
  | very_coarse | 0.150 / 0.226 | 0.900 | 0.250 |
  | coarse | 0.138 / 0.207 | 0.931 | 0.207 |
  | medium | 0.146 / 0.213 | 0.927 | 0.171 |
  | fine | 0.121 / 0.181 | 0.948 | 0.172 |
  | very_fine | 0.120 / 0.178 | 0.952 | 0.157 |

  The touching edge-node layout gives no boundary layer (`wall_arc_length` refuses it), and its inviscid
  solve goes sonic at the nose from coarse up.
- **Discretization of the laminar march** (a finding, not a defect). Marching the backward difference of
  Eq. 5-106 in ρ W² θ from θ = 0 gives θ² at the first step half the flat-plate value 2 C1 g μ x/(ρ W).
  This is a first-order start error. Refining the stations on the same edge distribution (our medium-grid
  Mach, interpolated linearly in f) moves the SS Re_θ = 250 crossing 0.135, 0.124, 0.115, 0.110, 0.106,
  0.105 for 24, 42, 83, 165, 329 and 657 stations. The squared form gives 0.104 to 0.102 on every
  station count. The separation is unaffected. Aungier's laminar H at his
  station 1 (SS 2.88, PS 3.21) needs Λ ≈ +4.7 and −6.6. Our march on his Mach numbers gives +1.9 and −3.6
  (H 2.99, 2.99), so his K is about 2 × ours there. That would fit a start without this error. Aungier
  only says "backward-difference approximations of the form of Eq. (5-94)" (2003 p. 113), so this is not
  confirmed. The transition and separation stations are the same either way.

**Two-way coupling, omega_BL** (medium grid, constant wake, tolerance 1e-4, up to 150 iterations):

| omega_BL | outcome |
|---|---|
| 0.1 | converges, monotone from below: β_out 23.57 → 25.87 deg, about 0.88 per iteration |
| 0.2 | converges: 25.87 deg (constant wake) / 26.47 deg (marched), about 0.77 per iteration |
| 0.3 | β_out 27.28, 24.74, then sonic at the SS TE node (91, 0) |
| 0.4 | 29.00, 26.63, 30.39, then sonic at node (89, 0) |
| 0.5 | 30.97, then sonic at node (91, 0) |
| 0.6 | sonic at node (91, 0) in the first coupled solve |

The other conditions behave the same way. β1 = 42 deg with omega 0.5 goes sonic at node (91, 0). M1 = 0.5
with omega 0.5 runs 28.64, 25.72, 35.76, then T ≤ 0. β1 = 40 deg fails already in the one-way layer, at
laminar station 2 (dW/W = −0.26: "no momentum thickness"); this is a separate open point.

From omega = 0.3 up, the iteration oscillates with growing amplitude. The oscillation shows only on the
separated suction-surface tail: the SS separation point and the loss alternate, while the pressure surface
and both transitions do not move. The failure is always at the suction-surface TE nodes.

A scratch probe (reverted, not committed) ran 25 iterations at omega = 0.1 and then switched to the value
under test. At the switch the change was 3.7e-2, close to the fixed point.

- omega 0.3: converged, in 34 iterations.
- omega 0.5: the change fell to 4.8e-4, then grew about 2 times per iteration (1.3e-3, 2.7e-3, ...,
  0.1 at k = 40). It failed at k = 49.
- omega 0.8 and 1.0: β_out alternated, the change grew 2 to 6 times per iteration, and the run failed at
  k = 35 and 34.

So the fixed point itself is unstable for omega ≥ 0.5, not only the first large step. The relaxed
iteration m ← m + omega (F(m) − m) has a mode with growth factor 1 − omega (1 − μ). From the growth at
0.5, 0.8 and 1.0, μ ≈ −5, so the limit is omega ≈ 2/(1 − μ) ≈ 0.3 (the independent check: about
0.33). At 0.3 that mode is barely damped (about
−0.8 per iteration). Without the ramp, the large first steps push it over. At 0.2 it is damped to about −0.2
per iteration: the early β_out alternation at omega 0.2 (25.79, 25.09, 25.74, 25.64 deg). The slow monotone
mode contracts at 0.88 (omega 0.1) and 0.77 (omega 0.2).

The mode lives on the separated tail; that it is caused by the separation is indicated, not proven. We
found no case without suction-surface separation to compare: β1 42 deg and M1 0.5 also separate, and β1
40 deg fails in the layer. A strongly negative μ (more displacement, more blockage, faster flow, much
thinner separated layer) fits the known limit of direct, edge-velocity-driven coupling near separation, for
which inverse and semi-inverse couplings exist. That is cited from general knowledge; no paper was opened.

**Defaults (changed 2026-11-24, the user's decision).** omega_BL = 0.2 and max_bl_iterations = 60. They
were 0.5 and 30. Both live in `CaseSettings` and `ViscousSettings`, and the window's Solver page takes
them from `CaseSettings`. At 0.2 the reference case converged in 35 iterations (constant wake) and 34
(marched), above the old limit of 30. These counts are coupled iterations after k = 0, which is what
max_bl_iterations counts; README.md's table counts k = 0 too (36 and 35). The independent check also
ran 0.2 at β1 42 deg (35 iterations, 25.70 deg) and at M1 0.5 (36 iterations, 25.75 deg). At 0.1 it took 74. The old default, 0.5, diverges on the reference case.
A test checks that this case converges with the defaults (`test_bl_validation.cpp`, `[slow]`).

The stability limit depends on the case. It comes from one eigenvalue estimate on one case, and 0.2 leaves
a margin of only about 1.6 below it. Another blade, incidence or Mach number can have a stronger
separated-tail response and a lower limit.

**Decisions (the user, 2026-11-24).**

1. **Station fraction.** The station-fraction reading of the axis is accepted. f is in
   `BoundaryLayerSurfaceData` (`f`, `transition_f`, `separation_f`, `laminar_separation_f`;
   `Viscous.hpp` `blade_station_fraction`) and in `bl_suction.csv` / `bl_pressure.csv` (the column after
   `x_c`). The Boundary layer window has a selectable abscissa (x/c or f).
2. **The inviscid SS study is later.** It is open item 1 below.
3. **omega_BL 0.2 and 60 BL iterations are the defaults.** The limit, about 0.3, depends on the case
   (above).
4. **The laminar start stays as it is.** It is the backward difference of Eq. 5-106 from θ = 0, with
   its first-order start error.
5. **The Re reading stays** (fixed P0, T0 from Re 5.5e5 at M 0.30).

**Open items.**

1. **Our inviscid SS distribution against Fig. 5-6** (a potential-flow study, for later).
   - Ruled out: the Kutta condition or the exit angle, and the station count (very_coarse with 21
     stations matches medium within 0.003).
   - The circulation check: Fig. 5-6's W_SS − W_PS over the interior stations gives Γ ≥ 12.65 m²/s,
     against s ΔW_u = 12.17 at 24.4 deg and 12.51 at our 23.57 deg. That is more circulation with less
     turning, which a potential flow of our blade, pitch and inlet cannot produce.
   - Candidates first:
     - the blade geometry: the camber-line distribution and Aungier's Ch. 4 representation;
     - the solidity or pitch: about 0.94 would reconcile them, a hypothesis only;
     - the inlet definition.
   - The LE: our SS/PS first-station Mach 0.77/0.34 against his 0.76/0.66. The stagnation point sits
     elsewhere.
2. **β1 = 40 deg, one way: the layer fails at laminar station 2.** The message is "no momentum
   thickness", with dW/W = −0.26 in one step. The backward-difference station equation has no root for
   such a step. Is it a refusal or a defect? Not yet decided.

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
