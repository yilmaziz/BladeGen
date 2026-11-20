# Changelog

## v1.0.1 — NACA 65 series and the section designer window

The first integration step of the BladeGen → ENIGMAgrid → TURBOlenza chain: TURBOlenza's blade geometry, the C2
NACA 65 cascade section, now lives in BladeGen, and a desktop window designs hub / mid / tip sections with live
curvature plots.

### New

**Parametric NACA 65 series** (`src/naca/Naca65`, migrated from TURBOlenza `b2b/geometry.py`)
- NACA 65-(Cl0 A10) t/c for Cl0 0–2.5 and t/c 0.02–0.20, both variants: **Kovach–Sandercock** (taper from
  station 0.60 to an R = 0.008 c TE circle) and **Emery** (sharp TE).
- Reviewed before porting. **Kept:**
  - Table 4-1 with the 2.5 % misprint corrected;
  - the closed-form A10 mean line and the arc-length thickness station;
  - the C4 quintic thickness in √station;
  - the nose circle;
  - the C2 **monotone-curvature leading edge** (Whewell-form fair, the curvature falls from 1/R_LE without
    overshoot or S-bend);
  - the equivalent-arc angles.

  **Dropped:**
  - the PCHIP thickness (curvature sawtooth);
  - the legacy LE circle with its polar fair;
  - kernel smoothing and the least-squares C2 refit.

  **Replaced by BladeGen's own constructions:**
  - a C2 blend at the Kovach start;
  - a G2 minimum-jerk TE connector through the TE circle's point with its tangent and curvature;
  - a C2 closing piece to the sharp Emery TE.
- Agreement with TURBOlenza's frozen C2 geometry on five reference sections:

  | variant | body | LE | TE | difference at the TE |
  |---|---|---|---|---|
  | Emery | ≤ 5·10⁻⁸ c | ≤ 4·10⁻⁸ c | ≤ 1.1·10⁻⁵ c | |
  | Kovach | ≤ 3·10⁻⁶ c | ≤ 5·10⁻⁸ c | ≤ 1.9·10⁻⁴ c | the different, G2 TE connector |
- About 0.25 s per section including the quality check (the same checks as custom profiles; the sharp Emery
  TE is reported as information).
- `math/Constrained.hpp`: Lawson–Hanson NNLS (with the anti-cycling safeguard; without it the LE fair took 7 s),
  least-distance programming, LSEI, Gauss–Legendre of any order.

**NACA 65 + Custom** (`fit/NacaFit`, `fit_naca65`)
- The custom BladeGen parameters are inverse-fitted to a NACA 65 section, starting from the 65-(12)10 fit scaled
  to Cl0 and t. The result is an ordinary `ProfileParams`, so every custom knob then edits the NACA geometry.
- Accuracy and time: RMS ≈ 0.011 % c for Kovach sections in 3–5 s.

**Section designer window** (`python -m bladegen.gui`; `tools/gui.cmd`, `tools/gui.sh`)
- Hub / mid / tip stations, each with its own construction method: **NACA 65**, **Custom** or **NACA 65 +
  Custom** (DCA and C4 shown as planned).
- Knobs grouped by family. The profile is shown with LE / TE zooms, the **curvature of the suction and
  pressure side** from the LE point to the TE point with LE / TE zooms, and a cascade view (stagger,
  solidity). Everything is redrawn on every change in a worker thread.
- The geometry is always built. S-shaped or looping LE / TE connectors and self-intersections are **errors**:
  red ✗ marks on the profile and the curvature, a red badge and the title "NOT MANUFACTURABLE AS DESIGNED".
- Overlay of all stations; save / open a design (JSON); export coordinates with curvature (CSV), **TURBOlenza
  CascadeSpec** entries for the NACA 65 stations (checked with `CascadeSpec(**entry)` and `make_cascade`), and
  the figure.
- Runs wherever the repository and its conda environment are: Tkinter + matplotlib, nothing else to install.
  The launchers build the Python module when it is missing or older than the GUI. High-DPI aware on Windows.
- By design there is no flowpath, span or 3-D in the window; the 3-D blade stays in the notebooks.

**Other**
- CLI: `bladegen naca65 <cl0> <t> [--variant …] [--csv sides.csv] [--fit custom.json]`; `bladegen version` now
  reports the CMake project version.
- Python: `Naca65Params`, `Naca65Section` (metadata, samples, camber line), `check_quality(section)`,
  `naca65_start`, `fit_naca65`.
- Tests:
  - Catch2 (NACA 65 and the solvers): Table 4-1, A10, the thickness law, the LE fair (monotone, G2), both TEs,
    the quality check over a Cl0 × t × variant grid, the TURBOlenza reference, NACA 65 + Custom, the solvers,
    invalid input.
  - pytest: the bindings, every GUI method, and the window driven headless.
  - CLI tests.
  - Linux CI runs the GUI test under Xvfb.

### Known limitations (v1.0.1)
- **Kovach station 0.60:** the table/taper corner (0.08°) is blended over ±0.005 in station. This keeps the
  surface within 3·10⁻⁶ c of TURBOlenza's, but the curvature has a short excursion there (TURBOlenza's frozen
  C2 surface has the same, κ −0.8…1.9 on the suction side), reported as suction-side inflection warnings. A
  wider blend removes it but moves away from TURBOlenza's geometry.
- **Pressure-side waviness:** the curvature of the thin NACA 65 pressure side changes sign a few times. This is
  inherited from the tabulated thickness and reported as a warning, not an error.
- **Emery in NACA 65 + Custom:** a sharp Emery TE cannot be matched by the round custom TE (RMS ≈ 0.07 % c,
  locally up to 0.36 % c at the TE).
- **The window:**
  - sections are 2-D (cascade plane), so there are no flowpath, span or 3-D in it;
  - a design has exactly three stations;
  - the TURBOlenza export covers NACA 65 stations only (TURBOlenza has no custom-profile input yet);
  - the fit runs 3–5 s and is not cancellable.
- **DCA and C4** profiles are not implemented yet (listed as planned in the window).
- All v1.0.0 limitations still apply.

### Next steps
- DCA and C4 sections (the next construction methods in the window), and the NACA 65 section as a 3-D
  station profile.
- An option for a curvature-fair Kovach start (wider blend) and a sharp / cusped TE in the custom profile, so
  NACA 65 + Custom can match the Emery TE.
- Export of stations in ENIGMAgrid's geometry format, and a TURBOlenza run from the window (solver option).

## v1.0.0 — first release

The C++23 / OpenCASCADE successor of the `blade2d` Python blade generator: a parametric C2 blade engine
with curvature-quality checks, inverse design from published coordinates, C2 NURBS solids, a rolling-ball
hub fillet, Python bindings and interactive notebooks. Tested on Windows, Linux and macOS.

### What is in it

**2-D profile engine** (`src/profile`)
- Port of the blade2d construction: DLR-style characteristic points at LE and TE, cubic B-spline suction
  side, pressure side coupled to it through a thickness law. The Python-compatible mode (`G2Legacy`)
  reproduces blade2d to 2·10⁻⁶ chord (regression test against the Python engine).
- New default LE / TE connector `Fair`: curvature-continuous (G2) at both junctions and the fairest such
  curve (minimum jerk energy, one KKT solve). Also `C2Parametric` (one parametrically C2 loop).
- 26 scalar parameters in six families (leading edge, trailing edge, suction side, thickness distribution,
  pressure side, connector) with generic access, spanwise blending and JSON I/O; blade2d
  `*_blade_params.json` files load directly.

**Geometry quality check** (`check_quality`) — the curvature plots the Python engine lacked
- Curvature of both sides from the LE point to the TE point; S-shaped / looping edges, curvature dips and
  double peaks, inflections, curvature oscillation, self-intersection, negative thickness, orientation,
  junction continuity. Problems are reported as errors / warnings with their location. The geometry is
  always built.

**3-D blades** (`src/blade`, `src/occt`)
- Meridional channel, exact arc length m and conformal m′ = ∫ dm / r (closed-form straight extrapolation
  past the ends); conformal placement with stagger, lean and mirror; natural-spline interpolation of every
  parameter in span.
- C2 NURBS skin (Piegl–Tiller skinning) that is **C2 across the periodic seam**: OCCT's periodic
  interpolation has a double knot there (C1). BladeGen solves the periodic system itself.
- Blade solid trimmed exactly on the hub and on the casing offset by the tip clearance; STEP AP214, IGES,
  STL, BREP export; tessellation for viewing.
- **Rolling-ball hub fillet** with an explicit B-rep (`BRepFilletAPI` fails on this junction). Ball error
  ~10⁻¹⁶ m, tangency 0.2–0.4°.

**Inverse design** (`src/fit`)
- Side-aware residuals, bounded Levenberg–Marquardt (Eigen), optional fairness penalty on edge defects.
- NACA 65-(12)10: RMS 0.0095 % chord in about 5 s.
- **NASA Rotor 67**: all 14 published sections (TP-2879) rebuilt on their stream surfaces with 3-D RMS
  1.6–16 µm and no geometry errors; C2 blade, 1.016 mm tip clearance, 1.78 mm hub fillet, STEP. The fitted
  definition is in `data/rotor67/rotor67_bladegen.json`.

**Tools**
- `bladegen` CLI (profile check, build, fit, demo) and `rotor67` (the whole Rotor 67 reconstruction).
- Python module `bladegen` (pybind11) and four notebooks: profile designer with knobs grouped by family and
  live curvature plots, 3-D blade with an **Add hub fillet** button, inverse fit, Rotor 67.
- `tools/jupyter.cmd` / `tools/jupyter.sh`: one command to build the module and open the notebooks.

**Quality assurance**
- 56 tests (Catch2, pytest, notebooks executed headless, CLI), each against a closed form, a published value
  or its defining equation, plus invalid input.
- CI on Windows (MSVC, OCCT 7.9), Linux (GCC, OCCT 7.6 from apt and an ASan/UBSan job) and macOS (AppleClang,
  OCCT 7.9). Zero warnings with warnings-as-errors everywhere.

### Findings on the Python engine and on OCCT
- The blade2d NACA 65-(12)10 fit has an S-shaped TE connector a few µm long. It is caused by fitting a round
  TE to the cusped NACA 65 trailing edge and is visible only in the curvature plot. The default
  `naca65_1210()` is a fairness-weighted refit with no geometry errors; the blade2d fit remains available as
  `naca65_1210_blade2d()`.
- blade2d's `enforce_c1` makes the thickness slopes at the peak equal in magnitude only; the law keeps a
  corner there (kept for compatibility, documented).
- OCCT periodic interpolation is C1 at the seam (this affects the STEP files of the Python Rotor 67
  scripts); `BRepFilletAPI` fails on the blade/hub edge.

### Known limitations
- **Fillet:** constant radius, hub only (no tip or casing fillet, no variable radius). The hub contour's
  curvature radius must exceed the fillet radius. Tangency to the walls is 0.2–0.4°, not exact, because the
  arcs are interpolated. About 5–10 s per blade.
- **Profile model:** the edges are always rounded, so cusped or sharp trailing edges (NACA 65) can only be
  approximated by a small round edge. The pressure side is shaped through three thickness stations. The
  power-law thickness distribution has a corner at its peak (it only places PS control points, so the PS
  itself stays smooth).
- **Continuity:** the profile is G2 (curvature-continuous) at the four junctions in the default mode;
  parametric C2 of the whole loop requires `C2Parametric`. The 3-D skin interpolates sampled section points,
  so it reproduces the analytic sections to interpolation accuracy (sub-µm), not exactly.
- **Inverse fit:** a local optimiser with a finite-difference Jacobian, so the result depends on the start
  point. No similarity transform is fitted, the PS stations are fixed, and a section takes 5–60 s.
- **Quality check:** the thresholds are heuristic. Circular published LEs give "double peak" warnings
  because a polynomial cannot have constant curvature.
- **Scope:** single blade rows. No splitter or tandem blades, squealer tips or blade-to-blade passage / fluid
  domain export (the Python scripts had one). Spanwise stacking is set by the station placements; there are
  no lean or sweep law curves.
- **Platforms:** macOS is tested in CI on Apple Silicon only. The notebooks need JupyterLab with ipywidgets.
- **Dependencies:** CppAD is not used (the fit's finite differences are fast enough at 24 parameters).

### Next steps
- **ENIGMAgrid coupling:** export section loops in (m′, θ), stream lines and the periodic passage boundaries
  directly in ENIGMAgrid's geometry format, so the mesher's blade-to-blade topologies read BladeGen output.
- **TURBOlenza as a solver option:** design → mesh (ENIGMAgrid) → flow solution (TURBOlenza) from one
  parameter set, then aerodynamic optimisation. Exact gradients (CppAD or an adjoint) become worthwhile at
  that point.
- Geometry: sharp / cusped TE and elliptic LE options, lean / sweep / bow stacking laws, variable-radius and
  tip fillets, passage and fluid-domain export with named faces, splitter blades.
- Inverse design: global start search, fitting the PS stations and a similarity transform, a full 3-D
  surface fit.

### Relation to ENIGMAgrid and TURBOlenza
BladeGen is the geometry stage of a three-part toolchain:

- **BladeGen** produces the blade: sections, NURBS skin, hub and casing, STEP / JSON.
- **ENIGMAgrid** (C++20 structured multi-block mesher) grids the blade-to-blade passage.
- **TURBOlenza** (flow solver, to be added as a solver option to these modules) computes the flow.

All three share the conformal (m′, θ) frame and the hub / casing contours, so a blade designed or
inverse-designed here can be meshed and analysed without a CAD round trip.
