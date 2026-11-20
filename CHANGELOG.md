# Changelog

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
