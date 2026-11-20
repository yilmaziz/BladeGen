"""Generate the BladeGen notebooks (python/notebooks/*.ipynb) from this script.

    python python/notebooks/build_notebooks.py

The notebooks are committed without outputs; ``python/run_notebooks.py`` executes
them headless (BLADEGEN_NB_STATIC=1) in ctest / CI.
"""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent

SETUP = r'''%matplotlib inline
import os, sys, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# the bladegen package lives in <build>/python (cmake --preset python): found on PYTHONPATH,
# else in $BLADEGEN_PY, else in the default build folder of this checkout
REPO = Path.cwd().resolve().parents[1]
try:
    import bladegen as bg
except ImportError:
    for cand in [os.environ.get("BLADEGEN_PY", ""), REPO / "build" / "python" / "python"]:
        if cand and Path(cand).is_dir():
            sys.path.insert(0, str(cand))
    import bladegen as bg
import bladegen.notebook as nb
OUT = Path(os.environ.get("BLADEGEN_TEST_OUTPUT", REPO / "build" / "notebook_output"))
OUT.mkdir(parents=True, exist_ok=True)
print("BladeGen", bg.__version__, "| widgets:", nb.HAVE_WIDGETS)'''


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


def write(name, cells):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nbf.write(nb, HERE / name)
    print("wrote", name)


# ─────────────────────────────────────────────────────────────────────────────
write("01_profile_designer.ipynb", [
    md("""# 01 · Profile designer — every knob under its family

The 2-D profile of BladeGen lives in the **unstaggered (m′, θ) plane**: LE point at (0, 0), TE point at
(1, 0), suction side (SS) above. It is built by the C++ engine exactly as the blade2d Python engine did
(DLR-style characteristic points, cubic B-spline SS, thickness-coupled PS), with **C2 / G2 connectors** at
the leading and trailing edges.

**What you see on the dashboard**

| panel | content |
|---|---|
| top | the profile: SS blue, PS orange, LE connector green, TE connector red; control polygons dashed; characteristic points black |
| zooms | LE and TE close up (where the connectors meet the sides) |
| thickness | measured normal thickness (SS normal → PS) and the thickness law t(ξ) |
| curvature | **curvature along SS and PS, from the LE point to the TE point** (the plot the Python engine was missing) |
| badge | the geometry check: **✗ red = ERROR**, ⚠ orange = warning — the geometry is still built and drawn |

Errors are marked where they occur (red ✗ on the profile and on the curvature plot)."""),
    code(SETUP),
    md("""## Interactive designer

The knobs are grouped by family — **Leading edge, Trailing edge, Suction side, Thickness distribution,
Pressure side, Connector** — each knob sits under its family heading. Every change rebuilds the profile in
C++ (≈ 1 ms) and redraws the dashboard. Try `LE connector tension` → 6 to see a loop flagged as an
error, or pull `SS CP2 theta` down until the pressure side crosses the suction side."""),
    code("knobs = nb.interactive_profile(bg.ProfileParams.naca65_1210())"),
    md("""## What the curvature check found in the Python engine's NACA 65 fit

`naca65_1210_blade2d()` is the **blade2d fit** (its `naca65_1210_blade_params.json`) with the blade2d
connector (`G2Legacy`): the geometry the Python engine produced. Its TE radius sits on the lower bound
(0.000255 chord) because NACA 65 has a *cusped* trailing edge; the TE connector then has to form an
**S-bend** a few micrometres long — invisible in a profile plot, obvious in the curvature plot, and flagged
as an error. `naca65_1210()` is BladeGen's own refit with a fairness penalty: same accuracy, no defect."""),
    code("""for p, t in [(bg.ProfileParams.naca65_1210_blade2d(), "blade2d fit + blade2d connector (Python geometry)"),
             (bg.ProfileParams.naca65_1210(), "BladeGen fairness-weighted refit (default)")]:
    fig, prof, q = nb.profile_dashboard(p, title=t)
    plt.show()
    print(q["summary"])"""),
    md("""## The three connector modes

* **Fair** (default): G2 at both junctions (position, unit tangent and curvature equal to the side's) and
  the *fairest* such curve — a degree-5 B-spline with 4 interior knots through the edge point that
  minimises the jerk energy ∫|C‴|² (one KKT solve).
* **C2Parametric**: copies the sides' first *and second* derivatives, so the closed loop is one C2 curve
  in a single parameter — exact, but the copied tangential acceleration can make a small edge spike.
* **G2Legacy**: blade2d's two-span quintic (the Python engine).

Below: the LE curvature of the same profile with each mode."""),
    code("""fig, ax = plt.subplots(1, 2, figsize=(14, 4.5))
base = bg.ProfileParams.naca65_1210()
for mode, ls in [("Fair", "-"), ("C2Parametric", "--"), ("G2Legacy", ":")]:
    p = base.copy(); p.connector = getattr(bg.ConnectorMode, mode)
    if mode != "Fair":
        p.le.tension, p.te.tension = 0.28, 0.20
    q = bg.check_quality(bg.Profile2D(p))
    for k, side, c in ((0, "ss", nb.C_SS), (1, "ps", nb.C_PS)):
        s = q[side]; x = s["s"] / s["s"][-1]
        ax[k].plot(x, s["kappa"], ls, color=c, label=f"{mode}: {q['errors']} err / {q['warnings']} warn")
for a, t in zip(ax, ["suction side", "pressure side"]):
    a.set_xlim(0, 0.05); a.set_yscale("symlog", linthresh=1); a.grid(alpha=.3); a.legend(fontsize=8)
    a.set_title(f"{t}: curvature near the LE"); a.set_xlabel("arc length / side length")
plt.show()"""),
    md("""## Saving and loading

Profiles are flat JSON files. `read_profile_json` also reads blade2d's `<name>_blade_params.json`
directly (its keys are recognised), so every profile fitted with the Python tools can be loaded."""),
    code("""p = knobs.params
bg.write_profile_json(p, str(OUT / "my_profile.json"))
q = bg.read_profile_json(str(OUT / "my_profile.json"))
old = bg.read_profile_json(str(REPO / "tests" / "data" / "naca65_1210_blade_params.json"))  # a blade2d file
print("round trip ok:", all(abs(p.get(k) - q.get(k)) < 1e-12 for k in p.to_dict() if isinstance(p.to_dict()[k], float)))
print("blade2d file -> connector", old.connector, "| LE radius", old.le.radius)"""),
])

# ─────────────────────────────────────────────────────────────────────────────
write("02_blade_3d.ipynb", [
    md("""# 02 · 3-D blade — C2 NURBS skin, solid, hub fillet

A blade is a set of **stations** (span, profile, placement) on stream surfaces of a **meridional channel**.
Between stations every profile parameter is interpolated in span (natural cubic splines), each section is
placed in the conformal plane (m′, θ) of its stream line (angles are the true flow angles) and mapped to 3-D
(x axial, y = r sin θ, z = r cos θ).

OpenCASCADE then builds:

1. the **skin**: every section interpolated by a *periodic C2* cubic B-spline with common parameters, the
   poles interpolated spanwise (Piegl & Tiller skinning) → one B-spline surface, periodic around the profile,
   C2 in both directions. (OCCT's own periodic interpolation puts a double knot at the seam — only C1 —
   so BladeGen solves the periodic system itself.)
2. the **solid**: the skin, extended past hub and casing, splits an annulus sector → caps exactly on the hub
   and casing surfaces of revolution;
3. on request the **hub fillet**: a rolling-ball fillet built face by face (OCCT's BRepFilletAPI fails on
   this junction)."""),
    code(SETUP),
    code("""blade = bg.Blade(bg.BladeDefinition.demo_rotor())
secs = blade.sections(30, 100)
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
nb.meridional_view(ax[0], blade, secs); ax[0].set_title("meridional view: channel and sections")
for s in secs:
    ax[1].plot(*s["conformal"].T, lw=1, color=plt.cm.viridis(min(max(s["span"], 0), 1)))
ax[1].set_aspect("equal"); ax[1].grid(alpha=.3); ax[1].set_xlabel("m'"); ax[1].set_ylabel("θ [rad]")
ax[1].set_title("sections in their conformal planes (colour = span)")
plt.show()"""),
    md("""## Station editor

Pick a station; its profile knobs (grouped by family) and placement knobs change only that station. The
plot shows the section at that span with its quality badge and all sections in 3-D."""),
    code("""definition = bg.BladeDefinition.demo_rotor()

def station_view(i=1):
    st = definition.stations[i]
    fig, prof, q = nb.profile_dashboard(st.profile, title=f"station {i} (span {st.span:.2f})")
    plt.show()

if nb.HAVE_WIDGETS:
    import ipywidgets as W
    from IPython.display import display
    pick = W.Dropdown(options=[(f"station {i}: span {s.span:.2f}", i) for i, s in enumerate(definition.stations)],
                      value=1, description="station")
    pl = {k: W.FloatSlider(value=getattr(definition.stations[1].placement, k), min=lo, max=hi, step=(hi - lo) / 400,
                           description=k, continuous_update=False)
          for k, lo, hi in [("stagger_deg", -10, 75), ("m_le", 0.05, 0.6), ("m_te", 0.4, 0.95), ("theta_le_deg", -10, 10)]}
    out = W.Output()
    knobs = nb.ProfileKnobs(definition.stations[1].profile)
    def refresh(*_):
        i = pick.value
        st = definition.stations[i]
        st.profile = knobs.params.copy()
        p = st.placement
        for k, s in pl.items():
            setattr(p, k, s.value)
        st.placement = p
        stations = list(definition.stations); stations[i] = st; definition.stations = stations
        with out:
            out.clear_output(wait=True)
            b = bg.Blade(definition)
            fig = plt.figure(figsize=(15, 5))
            a = fig.add_subplot(1, 2, 1)
            sec = b.section(st.span, 30, 100)
            q = bg.check_quality(bg.Profile2D(st.profile))
            a.plot(*sec["conformal"].T, color="k"); a.set_aspect("equal"); a.grid(alpha=.3)
            a.set_title(f"station {i}: {q['summary'].splitlines()[0]}", color="red" if q["errors"] else "green")
            a3 = fig.add_subplot(1, 2, 2, projection="3d")
            for s in b.sections(20, 60):
                a3.plot(*(s["xyz"].T * 1000), color="steelblue", lw=0.8)
            a3.set_box_aspect((1, 1, 1)); a3.set_title("all sections [mm]")
            plt.show()
    def on_pick(change):
        st = definition.stations[pick.value]
        knobs.set_params(st.profile)
        for k, s in pl.items():
            s.value = getattr(st.placement, k)
    pick.observe(on_pick, names="value")
    knobs.observe(lambda p: refresh())
    for s in pl.values():
        s.observe(refresh, names="value")
    display(W.VBox([W.HBox([pick] + list(pl.values())), knobs.widget, out]))
    refresh()
else:
    station_view(1)"""),
    md("""## The NURBS skin: degree, periodicity, C2 across the seam, curvature

`seam_jumps` evaluates the first and second u-derivatives on both sides of the periodic seam: relative
jumps ≈ 1e-10 mean C2. The mean-curvature map shows the surface fairness over (u, v)."""),
    code("""blade = bg.Blade(definition)
shapes = bg.build_blade(blade, 0.0, 40, 160)
S = shapes.surface
print(f"skin: degree ({S.u_degree}, {S.v_degree}), poles {S.n_poles}, u-periodic {S.u_periodic}, OCCT continuity {S.continuity}")
print(f"interpolation error at the section points: {shapes.skin_error:.2e} m")
for v in (0.1, 0.5, 0.9):
    gap, d1, d2 = S.seam_jumps(v)
    print(f"seam at v = {v}: gap {gap:.1e} m, |Δ d1| / |d1| = {d1:.1e}, |Δ d2| / |d2| = {d2:.1e}")
print(f"blade solid: valid {shapes.blade_valid}, volume {1e6 * shapes.blade_volume:.2f} cm³, area {1e4 * shapes.blade_area:.1f} cm²")
K, H = S.curvatures(240, 60)
fig, ax = plt.subplots(figsize=(14, 3.5))
im = ax.imshow(H.T, origin="lower", aspect="auto", extent=[0, 1, 0, 1], cmap="coolwarm",
               vmin=-np.nanpercentile(abs(H), 98), vmax=np.nanpercentile(abs(H), 98))
plt.colorbar(im, label="mean curvature [1/m]"); ax.set_xlabel("u (around the profile, from the LE point)"); ax.set_ylabel("v (span)")
ax.set_title("mean curvature of the C2 skin"); plt.show()"""),
    md("""## 3-D view with the **Add hub fillet** button

The button builds the rolling-ball fillet (radius in mm) in C++: ball centres on the hub contour offset by r
and on the blade surface offset by r along its fluid-side normal, circular arcs between the two contact
curves, a periodic C2 fillet surface, and a valid B-rep with exact shared edges. **Export** writes STEP and STL
(mm) of the blade and of the hub + blade (+ fillet)."""),
    code("""state = nb.fillet_panel(bg.Blade(definition), radius_mm=3.0, out_dir=str(OUT), name="demo_rotor")
if not nb.HAVE_WIDGETS:
    fig = state if not isinstance(state, dict) else None
    if fig is not None:
        try:
            fig.show()
        except Exception:
            plt.show()"""),
])

# ─────────────────────────────────────────────────────────────────────────────
write("03_inverse_fit.ipynb", [
    md("""# 03 · Inverse design — profile parameters from published coordinates

Port of blade2d's `tools/profile_inverse_fit.py` to C++: **side-aware residuals** (suction-side points are
measured against LE connector + SS + TE connector, pressure-side points against LE connector + PS + TE
connector), a **bounded Levenberg–Marquardt** solver (Eigen) in stages, and an optional **fairness penalty**
on concave dips of the edge-connector curvature, so a fit cannot hide an S-shaped edge.

Target: the published NACA 65-(12)10 section (blade2d `naca65/naca65_1210_XY_dense241.csv`, built from
NACA Table 6 with C_Li = 1.2, t/c = 0.10). blade2d reached RMS 0.0095 % chord on it."""),
    code(SETUP),
    code("""ss, ps = nb.read_sides_csv(str(REPO / "tests" / "data" / "naca65_1210_XY_dense241.csv"))
start = bg.ProfileParams.cda_start()            # a generic start, far from the answer
results = {}
for label, mode, w in [("blade2d connector, no fairness", "G2Legacy", 0.0),
                       ("Fair connector, no fairness", "Fair", 0.0),
                       ("Fair connector + fairness penalty", "Fair", 1.0)]:
    s0 = start.copy(); s0.connector = getattr(bg.ConnectorMode, mode)
    if mode == "G2Legacy":
        s0.le.tension, s0.te.tension = 0.28, 0.20
    opt = bg.FitOptions(); opt.params = nb.naca_search_box(s0); opt.fairness_weight = w
    t0 = time.perf_counter(); r = bg.fit_profile(ss, ps, s0, opt); dt = time.perf_counter() - t0
    q = bg.check_quality(bg.Profile2D(r["params"]))
    results[label] = (r, q)
    print(f"{label:<36} RMS {100 * r['rms']:.4f} % chord, max {100 * r['max']:.4f} %, "
          f"{r['evaluations']} evaluations, {dt:.1f} s, {q['errors']} errors / {q['warnings']} warnings")"""),
    md("""The C++ fit needs a few seconds where the Python fit needed minutes. Below: data against the fitted
profile, the LE / TE zooms and the residual of every published point."""),
    code("""r, q = results["Fair connector + fairness penalty"]
prof = bg.Profile2D(r["params"])
fig = plt.figure(figsize=(15, 9)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])
ax = fig.add_subplot(gs[0, :]); nb.draw_profile(ax, prof, q, control=False)
ax.plot(*ss.T, "o", mfc="none", mec="k", ms=3.5, label="published"); ax.plot(*ps.T, "o", mfc="none", mec="k", ms=3.5)
nb._badge(ax, q); ax.legend(fontsize=8, loc="lower right")
for k, (c, lo, hi) in enumerate([("LE", -0.01, 0.05), ("TE", 0.95, 1.01)]):
    a = fig.add_subplot(gs[1, k]); nb.draw_profile(a, prof, q, control=False)
    a.plot(*ss.T, "o", mfc="none", mec="k", ms=4); a.plot(*ps.T, "o", mfc="none", mec="k", ms=4)
    a.set_xlim(lo, hi); a.set_ylim(-0.03, 0.04) if c == "LE" else a.set_ylim(-0.02, 0.02); a.set_title(c + " zoom")
a = fig.add_subplot(gs[1, 2])
a.vlines(np.arange(len(r["d_ss"])), 0, 100 * r["d_ss"], color=nb.C_SS, label="suction side")
a.vlines(len(r["d_ss"]) + 2 + np.arange(len(r["d_ps"])), 0, 100 * r["d_ps"], color=nb.C_PS, label="pressure side")
a.set_ylabel("distance [% chord]"); a.set_title(f"residuals (RMS {100 * r['rms']:.4f} %)"); a.legend(fontsize=8)
plt.show()"""),
    md("""## Curvature quality of the three fits

Same accuracy, different geometry quality: without the penalty the fit reproduces the cusped NACA TE with an
S-shaped TE connector (error); with it the TE stays convex."""),
    code("""fig, ax = plt.subplots(1, 3, figsize=(16, 4))
for a, (label, (r, q)) in zip(ax, results.items()):
    nb.draw_curvature(a, q, bg.Profile2D(r["params"]))
    a.set_title(f"{label}\\n{q['errors']} errors, {q['warnings']} warnings", fontsize=9,
                color="red" if q["errors"] else "green")
plt.show()
bg.write_profile_json(results["Fair connector + fairness penalty"][0]["params"], str(OUT / "naca65_1210_fit.json"))"""),
])

# ─────────────────────────────────────────────────────────────────────────────
write("04_rotor67.ipynb", [
    md("""# 04 · NASA Rotor 67 rebuilt with BladeGen

The 14 blade sections and the flowpath of NASA Rotor 67 (Strazisar et al., NASA TP-2879, 1989; tables in
`data/rotor67`) are

1. **constructed as published** (surface splines with the tabulated end angles, edge circles, tangency points —
   the C++ port of `rotor67_sections.py`),
2. **inverse-designed** section by section: every section is brought to the unstaggered frame of its stream
   surface and the BladeGen profile parameters are fitted (C++, in parallel threads),
3. **rebuilt in 3-D**: a C2 NURBS skin through the 14 fitted sections, trimmed on the hub and on the casing
   contour of surface 14 offset by the 1.016 mm tip clearance (TP-2879), and
4. **filleted** at the hub with the TP-2879 radius of 1.78 mm (button below)."""),
    code(SETUP),
    code("""data = bg.read_rotor67(str(REPO / "data" / "rotor67" / "Rotor_14_Sections_and_Flowpath_EN.txt"))
print(f"{len(data.sections)} sections, {data.n_blades} blades, tip clearance {1e3 * data.tip_clearance:.4f} mm")
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
ax[0].plot(*data.hub.T * 1000, "k", lw=1.5); ax[0].plot(*data.shroud.T * 1000, "k", lw=1.5)
for s in data.sections:
    ax[0].plot(*s.streamline.sample(60).T * 1000, color="0.5", lw=0.8)
ax[0].set_aspect("equal"); ax[0].set_xlim(-30, 120); ax[0].grid(alpha=.3); ax[0].set_title("flowpath and section stream lines [mm]")
for s in data.sections:
    ax[1].plot(*s.profile_m_theta().T, lw=0.9)
ax[1].set_xlabel("m [m]"); ax[1].set_ylabel("θ [rad]"); ax[1].grid(alpha=.3); ax[1].set_title("published sections in (m, θ)")
plt.show()"""),
    md("""## Fit all 14 sections (C++, threads)

The fit releases the GIL, so the 14 sections run in parallel. `fairness_weight` keeps the fitted edges free of
S-bends. The table reports the 2-D deviation (% chord) and the 3-D deviation of the published points from the
rebuilt sections (µm)."""),
    code("""from concurrent.futures import ThreadPoolExecutor
def fit(k):
    s = data.sections[k]
    opt = bg.FitOptions(); opt.fairness_weight = 0.3
    return bg.fit_section(s.profile_conformal(), s.streamline, k / 13.0, bg.ProfileParams.cda_start(), opt)
t0 = time.perf_counter()
with ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
    fits = list(ex.map(fit, range(len(data.sections))))
print(f"14 sections fitted in {time.perf_counter() - t0:.0f} s\\n")
print(" sec  stagger  RMS 2D [%c]  max 2D [%c]   RMS 3D [um]  max 3D [um]   quality")
for s, f in zip(data.sections, fits):
    q = bg.check_quality(bg.Profile2D(f["fit"]["params"]))
    print(f" {s.index:3d}  {np.degrees(f['frame']['stagger']):6.2f}   {100 * f['fit']['rms']:9.4f}   {100 * f['fit']['max']:9.4f}"
          f"   {1e6 * f['rms_3d']:10.2f}  {1e6 * f['max_3d']:10.2f}   {q['errors']} err / {q['warnings']} warn")"""),
    md("""## Published vs rebuilt sections (hub, mid, tip) and their curvature"""),
    code("""fig = plt.figure(figsize=(17, 12))
gs = fig.add_gridspec(3, 3, width_ratios=[2.6, 1.0, 1.5], hspace=0.45, wspace=0.25)
for row, k in enumerate([0, 6, 13]):
    f = fits[k]; prof = bg.Profile2D(f["fit"]["params"]); q = bg.check_quality(prof)
    fr = f["frame"]; color = "red" if q["errors"] else ("darkorange" if q["warnings"] else "green")
    a = fig.add_subplot(gs[row, 0]); nb.draw_profile(a, prof, q, control=False, points=False, adjustable="box")
    a.plot(*fr["ss"].T, "o", mfc="none", mec="k", ms=2.5); a.plot(*fr["ps"].T, "o", mfc="none", mec="k", ms=2.5)
    a.set_title(f"section {data.sections[k].index}: published (o) vs BladeGen — RMS {1e6 * f['rms_3d']:.1f} µm; "
                f"{q['errors']} errors, {q['warnings']} warnings", fontsize=9, color=color)
    a = fig.add_subplot(gs[row, 1]); nb.draw_profile(a, prof, q, control=False, points=False, adjustable="box")
    a.plot(*fr["ss"].T, "o", mfc="none", mec="k", ms=4); a.plot(*fr["ps"].T, "o", mfc="none", mec="k", ms=4)
    r = f["fit"]["params"].le.radius; a.set_xlim(-2 * r, 8 * r); a.set_ylim(-5 * r, 5 * r); a.set_title("LE zoom", fontsize=9)
    nb.draw_curvature(fig.add_subplot(gs[row, 2]), q, prof)
plt.show()
for k in (0, 6, 13):
    print(f"section {data.sections[k].index}:", bg.check_quality(bg.Profile2D(fits[k]["fit"]["params"]))["summary"])"""),
    md("""## The 3-D blade: C2 NURBS skin through the 14 fitted sections

The definition is written to JSON (exact B-spline stream lines) so the CLI (`bladegen build`) reproduces the
same blade."""),
    code("""definition = bg.BladeDefinition()
definition.name = "NASA_Rotor_67"; definition.n_blades = data.n_blades; definition.tip_clearance = data.tip_clearance
definition.stations = [f["station"] for f in fits]
definition.channel = bg.MeridionalChannel(data.hub_curve(), data.shroud_curve())
bg.write_blade_json(definition, str(OUT / "rotor67_bladegen.json"))
blade = bg.Blade(definition)
loops = [s["xyz"] for s in blade.sections(60, 200)]
inside = None
mid = blade.section(0.5, 60, 200)["xyz"]; inside = 0.5 * (mid[len(mid) // 4] + mid[3 * len(mid) // 4])
t0 = time.perf_counter()
shapes = bg.build_shapes(loops, data.hub_curve(), data.sections[-1].streamline, inside, data.tip_clearance, 0.0, 1.0)
print(f"skin + solid in {time.perf_counter() - t0:.1f} s: valid {shapes.blade_valid}, volume {1e6 * shapes.blade_volume:.3f} cm³")
S = shapes.surface
print(f"skin degree ({S.u_degree},{S.v_degree}), poles {S.n_poles}, seam jumps (gap, d1, d2) = {tuple(f'{x:.1e}' for x in S.seam_jumps(0.5))}")
pub = np.vstack([f["published_xyz"][::6] for f in fits])
d = shapes.blade_face.distances(pub)
print(f"published points → NURBS skin: RMS {1e6 * np.sqrt(np.mean(d ** 2)):.2f} um, max {1e6 * d.max():.2f} um ({len(d)} points)")"""),
    code("""fig = nb.show_shapes([(shapes.blade, "lightsteelblue", "Rotor 67 blade (BladeGen)"), (shapes.hub, "0.8", "hub sector")],
                     title="NASA Rotor 67 — C2 NURBS blade, tip clearance 1.016 mm")
try:
    import plotly.graph_objects as go
    fig.add_trace(go.Scatter3d(x=pub[:, 0] * 1000, y=pub[:, 1] * 1000, z=pub[:, 2] * 1000, mode="markers",
                               marker=dict(size=1.5, color="black"), name="published points"))
    fig.show()
except ImportError:
    plt.show()"""),
    md("""## Hub fillet, r = 1.78 mm (TP-2879) — and export

Press **Add hub fillet**. The rolling-ball construction reports its own checks: the ball touches both walls at
exactly r, and the fillet is tangent to both (G1 angle)."""),
    code("""sf = bg.build_shapes(loops, data.hub_curve(), data.sections[-1].streamline, inside, data.tip_clearance, 1.78e-3, 1.0)
print(sf.fillet_message, f"| G1 {sf.fillet_g1_deg:.3f}°, ball error {sf.fillet_ball_error:.1e} m")
if sf.fillet_ok:
    sf.blade_on_hub.write(str(OUT / "rotor67_hub_blade_fillet.step"))
    sf.blade_with_fillet.write(str(OUT / "rotor67_blade_fillet.step"))
shapes.blade.write(str(OUT / "rotor67_blade.step"))
fig = nb.show_shapes([(sf.blade_on_hub if sf.fillet_ok else shapes.blade, "lightsteelblue", "hub + blade + fillet")],
                     title="Rotor 67 with the 1.78 mm hub fillet")
try:
    fig.show()
except Exception:
    plt.show()
if nb.HAVE_WIDGETS:
    import ipywidgets as W
    from IPython.display import display
    radius = W.BoundedFloatText(value=1.78, min=0.2, max=6.0, step=0.01, description="fillet [mm]")
    btn = W.Button(description="Add hub fillet", button_style="primary", icon="magic")
    msg = W.HTML(); view = W.Output()
    def add(_):
        msg.value = "<i>building ...</i>"
        s2 = bg.build_shapes(loops, data.hub_curve(), data.sections[-1].streamline, inside, data.tip_clearance,
                             radius.value * 1e-3, 1.0)
        msg.value = (f"<b style='color:{'green' if s2.fillet_ok else 'red'}'>{s2.fillet_message}</b> "
                     f"G1 {s2.fillet_g1_deg:.3f}°")
        if s2.fillet_ok:
            s2.blade_on_hub.write(str(OUT / f"rotor67_hub_blade_fillet_{radius.value:.2f}mm.step"))
            with view:
                view.clear_output(wait=True)
                nb.show_shapes([(s2.blade_on_hub, "lightsteelblue", "hub + blade + fillet")]).show()
    btn.on_click(add)
    display(W.VBox([W.HBox([radius, btn]), msg, view]))"""),
    md("""## Compare with the STEP files of the Python scripts (if present)

`newblade/rotor67/rotor67_blade.step` was written by `rotor67_passage_solid.py` (blade to the casing, no tip
clearance by default); volumes are compared for orientation only."""),
    code("""py_step = Path(os.environ.get("ROTOR67_PY_STEP", r"C:\\Users\\ASUS\\Documents\\newblade\\rotor67\\rotor67_blade.step"))
if py_step.is_file():
    py = bg.read_step(str(py_step))
    print(f"Python blade (no tip clearance): {1e6 * py.volume:.3f} cm³   BladeGen (tip clearance 1.016 mm): {1e6 * shapes.blade_volume:.3f} cm³")
else:
    print("Python STEP not found; skipped")"""),
])
