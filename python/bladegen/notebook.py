"""Notebook helpers for BladeGen: family-grouped knobs, profile dashboard, 3-D views.

Every geometric quantity shown here comes from the C++ engine (``bladegen._core``):
this module only builds widgets and draws.  Conventions of every figure:

* suction side blue, pressure side orange, LE connector green, TE connector red;
* curvature plots run from the LE point to the TE point along each side
  (convexity-signed: > 0 where the surface bulges outward);
* geometry problems found by ``check_quality`` are drawn on the plots
  (ERROR: red cross, WARNING: orange triangle) and summarised in a badge; the
  geometry is still generated and drawn.

Set ``BLADEGEN_NB_STATIC=1`` to draw every figure once with its defaults (no
widgets): that is how the notebooks are executed headless in CI.
"""
from __future__ import annotations

import os
from typing import Callable

import numpy as np

from . import _core as bg

STATIC = os.environ.get("BLADEGEN_NB_STATIC") == "1"

C_SS, C_PS, C_LE, C_TE = "tab:blue", "tab:orange", "tab:green", "tab:red"
PIECE_COLOR = {0: C_LE, 1: C_SS, 2: C_TE, 3: C_PS}
FAMILY_ORDER = ["Leading edge", "Trailing edge", "Suction side", "Thickness distribution", "Pressure side"]
FAMILY_COLOR = {"Leading edge": C_LE, "Trailing edge": C_TE, "Suction side": C_SS,
                "Thickness distribution": "tab:purple", "Pressure side": C_PS,
                "Connector": "0.3", "Placement": "tab:brown"}

try:  # widgets are optional: without them every figure is drawn once
    import ipywidgets as W
    from IPython.display import display

    HAVE_WIDGETS = not STATIC
except ImportError:  # pragma: no cover
    W = None
    HAVE_WIDGETS = False

    def display(*a, **k):  # noqa: D401 - fallback
        pass


# ─────────────────────────────────────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────────────────────────────────────

def read_sides_csv(path, ss_cols=("upper_x", "upper_y"), ps_cols=("lower_x", "lower_y"), negate_y=False):
    """Suction / pressure side points from a '#'-commented CSV with a header line
    (blade2d profile_inverse_fit.load_target format; NaN-padded rows allowed)."""
    lines = [ln for ln in open(path, encoding="utf-8").read().splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    hdr = [h.strip() for h in lines[0].split(",")]
    rows = [[float(c) if c.strip() not in ("", "nan", "NaN") else np.nan for c in ln.split(",")] for ln in lines[1:]]
    tab = {h: np.array([r[i] for r in rows]) for i, h in enumerate(hdr)}

    def grab(cols):
        a = np.column_stack([tab[cols[0]], tab[cols[1]]])
        a = a[~np.isnan(a).any(axis=1)]
        if negate_y:
            a[:, 1] = -a[:, 1]
        return a[np.r_[True, np.any(np.diff(a, axis=0) != 0.0, axis=1)]]

    return grab(ss_cols), grab(ps_cols)


def naca_search_box(start):
    """The blade2d naca65/fit_naca65_1210.py parameter box (wider edges, cusped TE)."""
    box = {"le.radius": (0.0015, 0.025), "le.stretch_ab": (0.3, 4.0), "le.asym": (-0.008, 0.008),
           "le.beta_deg": (10.0, 60.0), "le.d_beta_ss_deg": (5.0, 75.0), "le.d_beta_ps_deg": (5.0, 80.0),
           "te.radius": (0.0002, 0.012), "te.beta_deg": (5.0, 40.0), "thickness.t_max": (0.05, 0.16),
           "thickness.t_te": (0.0005, 0.02), "ps.handle_le_frac": (0.002, 0.5), "ps.handle_te_frac": (0.005, 0.5)}
    spec = bg.default_fit_params(start)
    for p in spec:
        if p.name in box:
            p.lo, p.hi = box[p.name]
        p.x0 = min(max(start.get(p.name), p.lo), p.hi)
    return spec


# ─────────────────────────────────────────────────────────────────────────────
# Knobs
# ─────────────────────────────────────────────────────────────────────────────

def families() -> dict[str, list[dict]]:
    """The C++ parameter table grouped by family, in display order."""
    out: dict[str, list[dict]] = {f: [] for f in FAMILY_ORDER}
    for p in bg.param_table():
        out.setdefault(p["family"], []).append(p)
    return out


class ProfileKnobs:
    """Sliders for every profile parameter, each one under its family heading.

    ``knobs.params`` is the current :class:`bladegen.ProfileParams`;
    ``knobs.observe(fn)`` calls ``fn(params)`` after every change.
    """

    def __init__(self, params: "bg.ProfileParams", columns: int = 3, width: str = "330px"):
        self.params = params.copy()
        self._callbacks: list[Callable] = []
        self._mute = False
        self.sliders: dict[str, object] = {}
        if not HAVE_WIDGETS:
            self.widget = None
            return
        boxes = []
        for fam, items in families().items():
            color = FAMILY_COLOR.get(fam, "black")
            rows = [W.HTML(f"<b style='color:{color};font-size:13px'>{fam}</b>")]
            for it in items:
                lo, hi = it["lo"], it["hi"]
                v = self.params.get(it["name"])
                lo, hi = min(lo, v), max(hi, v)
                s = W.FloatSlider(value=v, min=lo, max=hi, step=(hi - lo) / 1000.0, description=it["label"],
                                  readout_format=".4g", continuous_update=False,
                                  style={"description_width": "165px"}, layout=W.Layout(width=width))
                s.observe(self._on_change(it["name"]), names="value")
                self.sliders[it["name"]] = s
                rows.append(s)
            if fam == "Pressure side":
                for i, mp in enumerate(self.params.ps.mprime_ctrl):
                    s = W.FloatSlider(value=mp, min=0.05, max=0.95, step=0.005, description=f"PS station {i + 1} m'",
                                      continuous_update=False, style={"description_width": "165px"},
                                      layout=W.Layout(width=width))
                    s.observe(self._on_station(i), names="value")
                    rows.append(s)
            box = W.VBox(rows, layout=W.Layout(border=f"1px solid {color}", padding="4px", margin="2px"))
            boxes.append(box)
        conn = W.Dropdown(options=["Fair", "C2Parametric", "G2Legacy"], value=str(self.params.connector).split(".")[-1],
                          description="connector", style={"description_width": "165px"}, layout=W.Layout(width=width))
        conn.observe(self._on_connector, names="value")
        reset = W.Button(description="reset", layout=W.Layout(width="90px"))
        reset.on_click(lambda _b: self.set_params(params))
        boxes.append(W.VBox([W.HTML("<b>Connector</b>"), conn, reset],
                            layout=W.Layout(border="1px solid gray", padding="4px", margin="2px")))
        self.connector = conn
        self.widget = W.GridBox(boxes, layout=W.Layout(grid_template_columns=f"repeat({columns}, {width})",
                                                       grid_gap="4px"))

    def _on_change(self, name):
        def cb(change):
            self.params.set(name, float(change["new"]))
            self._fire()
        return cb

    def _on_station(self, i):
        def cb(change):
            st = list(self.params.ps.mprime_ctrl)
            st[i] = float(change["new"])
            if all(a < b for a, b in zip(st, st[1:])):
                self.params.ps.mprime_ctrl = st
                self._fire()
        return cb

    def _on_connector(self, change):
        self.params.connector = getattr(bg.ConnectorMode, change["new"])
        self._fire()

    def _fire(self):
        if self._mute:
            return
        for cb in self._callbacks:
            cb(self.params)

    def set_params(self, p):
        """Move every slider to the values of `p` (one redraw at the end)."""
        self._mute = True
        try:
            for name, s in self.sliders.items():
                v = p.get(name)
                s.min, s.max = min(s.min, v), max(s.max, v)
                s.value = v
        finally:
            self._mute = False
        self.params = p.copy()
        self._fire()

    def observe(self, fn: Callable):
        self._callbacks.append(fn)


# ─────────────────────────────────────────────────────────────────────────────
# Profile dashboard
# ─────────────────────────────────────────────────────────────────────────────

def _badge(ax, q: dict, loc=(0.01, 0.98)):
    if q["errors"]:
        txt, fc = f"✗ {q['errors']} GEOMETRY ERROR(S), {q['warnings']} warning(s)", "#ffd6d6"
        ec = "red"
    elif q["warnings"]:
        txt, fc, ec = f"⚠ {q['warnings']} warning(s) — geometry built", "#fff1d6", "darkorange"
    else:
        txt, fc, ec = "✓ geometry OK (C2 / G2 junctions, no defects found)", "#dcf5dc", "green"
    lines = [txt] + [f"[{i['severity']}] {i['message'][:110]}" for i in q["issues"] if i["severity"] != "INFO"][:6]
    ax.text(*loc, "\n".join(lines), transform=ax.transAxes, va="top", ha="left", fontsize=8, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", fc=fc, ec=ec, lw=1.5), zorder=50)


def _issue_markers(ax, q: dict, where_key="where"):
    for i in q["issues"]:
        if i["severity"] == "ERROR":
            ax.plot(*i[where_key], "x", color="red", ms=12, mew=3, zorder=40)
        elif i["severity"] == "WARNING":
            ax.plot(*i[where_key], "^", color="darkorange", ms=9, mfc="none", mew=2, zorder=40)


def draw_profile(ax, prof: "bg.Profile2D", q: dict | None = None, control=True, points=True, lw=2.0):
    pcs = prof.pieces(300)
    ax.plot(*pcs["ss"].T, color=C_SS, lw=lw, label="suction side")
    ax.plot(*pcs["ps"].T, color=C_PS, lw=lw, label="pressure side")
    ax.plot(*pcs["le"].T, color=C_LE, lw=lw, label="LE connector")
    ax.plot(*pcs["te"].T, color=C_TE, lw=lw, label="TE connector")
    if control:
        cp = prof.control_polygons()
        for k, c in (("ss", C_SS), ("ps", C_PS), ("le", C_LE), ("te", C_TE)):
            ax.plot(*cp[k].T, "--o", color=c, lw=0.6, ms=2.5, alpha=0.45)
    if points:
        for name, p in prof.characteristic_points().items():
            if name.startswith("P_"):
                ax.plot(*p, "ko", ms=3.5)
    if q is not None:
        _issue_markers(ax, q)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(alpha=0.25)
    ax.set_xlabel("m' / chord")
    ax.set_ylabel(r"$\theta$ / chord")


def _zoom(ax, center, half):
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)


def draw_curvature(ax, q: dict, prof: "bg.Profile2D"):
    """Curvature along both sides, LE point -> TE point, vs normalised arc length."""
    for side, c, lbl in (("ss", C_SS, "suction side"), ("ps", C_PS, "pressure side")):
        s = q[side]
        x = s["s"] / s["s"][-1]
        ax.plot(x, s["kappa"], color=c, lw=1.6, label=lbl)
        # piece boundaries (junctions)
        pc = s["piece"]
        for j in np.flatnonzero(np.diff(pc) != 0):
            ax.axvline(x[j + 1], color=c, ls=":", lw=0.8, alpha=0.7)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.axhline(0.0, color="k", lw=0.6)
    for i in q["issues"]:
        if i["side"] in (0, 1) and i["severity"] in ("ERROR", "WARNING"):
            s = q["ss" if i["side"] == 0 else "ps"]
            k = int(np.argmin(np.hypot(*(s["xy"] - i["where"]).T)))
            x = s["s"][k] / s["s"][-1]
            ax.plot(x, s["kappa"][k], "x" if i["severity"] == "ERROR" else "^",
                    color="red" if i["severity"] == "ERROR" else "darkorange", ms=10, mew=2, mfc="none")
    ax.set_xlabel("arc length from the LE point / side length")
    ax.set_ylabel(r"curvature $\kappa$ [1/chord] (symlog)")
    ax.set_title("Curvature: LE point → TE point (dotted: junctions)", fontsize=10)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, loc="upper center")


def draw_thickness(ax, prof: "bg.Profile2D"):
    t = prof.thickness(200)
    ax.plot(t["xi"], t["measured"], color="k", lw=1.6, label="measured (SS normal → PS)")
    ax.plot(t["xi"], t["law"], color="tab:purple", lw=1.2, ls="--", label="thickness law t(ξ)")
    for x in t["xi_ctrl"]:
        ax.axvline(x, color=C_PS, ls=":", lw=0.9)
    ax.axvline(t["xi_max"], color="tab:purple", ls=":", lw=0.9)
    ax.set_xlabel("ξ (SS arc length)")
    ax.set_ylabel("thickness / chord")
    ax.set_title("Thickness (dotted: PS stations, peak)", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)


def profile_dashboard(params: "bg.ProfileParams", fig=None, title: str | None = None):
    """Profile, LE / TE zooms, curvature of both sides, thickness, quality badge."""
    import matplotlib.pyplot as plt

    prof = bg.Profile2D(params)  # never throws for a bad shape; only for invalid input
    q = bg.check_quality(prof)
    if fig is None:
        fig = plt.figure(figsize=(15, 10.5))
    fig.clf()
    gs = fig.add_gridspec(3, 3, height_ratios=[1.15, 1.0, 1.0], hspace=0.42, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    draw_profile(ax, prof, q)
    _badge(ax, q)
    ax.legend(loc="lower right", fontsize=8, ncol=4)
    ax.set_title(title or f"Profile in the unstaggered (m', θ) plane — connector {str(params.connector).split('.')[-1]}",
                 fontsize=11)
    cp = prof.characteristic_points()
    for col, (name, c, center) in enumerate((("LE", C_LE, cp["P_LE"]), ("TE", C_TE, cp["P_TE"]))):
        a = fig.add_subplot(gs[1, col])
        draw_profile(a, prof, q, lw=1.8)
        r = (params.le if name == "LE" else params.te).radius
        _zoom(a, center + np.array([0.6 * (1 if name == "LE" else -1) * r * 3, 0.0]), max(6 * r, 0.004))
        a.set_title(f"{name} zoom", fontsize=10)
    a = fig.add_subplot(gs[1, 2])
    draw_thickness(a, prof)
    a = fig.add_subplot(gs[2, :2])
    draw_curvature(a, q, prof)
    a = fig.add_subplot(gs[2, 2])
    a.axis("off")
    j = prof.junctions()
    rows = [f"{x['name']:<20} gap {x['gap']:.1e}  angle {np.degrees(x['tangent_angle']):.1e}°  Δκ {x['curvature_jump']:.1e}"
            for x in j]
    a.text(0.0, 1.0, "Junction continuity (G2 = gap, angle, Δκ ≈ 0)\n\n" + "\n".join(rows) +
           f"\n\nmax curvature    {q['max_curvature']:.1f}\nmin thickness    {q['min_thickness']:.4f}"
           f"\nSS / PS inflections {q['ss_inflections']} / {q['ps_inflections']}"
           f"\nself-intersections  {q['self_intersections']}",
           va="top", ha="left", family="monospace", fontsize=8)
    return fig, prof, q


def interactive_profile(params: "bg.ProfileParams"):
    """Knobs (grouped by family) + a live dashboard.  Returns the knobs."""
    import matplotlib.pyplot as plt

    knobs = ProfileKnobs(params)
    if not HAVE_WIDGETS:
        profile_dashboard(params)
        plt.show()
        return knobs
    out = W.Output()

    def redraw(p):
        with out:
            out.clear_output(wait=True)
            fig = plt.figure(figsize=(15, 10.5))
            try:
                profile_dashboard(p, fig)
            except ValueError as e:  # invalid input (e.g. a radius <= 0): say why
                fig.clf()
                fig.text(0.05, 0.5, f"invalid parameters: {e}", color="red", fontsize=14)
            plt.show()

    knobs.observe(redraw)
    display(W.VBox([knobs.widget, out]))
    redraw(knobs.params)
    return knobs


# ─────────────────────────────────────────────────────────────────────────────
# 3-D
# ─────────────────────────────────────────────────────────────────────────────

def _have_plotly() -> bool:
    try:
        import plotly  # noqa: F401
        return True
    except ImportError:
        return False


def mesh_traces(shape, color="lightsteelblue", name="blade", deflection=2e-4, opacity=1.0, scale=1000.0):
    """Plotly Mesh3d of an OCCT shape (coordinates in mm)."""
    import plotly.graph_objects as go

    v, f, _ = shape.tessellate(deflection)
    v = v * scale
    return go.Mesh3d(x=v[:, 0], y=v[:, 1], z=v[:, 2], i=f[:, 0], j=f[:, 1], k=f[:, 2], color=color, name=name,
                     opacity=opacity, flatshading=False, showlegend=True,
                     lighting=dict(ambient=0.35, diffuse=0.8, specular=0.3, roughness=0.6, fresnel=0.2))


def show_shapes(shapes: list[tuple], title="", height=650, deflection=2e-4):
    """[(Shape, colour, name), ...] in one 3-D view (plotly, or matplotlib without it)."""
    if _have_plotly():
        import plotly.graph_objects as go

        fig = go.Figure([mesh_traces(s, c, n, deflection) for s, c, n in shapes if not s.is_null])
        fig.update_layout(title=title, height=height, scene=dict(aspectmode="data", xaxis_title="x (axial) [mm]",
                                                                yaxis_title="y [mm]", zaxis_title="z [mm]"),
                          margin=dict(l=0, r=0, t=40, b=0))
        return fig
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(projection="3d")
    allv = []
    for s, c, _n in shapes:
        if s.is_null:
            continue
        v, f, _ = s.tessellate(deflection * 2)
        v = v * 1000.0
        ax.add_collection3d(Poly3DCollection(v[f], facecolor=c, edgecolor="none", alpha=0.9))
        allv.append(v)
    v = np.vstack(allv)
    lo, hi = v.min(0), v.max(0)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
    ax.set_box_aspect(hi - lo)
    ax.set_title(title)
    return fig


def fillet_panel(blade: "bg.Blade", radius_mm: float = 2.0, out_dir: str = ".", name: str = "blade",
                 n_edge: int = 40, n_body: int = 160):
    """A 3-D view with an 'Add hub fillet' button (rolling-ball fillet in C++) and export buttons."""
    state = {"shapes": bg.build_blade(blade, 0.0, n_edge, n_body)}

    def figure():
        s = state["shapes"]
        if s.fillet_ok:
            return show_shapes([(s.blade_on_hub, "lightsteelblue", "hub + blade + fillet")],
                               title=f"{name}: hub fillet r = {radius.value if HAVE_WIDGETS else radius_mm} mm")
        return show_shapes([(s.blade, "lightsteelblue", "blade"), (s.hub, "0.75", "hub sector")], title=f"{name}: blade on hub")

    if not HAVE_WIDGETS:
        s = state["shapes"]
        print(f"blade solid valid {s.blade_valid}, volume {1e6 * s.blade_volume:.3f} cm^3, skin error {s.skin_error:.1e} m")
        sf = bg.build_blade(blade, radius_mm * 1e-3, n_edge, n_body)
        print(f"fillet r = {radius_mm} mm: {sf.fillet_message}; G1 {sf.fillet_g1_deg:.3f} deg")
        state["shapes"] = sf if sf.fillet_ok else s
        return figure()

    radius = W.BoundedFloatText(value=radius_mm, min=0.05, max=50.0, step=0.1, description="radius [mm]")
    btn = W.Button(description="Add hub fillet", button_style="primary", icon="magic")
    clear = W.Button(description="Remove fillet")
    exp = W.Button(description="Export STEP / STL", button_style="success", icon="save")
    status = W.HTML()
    view = W.Output()

    def redraw():
        with view:
            view.clear_output(wait=True)
            f = figure()
            f.show() if hasattr(f, "show") and _have_plotly() else display(f)

    def add(_b):
        status.value = "<i>building the rolling-ball fillet in C++ ...</i>"
        sf = bg.build_blade(blade, radius.value * 1e-3, n_edge, n_body)
        if sf.fillet_ok:
            state["shapes"] = sf
            status.value = (f"<b style='color:green'>✓ {sf.fillet_message}</b> — G1 {sf.fillet_g1_deg:.3f}°, "
                            f"volume {1e6 * sf.blade_on_hub.volume:.3f} cm³")
        else:
            status.value = f"<b style='color:red'>✗ {sf.fillet_message}</b>"
        redraw()

    def remove(_b):
        state["shapes"] = bg.build_blade(blade, 0.0, n_edge, n_body)
        status.value = "fillet removed"
        redraw()

    def export(_b):
        s = state["shapes"]
        os.makedirs(out_dir, exist_ok=True)
        files = []
        for shp, suffix in ((s.blade, "blade"), (s.blade_on_hub, "hub_blade" + ("_fillet" if s.fillet_ok else ""))):
            for ext in (".step", ".stl"):
                p = os.path.join(out_dir, f"{name}_{suffix}{ext}")
                shp.write(p)
                files.append(p)
        status.value = "<b>written:</b> " + ", ".join(os.path.basename(f) for f in files)

    btn.on_click(add)
    clear.on_click(remove)
    exp.on_click(export)
    display(W.VBox([W.HBox([radius, btn, clear, exp]), status, view]))
    redraw()
    return state


def meridional_view(ax, blade: "bg.Blade", sections: list | None = None):
    d = blade.definition
    try:
        h = d.channel.hub.sample(200)
        s = d.channel.shroud.sample(200)
        ax.plot(*h.T * 1000, "k", lw=1.5, label="hub / shroud")
        ax.plot(*s.T * 1000, "k", lw=1.5)
    except Exception:  # pragma: no cover - explicit stream lines only
        pass
    for sec in sections or blade.sections(20, 60):
        ax.plot(*sec["xr"].T * 1000, color="0.5", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("x (axial) [mm]")
    ax.set_ylabel("r [mm]")
    ax.grid(alpha=0.25)
