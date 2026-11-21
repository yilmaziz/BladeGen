"""BladeGen designer: a desktop window for 2-D blade sections (Tk + matplotlib).

    python -m bladegen.gui            (or tools/gui.cmd on Windows, tools/gui.sh on Linux / macOS)

Three stations (hub, mid, tip), each with its own construction method:

* **NACA 65** - the parametric NACA 65-(Cl0 A10) t/c section (A10 mean line, Kovach-Sandercock
  or Emery variant, C2 monotone-curvature leading edge), the blade geometry of TURBOlenza;
* **Custom** - the BladeGen parametric profile, every knob grouped by family;
* **NACA 65 + Custom** - custom parameters inverse-fitted to the NACA 65 section, then editable.

The geometry is always built; the curvature of both sides (LE point -> TE point, with LE / TE
zooms) is redrawn after every change, and the quality check marks S-shaped / looping edges and
other non-manufacturable defects as errors on the plot.  No flowpath, span or 3-D here: the
stations are 2-D sections (cascade plane), saved / exported for the BladeGen -> ENIGMAgrid ->
TURBOlenza chain (design JSON, coordinates, TURBOlenza CascadeSpec).
"""
from __future__ import annotations

import json
import math
import sys
import queue
import threading
import time
import tkinter as tk
from dataclasses import dataclass, field
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import numpy as np

import matplotlib
import matplotlib.patches  # noqa: E402
import matplotlib.ticker  # noqa: E402  (no backend switch: the canvas is embedded explicitly)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

import bladegen as bg  # noqa: E402

C_SS, C_PS, C_LE, C_TE = "tab:blue", "tab:orange", "tab:green", "tab:red"
PIECE_COLOR = {0: C_LE, 1: C_SS, 2: C_TE, 3: C_PS}
FAMILY_ORDER = ["Leading edge", "Trailing edge", "Suction side", "Thickness distribution", "Pressure side"]
METHODS = ["NACA 65", "Custom", "NACA 65 + Custom"]
PLANNED = ["DCA (planned)", "C4 (planned)"]
VARIANTS = ["kovach_sandercock", "emery_base"]
DESIGN_FORMAT = "bladegen-design-1"
UI_SCALE = 1.0  # set from the Tk DPI scaling when the window opens


# ─────────────────────────────────────────────────────────────────────────────
# Stations
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Station:
    name: str
    span: float
    method: str = "NACA 65"
    cl0: float = 1.2
    t_max_over_c: float = 0.10
    variant: str = "kovach_sandercock"
    custom: dict = field(default_factory=lambda: bg.ProfileParams.naca65_1210().to_dict())
    fitted_to: tuple | None = None        # the NACA (cl0, t, variant) the custom parameters were fitted to
    fit_rms: float | None = None
    stagger_deg: float = 30.0
    solidity: float = 1.25

    def naca_key(self) -> tuple:
        return (round(self.cl0, 6), round(self.t_max_over_c, 6), self.variant)

    def to_json(self) -> dict:
        d = {k: getattr(self, k) for k in ("name", "span", "method", "cl0", "t_max_over_c", "variant", "custom",
                                           "fit_rms", "stagger_deg", "solidity")}
        d["fitted_to"] = list(self.fitted_to) if self.fitted_to else None
        return d

    @classmethod
    def from_json(cls, d: dict) -> "Station":
        s = cls(d["name"], float(d.get("span", 0.0)))
        for k in ("method", "cl0", "t_max_over_c", "variant", "custom", "fit_rms", "stagger_deg", "solidity"):
            if k in d and d[k] is not None:
                setattr(s, k, d[k])
        s.fitted_to = tuple(d["fitted_to"]) if d.get("fitted_to") else None
        if s.method not in METHODS:
            raise ValueError(f"unknown construction method {s.method!r}")
        return s


def default_stations() -> list[Station]:
    """A compressor-like hub / mid / tip set: camber and thickness fall, stagger rises toward the tip."""
    return [Station("Hub", 0.0, cl0=1.5, t_max_over_c=0.10, stagger_deg=25.0, solidity=1.6),
            Station("Mid", 0.5, cl0=1.2, t_max_over_c=0.08, stagger_deg=37.0, solidity=1.3),
            Station("Tip", 1.0, cl0=0.8, t_max_over_c=0.06, stagger_deg=48.0, solidity=1.0)]


# ─────────────────────────────────────────────────────────────────────────────
# Evaluation (worker thread: C++ calls only, no Tk / matplotlib)
# ─────────────────────────────────────────────────────────────────────────────

_naca_cache: dict[tuple, "bg.Naca65Section"] = {}
_naca_lock = threading.Lock()


def naca_section(key: tuple) -> "bg.Naca65Section":
    with _naca_lock:
        sec = _naca_cache.get(key)
    if sec is None:
        sec = bg.Naca65Section(bg.Naca65Params(key[0], key[1], key[2]))
        with _naca_lock:
            if len(_naca_cache) > 32:
                _naca_cache.clear()
            _naca_cache[key] = sec
    return sec


def custom_params(d: dict) -> "bg.ProfileParams":
    return bg.ProfileParams.from_dict(d)


def evaluate(st: Station) -> dict:
    """Geometry, curvature and quality of one station (raises ValueError on invalid input)."""
    out: dict = {"method": st.method, "name": st.name}
    if st.method in ("NACA 65", "NACA 65 + Custom"):
        sec = naca_section(st.naca_key())
        out["naca"] = {"q": bg.check_quality(sec), "meta": sec.metadata, "camber": sec.camber_line(200),
                       "sharp_te": sec.sharp_te}
    if st.method == "NACA 65":
        out["q"] = out["naca"]["q"]
        out["le_r"] = out["naca"]["meta"]["le_radius"]
        out["te_r"] = max(out["naca"]["meta"]["te_radius"], 0.004)
        return out
    if st.method == "NACA 65 + Custom" and st.fitted_to is None:
        out["q"] = out["naca"]["q"]
        out["le_r"] = out["naca"]["meta"]["le_radius"]
        out["te_r"] = max(out["naca"]["meta"]["te_radius"], 0.004)
        out["unfitted"] = True
        return out
    prm = custom_params(st.custom)
    prof = bg.Profile2D(prm)
    out["q"] = bg.check_quality(prof)
    out["control"] = prof.control_polygons()
    out["points"] = {k: v for k, v in prof.characteristic_points().items() if k.startswith("P_")}
    out["le_r"] = prm.le.radius * max(prm.le.stretch_ab, 1.0)
    out["te_r"] = max(prm.te.radius * max(prm.te.stretch_ab, 1.0), 0.004)
    if st.method == "NACA 65 + Custom":
        nq = out["naca"]["q"]
        dss, dps = bg.side_distances(prof, nq["ss"]["xy"][::3], nq["ps"]["xy"][::3])
        d = np.concatenate([dss, dps])
        out["deviation"] = {"rms": float(np.sqrt(np.mean(d * d))), "max": float(np.max(np.abs(d)))}
    return out


class Engine(threading.Thread):
    """One worker; a newer request replaces a pending one (sliders stay responsive)."""

    def __init__(self, results: "queue.Queue"):
        super().__init__(daemon=True)
        self.results = results
        self._job = None
        self._cv = threading.Condition()

    def submit(self, tag, station: Station):
        with self._cv:
            self._job = (tag, Station.from_json(station.to_json()))
            self._cv.notify()

    def run(self):
        while True:
            with self._cv:
                while self._job is None:
                    self._cv.wait()
                tag, st = self._job
                self._job = None
            t0 = time.perf_counter()
            try:
                res = evaluate(st)
                res["seconds"] = time.perf_counter() - t0
                self.results.put(("eval", tag, res))
            except Exception as e:  # invalid input: shown in the window
                self.results.put(("eval", tag, {"error": f"{type(e).__name__}: {e}", "method": st.method}))


# ─────────────────────────────────────────────────────────────────────────────
# Drawing
# ─────────────────────────────────────────────────────────────────────────────

def draw_sides(ax, q: dict, lw=1.8, alpha=1.0, ls="-", color=None, label=True):
    for side in ("ss", "ps"):
        s = q[side]
        xy, pc = s["xy"], s["piece"]
        cuts = np.flatnonzero(np.diff(pc) != 0) + 1
        start = 0
        for end in list(cuts) + [len(pc)]:
            seg = xy[max(start - 1, 0):end]
            ax.plot(seg[:, 0], seg[:, 1], ls, color=color or PIECE_COLOR.get(int(pc[start]), "k"), lw=lw, alpha=alpha)
            start = end
    if label and color is None:
        for c, lbl in ((C_SS, "suction side"), (C_PS, "pressure side"), (C_LE, "LE piece"), (C_TE, "TE piece")):
            ax.plot([], [], color=c, lw=lw, label=lbl)


def issue_markers(ax, q: dict):
    for i in q["issues"]:
        if i["severity"] == "ERROR":
            ax.plot(*i["where"], "x", color="red", ms=13, mew=3, zorder=40)
        elif i["severity"] == "WARNING":
            ax.plot(*i["where"], "^", color="darkorange", ms=9, mfc="none", mew=2, zorder=40)


def badge(ax, q: dict, extra: list[str] | None = None):
    if q["errors"]:
        head, fc, ec = f"✗ {q['errors']} GEOMETRY ERROR(S) — not manufacturable as designed", "#ffd6d6", "red"
    elif q["warnings"]:
        head, fc, ec = f"⚠ {q['warnings']} warning(s) — geometry built", "#fff1d6", "darkorange"
    else:
        head, fc, ec = "✓ geometry OK: C2 / G2 junctions, no edge defects", "#dcf5dc", "green"
    lines = [head] + [f"[{i['severity']}] {i['message'][:105]}" for i in q["issues"] if i["severity"] != "INFO"][:5]
    lines += extra or []
    ax.text(0.01, 0.98, "\n".join(lines), transform=ax.transAxes, va="top", ha="left", fontsize=7.5,
            family="monospace", bbox=dict(boxstyle="round,pad=0.4", fc=fc, ec=ec, lw=1.5), zorder=50)


def edge_window(q: dict, which: str) -> float:
    w = []
    for side in ("ss", "ps"):
        sd = q[side]
        x = sd["s"] / sd["s"][-1]
        j = np.flatnonzero(np.diff(sd["piece"]) != 0)
        if len(j):
            w.append(x[j[0] + 1] if which == "LE" else 1.0 - x[j[-1] + 1])
    return min(max(3.0 * max(w), 0.004), 0.5) if w else 0.05


def draw_curvature(ax, q: dict, window: str | None = None, ref: dict | None = None):
    for side, c, lbl in (("ss", C_SS, "suction side"), ("ps", C_PS, "pressure side")):
        s = q[side]
        x = s["s"] / s["s"][-1]
        ax.plot(x, s["kappa"], color=c, lw=1.5, label=lbl)
        for j in np.flatnonzero(np.diff(s["piece"]) != 0):
            ax.axvline(x[j + 1], color=c, ls=":", lw=0.8, alpha=0.7)
        if ref is not None:
            r = ref[side]
            ax.plot(r["s"] / r["s"][-1], r["kappa"], color="0.45", lw=1.0, ls="--", label="NACA 65" if side == "ss" else None)
    for i in q["issues"]:
        if i["side"] in (0, 1) and i["severity"] in ("ERROR", "WARNING"):
            s = q["ss" if i["side"] == 0 else "ps"]
            k = int(np.argmin(np.hypot(*(s["xy"] - i["where"]).T)))
            ax.plot(s["s"][k] / s["s"][-1], s["kappa"][k], "x" if i["severity"] == "ERROR" else "^",
                    color="red" if i["severity"] == "ERROR" else "darkorange", ms=11, mew=2.2, mfc="none", zorder=30)
    ax.set_yscale("symlog", linthresh=1.0)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.grid(alpha=0.25, which="both")
    if window is None:
        ax.set_title("Curvature κ, LE point → TE point (dotted: junctions; < 0: concave)", fontsize=9)
        ax.set_xlabel("arc length / side length", fontsize=8)
        ax.legend(fontsize=7, loc="upper center", ncol=3)
    else:
        w = edge_window(q, window)
        ax.set_xlim((0.0, w) if window == "LE" else (1.0 - w, 1.0))
        # y range from the visible samples
        ks = []
        for side in ("ss", "ps"):
            s = q[side]
            x = s["s"] / s["s"][-1]
            m = (x <= w) if window == "LE" else (x >= 1.0 - w)
            ks.append(s["kappa"][m])
        k = np.concatenate(ks) if ks else np.array([0.0, 1.0])
        if k.size:
            lo, hi = float(k.min()), float(k.max())
            pad = 0.15 * (hi - lo + 1.0)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(f"κ near the {window}", fontsize=9)
    ax.set_ylabel("κ [1/chord]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))


def zoom(ax, c, half):
    ax.set_xlim(c[0] - half, c[0] + half)
    ax.set_ylim(c[1] - half, c[1] + half)


def rotated(xy: np.ndarray, stagger_deg: float) -> np.ndarray:
    """Cascade frame: chord rotated by the stagger from the axial (x) direction, LE at the origin."""
    a = math.radians(stagger_deg)
    c, s = math.cos(a), math.sin(a)
    return np.column_stack([c * xy[:, 0] + s * xy[:, 1], -s * xy[:, 0] + c * xy[:, 1]])


# ─────────────────────────────────────────────────────────────────────────────
# Widgets
# ─────────────────────────────────────────────────────────────────────────────

class Knob(ttk.Frame):
    """Label, slider and entry for one float; calls on_change(value) on every change."""

    def __init__(self, parent, label, lo, hi, value, on_change, fmt="{:.5g}", width=24):
        super().__init__(parent)
        self.lo, self.hi, self.fmt, self.on_change = lo, hi, fmt, on_change
        self._mute = False
        self.var = tk.DoubleVar(value=value)
        self.text = tk.StringVar(value=fmt.format(value))
        ttk.Label(self, text=label, width=width, anchor="w").grid(row=0, column=0, sticky="w")
        self.scale = ttk.Scale(self, from_=lo, to=hi, variable=self.var, command=self._slid, length=int(150 * UI_SCALE))
        self.scale.grid(row=0, column=1, sticky="ew", padx=2)
        self.entry = ttk.Entry(self, textvariable=self.text, width=9)
        self.entry.grid(row=0, column=2)
        self.entry.bind("<Return>", self._typed)
        self.entry.bind("<FocusOut>", self._typed)
        self.columnconfigure(1, weight=1)

    def _slid(self, _v):
        v = float(self.var.get())
        self.text.set(self.fmt.format(v))
        if not self._mute:
            self.on_change(v)

    def _typed(self, _e=None):
        try:
            v = float(self.text.get())
        except ValueError:
            self.text.set(self.fmt.format(self.var.get()))
            return
        if v < self.lo or v > self.hi:  # widen the slider instead of clamping the value
            self.lo, self.hi = min(self.lo, v), max(self.hi, v)
            self.scale.configure(from_=self.lo, to=self.hi)
        self.set(v)
        self.on_change(v)

    def set(self, v):
        self._mute = True
        try:
            if v < self.lo or v > self.hi:
                self.lo, self.hi = min(self.lo, v), max(self.hi, v)
                self.scale.configure(from_=self.lo, to=self.hi)
            self.var.set(v)
            self.text.set(self.fmt.format(v))
        finally:
            self._mute = False

    def enable(self, on: bool):
        st = "!disabled" if on else "disabled"
        self.scale.state([st])
        self.entry.state([st])


class Scrolled(ttk.Frame):
    def __init__(self, parent, width=440):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, width=int(width * UI_SCALE), highlightthickness=0)
        sb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(win, width=e.width))
        self.canvas.configure(yscrollcommand=sb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.canvas.bind_all("<MouseWheel>", self._wheel)
        self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-2, "units"))
        self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(2, "units"))

    def _wheel(self, e):
        w = e.widget
        while w is not None:  # scroll only when the pointer is over the panel
            if w is self.canvas:
                step = -1 if e.delta > 0 else 1
                self.canvas.yview_scroll(step * (1 if abs(e.delta) < 120 else abs(e.delta) // 120), "units")
                return
            w = getattr(w, "master", None)


# ─────────────────────────────────────────────────────────────────────────────
# Main window
# ─────────────────────────────────────────────────────────────────────────────

class DesignerApp:
    def __init__(self, root: tk.Tk, stations: list[Station] | None = None):
        global UI_SCALE
        self.root = root
        try:
            UI_SCALE = max(1.0, float(root.tk.call("tk", "scaling")) / (96.0 / 72.0))
        except tk.TclError:
            UI_SCALE = 1.0
        root.title(f"BladeGen {bg.__version__} — blade section designer")
        self.stations = stations or default_stations()
        self.cur = 1  # mid
        self.results: queue.Queue = queue.Queue()
        self.engine = Engine(self.results)
        self.engine.start()
        self.last: dict | None = None
        self._tag = 0
        self._shown_tag = -1
        self._after = None
        self._fitting = False
        self.show_cascade = tk.BooleanVar(value=True)
        self.overlay = tk.BooleanVar(value=False)
        self.overlay_data: dict[str, dict] = {}
        self._overlay_keys: dict[str, str] = {}
        self._build_ui()
        self._load_station()
        self.root.after(40, self._poll)

    # ── layout ───────────────────────────────────────────────────────────────
    def _build_ui(self):
        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("Err.TLabel", foreground="#b00000")
        style.configure("Ok.TLabel", foreground="#006400")
        style.configure("Head.TLabel", font=("TkDefaultFont", 10, "bold"))
        pw = ttk.PanedWindow(self.root, orient="horizontal")
        pw.pack(fill="both", expand=True)
        left = Scrolled(pw)
        pw.add(left, weight=0)
        right = ttk.Frame(pw)
        pw.add(right, weight=1)
        L = left.inner

        # station
        box = ttk.LabelFrame(L, text="Station")
        box.pack(fill="x", padx=6, pady=4)
        self.station_var = tk.IntVar(value=self.cur)
        for i, s in enumerate(self.stations):
            ttk.Radiobutton(box, text=f"{s.name} (span {s.span:.1f})", value=i, variable=self.station_var,
                            command=self._switch_station).grid(row=0, column=i, padx=4, pady=2, sticky="w")
        ttk.Checkbutton(box, text="overlay all stations", variable=self.overlay,
                        command=self._toggle_overlay).grid(row=1, column=0, columnspan=2, sticky="w", padx=4)
        ttk.Button(box, text="copy section to the other stations", command=self._copy_to_all).grid(
            row=1, column=2, sticky="e", padx=4)

        # method
        box = ttk.LabelFrame(L, text="Construction method")
        box.pack(fill="x", padx=6, pady=4)
        self.method_var = tk.StringVar(value="NACA 65")
        for i, m in enumerate(METHODS):
            ttk.Radiobutton(box, text=m, value=m, variable=self.method_var, command=self._method_changed).grid(
                row=0, column=i, padx=4, pady=2, sticky="w")
        for i, m in enumerate(PLANNED):
            ttk.Radiobutton(box, text=m, value=m, state="disabled").grid(row=1, column=i, padx=4, sticky="w")

        # NACA 65 family
        self.naca_box = ttk.LabelFrame(L, text="NACA 65 series (A10 mean line) — TURBOlenza geometry")
        self.naca_box.pack(fill="x", padx=6, pady=4)
        self.k_cl0 = Knob(self.naca_box, "design lift Cl0", 0.0, 2.5, 1.2, lambda v: self._naca_changed("cl0", v), "{:.3f}")
        self.k_cl0.pack(fill="x", padx=4)
        self.k_t = Knob(self.naca_box, "thickness t_max / c", 0.02, 0.20, 0.10, lambda v: self._naca_changed("t_max_over_c", v),
                        "{:.4f}")
        self.k_t.pack(fill="x", padx=4)
        row = ttk.Frame(self.naca_box)
        row.pack(fill="x", padx=4, pady=2)
        ttk.Label(row, text="profile variant", width=24).pack(side="left")
        self.variant_var = tk.StringVar(value=VARIANTS[0])
        cb = ttk.Combobox(row, textvariable=self.variant_var, values=VARIANTS, state="readonly", width=20)
        cb.pack(side="left")
        cb.bind("<<ComboboxSelected>>", lambda e: self._naca_changed("variant", self.variant_var.get()))
        self.naca_info = ttk.Label(self.naca_box, text="", justify="left", font=("TkFixedFont", 8))
        self.naca_info.pack(fill="x", padx=4, pady=2)
        self.fit_row = ttk.Frame(self.naca_box)
        self.fit_btn = ttk.Button(self.fit_row, text="Fit custom parameters to this NACA 65  →", command=self._fit)
        self.fit_btn.pack(side="left")
        self.fit_info = ttk.Label(self.fit_row, text="")
        self.fit_info.pack(side="left", padx=6)

        # custom families
        self.custom_box = ttk.LabelFrame(L, text="Custom parameters (BladeGen profile), by family")
        self.custom_box.pack(fill="both", padx=6, pady=4)
        nbk = ttk.Notebook(self.custom_box)
        nbk.pack(fill="both", expand=True)
        self.knobs: dict[str, Knob] = {}
        fams: dict[str, list] = {f: [] for f in FAMILY_ORDER}
        for p in bg.param_table():
            fams.setdefault(p["family"], []).append(p)
        short = {"Leading edge": "LE", "Trailing edge": "TE", "Suction side": "SS", "Thickness distribution": "Thickness",
                 "Pressure side": "PS"}
        for fam, items in fams.items():
            fr = ttk.Frame(nbk)
            nbk.add(fr, text=short.get(fam, fam))
            ttk.Label(fr, text=fam, style="Head.TLabel").pack(anchor="w", padx=4, pady=(3, 1))
            for it in items:
                k = Knob(fr, it["label"], it["lo"], it["hi"], it["lo"], self._custom_cb(it["name"]))
                k.pack(fill="x", padx=4, pady=1)
                self.knobs[it["name"]] = k
            if fam == "Pressure side":
                self.ps_knobs = []
                for i in range(3):
                    k = Knob(fr, f"PS station {i + 1} m'", 0.05, 0.95, 0.2 + 0.3 * i, self._ps_cb(i), "{:.3f}")
                    k.pack(fill="x", padx=4, pady=1)
                    self.ps_knobs.append(k)
        fr = ttk.Frame(nbk)
        nbk.add(fr, text="Connector")
        self.conn_var = tk.StringVar(value="Fair")
        for i, m in enumerate(["Fair", "C2Parametric", "G2Legacy"]):
            ttk.Radiobutton(fr, text=m, value=m, variable=self.conn_var, command=self._conn_changed).pack(anchor="w", padx=6)
        ttk.Button(fr, text="reset to NACA 65-(12)10 fit", command=self._reset_custom).pack(anchor="w", padx=6, pady=6)

        # cascade
        box = ttk.LabelFrame(L, text="Cascade (view and TURBOlenza export)")
        box.pack(fill="x", padx=6, pady=4)
        self.k_stagger = Knob(box, "stagger [deg from axial]", -70.0, 70.0, 30.0, lambda v: self._cascade_changed("stagger_deg", v),
                              "{:.2f}")
        self.k_stagger.pack(fill="x", padx=4)
        self.k_sol = Knob(box, "solidity c / s", 0.3, 3.0, 1.25, lambda v: self._cascade_changed("solidity", v), "{:.3f}")
        self.k_sol.pack(fill="x", padx=4)
        ttk.Checkbutton(box, text="show cascade view", variable=self.show_cascade, command=self._redraw).pack(anchor="w", padx=4)

        # files
        box = ttk.LabelFrame(L, text="Files")
        box.pack(fill="x", padx=6, pady=4)
        for i, (txt, fn) in enumerate((("Open design…", self._open), ("Save design…", self._save),
                                       ("Export coordinates…", self._export_coords),
                                       ("Export TURBOlenza CascadeSpec…", self._export_turbolenza),
                                       ("Save figure…", self._save_figure))):
            ttk.Button(box, text=txt, command=fn).grid(row=i // 2, column=i % 2, sticky="ew", padx=3, pady=2)
        box.columnconfigure(0, weight=1)
        box.columnconfigure(1, weight=1)

        # quality
        box = ttk.LabelFrame(L, text="Geometry check")
        box.pack(fill="x", padx=6, pady=4)
        self.status = ttk.Label(box, text="", justify="left", wraplength=int(410 * UI_SCALE))
        self.status.pack(fill="x", padx=4, pady=2)
        self.issues = tk.Text(box, height=9, width=56, font=("TkFixedFont", 8), wrap="word")
        self.issues.pack(fill="x", padx=4, pady=2)
        self.issues.tag_configure("ERROR", foreground="#b00000")
        self.issues.tag_configure("WARNING", foreground="#a05a00")
        self.issues.tag_configure("INFO", foreground="#404040")

        # figure
        self.fig = Figure(figsize=(12, 9), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        NavigationToolbar2Tk(self.canvas, right).update()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # ── station state <-> widgets ────────────────────────────────────────────
    @property
    def st(self) -> Station:
        return self.stations[self.cur]

    def _load_station(self):
        s = self.st
        self.method_var.set(s.method)
        self.k_cl0.set(s.cl0)
        self.k_t.set(s.t_max_over_c)
        self.variant_var.set(s.variant)
        self.k_stagger.set(s.stagger_deg)
        self.k_sol.set(s.solidity)
        self._load_custom()
        self._update_enabled()
        self._request()

    def _load_custom(self):
        d = self.st.custom
        for name, k in self.knobs.items():
            k.set(float(d[name]))
        for i, k in enumerate(self.ps_knobs):
            k.set(float(d["ps.mprime_ctrl"][i]))
        self.conn_var.set(d.get("connector", "Fair"))

    def _update_enabled(self):
        m = self.st.method
        naca_on = m in ("NACA 65", "NACA 65 + Custom")
        custom_on = m == "Custom" or (m == "NACA 65 + Custom" and self.st.fitted_to is not None)
        for k in (self.k_cl0, self.k_t):
            k.enable(naca_on)
        for k in list(self.knobs.values()) + self.ps_knobs:
            k.enable(custom_on)
        if m == "NACA 65 + Custom":
            self.fit_row.pack(fill="x", padx=4, pady=3)
            if self.st.fitted_to is None:
                self.fit_info.configure(text="custom knobs unlock after the fit (3-5 s)")
            elif self.st.fitted_to != self.st.naca_key():
                self.fit_info.configure(text="NACA changed since the fit: refit to follow it")
            else:
                self.fit_info.configure(text=f"fitted: RMS {1e2 * (self.st.fit_rms or 0):.4f} % c")
        else:
            self.fit_row.pack_forget()
        self.fit_btn.state(["disabled" if self._fitting else "!disabled"])

    def _switch_station(self):
        self.cur = int(self.station_var.get())
        if self.overlay.get():
            self._overlay_requests()
        self._load_station()

    def _method_changed(self):
        self.st.method = self.method_var.get()
        self._update_enabled()
        self._request()

    def _naca_changed(self, key, v):
        setattr(self.st, key, v)
        self._update_enabled()
        self._request()

    def _cascade_changed(self, key, v):
        setattr(self.st, key, v)
        self._redraw()

    def _custom_cb(self, name):
        def cb(v):
            self.st.custom[name] = float(v)
            self._request()
        return cb

    def _ps_cb(self, i):
        def cb(v):
            st = list(self.st.custom["ps.mprime_ctrl"])
            st[i] = float(v)
            if all(a < b for a, b in zip(st, st[1:])):
                self.st.custom["ps.mprime_ctrl"] = st
                self._request()
            else:
                self.status.configure(text="PS stations must increase along the chord", style="Err.TLabel")
        return cb

    def _conn_changed(self):
        self.st.custom["connector"] = self.conn_var.get()
        self._request()

    def _reset_custom(self):
        self.st.custom = bg.ProfileParams.naca65_1210().to_dict()
        self._load_custom()
        self._request()

    def _copy_to_all(self):
        src = self.st.to_json()
        for i, s in enumerate(self.stations):
            if i != self.cur:
                keep = (s.name, s.span)
                self.stations[i] = Station.from_json(src)
                self.stations[i].name, self.stations[i].span = keep
        self.overlay_data.clear()
        self._overlay_keys.clear()
        if self.overlay.get():
            self._overlay_requests()
        self._request()

    # ── evaluation round trip ────────────────────────────────────────────────
    def _request(self):
        if self._after is not None:
            self.root.after_cancel(self._after)
        self._after = self.root.after(60, self._send)

    def _send(self):
        self._after = None
        self._tag += 1
        self.engine.submit(("main", self.cur, self._tag), self.st)
        self.status.configure(text="computing…", style="TLabel")

    def _overlay_requests(self):
        """The other stations, for the overlay (cached: recomputed only when a station changed)."""
        for i, s in enumerate(self.stations):
            if i == self.cur:
                continue
            key = json.dumps(s.to_json(), sort_keys=True)
            if self._overlay_keys.get(s.name) == key:
                continue
            try:
                self.overlay_data[s.name] = evaluate(s)
                self._overlay_keys[s.name] = key
            except Exception:
                self.overlay_data.pop(s.name, None)

    def _toggle_overlay(self):
        if self.overlay.get():
            self._overlay_requests()
        self._redraw()

    def _poll(self):
        try:
            while True:
                kind, tag, res = self.results.get_nowait()
                if kind == "eval" and tag[0] == "main" and tag[1] == self.cur and tag[2] == self._tag:
                    self.last = res
                    self._shown_tag = tag[2]
                    self._redraw()
                elif kind == "fit":
                    self._fit_done(tag, res)
        except queue.Empty:
            pass
        self.root.after(40, self._poll)

    def idle(self) -> bool:
        """True when no request or fit is pending and the latest result is drawn (tests, scripting)."""
        return self._after is None and not self._fitting and self._shown_tag == self._tag

    # ── NACA 65 + Custom ─────────────────────────────────────────────────────
    def _fit(self):
        if self._fitting:
            return
        self._fitting = True
        self._update_enabled()
        idx, key = self.cur, self.st.naca_key()
        self.fit_info.configure(text="fitting custom parameters to the NACA 65 section…")

        def work():
            try:
                sec = naca_section(key)
                r = bg.fit_naca65(sec)
                self.results.put(("fit", (idx, key), {"params": r["params"].to_dict(), "rms": r["rms"], "max": r["max"],
                                                       "seconds": r["seconds"]}))
            except Exception as e:
                self.results.put(("fit", (idx, key), {"error": f"{type(e).__name__}: {e}"}))

        threading.Thread(target=work, daemon=True).start()

    def _fit_done(self, tag, res):
        self._fitting = False
        idx, key = tag
        if "error" in res:
            messagebox.showerror("BladeGen", f"fit failed: {res['error']}")
        else:
            st = self.stations[idx]
            st.custom = res["params"]
            st.fitted_to = key
            st.fit_rms = res["rms"]
            if idx == self.cur:
                self._load_custom()
            self.status.configure(text=f"fit: RMS {1e2 * res['rms']:.4f} % c, max {1e2 * res['max']:.4f} % c "
                                       f"in {res['seconds']:.1f} s", style="Ok.TLabel")
        self._update_enabled()
        self._request()

    # ── drawing ──────────────────────────────────────────────────────────────
    def _redraw(self):
        res = self.last
        fig = self.fig
        fig.clf()
        if res is None:
            self.canvas.draw_idle()
            return
        if "error" in res:
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.5, 0.5, "invalid input — no geometry:\n" + res["error"], ha="center", va="center", color="red",
                    fontsize=12, transform=ax.transAxes)
            self.status.configure(text=res["error"], style="Err.TLabel")
            self.canvas.draw_idle()
            return
        q = res["q"]
        st = self.st
        casc = self.show_cascade.get()
        gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 1.0, 1.0], hspace=0.42, wspace=0.28,
                              left=0.075, right=0.985, top=0.90, bottom=0.06)
        ax = fig.add_subplot(gs[0, 0:2])
        ax_le = fig.add_subplot(gs[0, 2])
        ax_k = fig.add_subplot(gs[1, 0:2])
        ax_te = fig.add_subplot(gs[1, 2])
        ax_kle = fig.add_subplot(gs[2, 0])
        ax_kte = fig.add_subplot(gs[2, 1])
        ax_c = fig.add_subplot(gs[2, 2])
        ref = res.get("naca") if st.method == "NACA 65 + Custom" and not res.get("unfitted") else None

        for a in (ax, ax_le, ax_te):
            if self.overlay.get():
                for name, o in self.overlay_data.items():
                    if "q" in o:
                        draw_sides(a, o["q"], lw=1.0, alpha=0.6, ls="-", color="0.65", label=False)
            if ref is not None:
                draw_sides(a, ref["q"], lw=1.1, ls="--", color="0.35", label=False)
            draw_sides(a, q, label=a is ax)
            if "naca" in res and (st.method == "NACA 65" or res.get("unfitted")):
                cam = res["naca"]["camber"]
                a.plot(cam[:, 0], cam[:, 1], "-.", color="0.5", lw=0.8)
                m = res["naca"]["meta"]
                a.add_patch(matplotlib.patches.Circle(m["le_center"], m["le_radius"], fill=False, ls=":", color="0.4", lw=0.8))
                if m["te_radius"] > 0:
                    a.add_patch(matplotlib.patches.Circle(m["te_center"], m["te_radius"], fill=False, ls=":", color="0.4", lw=0.8))
            if "control" in res:
                cp = res["control"]
                for kk, c in (("ss", C_SS), ("ps", C_PS), ("le", C_LE), ("te", C_TE)):
                    a.plot(cp[kk][:, 0], cp[kk][:, 1], "--o", color=c, lw=0.5, ms=2.2, alpha=0.4)
                for p in res["points"].values():
                    a.plot(*p, "ko", ms=3)
            issue_markers(a, q)
            a.set_aspect("equal", adjustable="datalim" if a is ax else "box")
            a.grid(alpha=0.25)
            a.tick_params(labelsize=7)
        ax.set_xlabel("x / chord", fontsize=8)
        ax.set_ylabel("y / chord", fontsize=8)
        if self.overlay.get():
            ax.plot([], [], color="0.65", label="other stations")
        if ref is not None:
            ax.plot([], [], "--", color="0.35", label="NACA 65 target")
        ax.legend(fontsize=7, loc="lower right", ncol=2)
        extra = []
        if "deviation" in res:
            extra.append(f"deviation from the NACA 65: RMS {1e2 * res['deviation']['rms']:.4f} % c, "
                         f"max {1e2 * res['deviation']['max']:.4f} % c")
        if res.get("unfitted"):
            extra.append("NACA 65 + Custom: press 'Fit custom parameters' to start editing")
        badge(ax, q, extra)
        zoom(ax_le, (res["le_r"] * 0.8, 0.0), 3.0 * res["le_r"] + 0.002)
        ax_le.set_title("leading edge", fontsize=9)
        zoom(ax_te, (1.0 - res["te_r"] * 0.8, 0.0), 3.0 * res["te_r"] + 0.002)
        ax_te.set_title("trailing edge", fontsize=9)

        nref = ref["q"] if ref is not None else None
        draw_curvature(ax_k, q, None, nref)
        draw_curvature(ax_kle, q, "LE", nref)
        draw_curvature(ax_kte, q, "TE", nref)

        if casc:
            loop = np.vstack([q["ss"]["xy"], q["ps"]["xy"][::-1]])
            pitch = 1.0 / max(st.solidity, 1e-6)
            r = rotated(loop, st.stagger_deg)
            for k in (-1, 0, 1):
                ax_c.fill(r[:, 0], r[:, 1] + k * pitch, color="lightsteelblue" if k else "steelblue", ec="k", lw=0.6)
            ax_c.annotate("", xy=(-0.15, 0.0), xytext=(-0.55, 0.0), arrowprops=dict(arrowstyle="->", color="0.3"))
            ax_c.set_aspect("equal", adjustable="datalim")
            ax_c.set_title(f"cascade: stagger {st.stagger_deg:.1f}°, σ = {st.solidity:.2f}", fontsize=9)
            ax_c.tick_params(labelsize=7)
            ax_c.grid(alpha=0.2)
        else:
            ax_c.axis("off")

        name = self._section_name(res)
        title_col = "#b00000" if q["errors"] else "black"
        fig.suptitle(f"{st.name} station — {name}" + ("   ✗ NOT MANUFACTURABLE AS DESIGNED" if q["errors"] else ""),
                     fontsize=12, color=title_col, fontweight="bold" if q["errors"] else "normal")
        self.canvas.draw_idle()
        self._update_panel(res)

    def _section_name(self, res) -> str:
        st = self.st
        if st.method == "NACA 65":
            return f"{res['naca']['meta']['designation']} ({st.variant})"
        if st.method == "Custom":
            return "custom BladeGen profile"
        return f"custom fitted to {res['naca']['meta']['designation']} ({st.variant})"

    def _update_panel(self, res):
        q = res["q"]
        if "naca" in res:
            m = res["naca"]["meta"]
            self.naca_info.configure(text=(
                f"{m['designation']}: t_max {m['max_thickness']:.5f} c at station {m['max_thickness_station']:.4f}\n"
                f"R_LE {m['le_radius']:.5f} c, R_TE {m['te_radius']:.4f} c"
                f"{'  (sharp TE)' if res['naca']['sharp_te'] else ''}\n"
                f"equivalent arc: camber {m['camber_angle_deg']:.2f}°, chi_LE {m['chi_le_deg']:.2f}°, "
                f"chi_TE {m['chi_te_deg']:.2f}°\n"
                f"LE fair (SS/PS): anchor x {m['le_anchor_x'][0]:.4f}/{m['le_anchor_x'][1]:.4f}, "
                f"fairness {m['le_fairness'][0]:.2f}/{m['le_fairness'][1]:.2f} (limit 4)"))
        else:
            self.naca_info.configure(text="")
        if q["errors"]:
            self.status.configure(text=f"✗ {q['errors']} error(s): the shape is built but is not manufacturable "
                                       f"as designed (see the red marks).", style="Err.TLabel")
        else:
            self.status.configure(text=f"✓ {q['summary']} — {1e3 * res.get('seconds', 0):.0f} ms", style="Ok.TLabel")
        self.issues.configure(state="normal")
        self.issues.delete("1.0", "end")
        if not q["issues"]:
            self.issues.insert("end", "no issues\n")
        for i in q["issues"]:
            self.issues.insert("end", f"[{i['severity']}] {i['message']}\n", i["severity"])
        self.issues.configure(state="disabled")

    # ── files ────────────────────────────────────────────────────────────────
    def design_json(self) -> dict:
        return {"format": DESIGN_FORMAT, "bladegen_version": bg.__version__,
                "stations": [s.to_json() for s in self.stations]}

    def _save(self):
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("BladeGen design", "*.json")],
                                            initialfile="bladegen_design.json")
        if path:
            Path(path).write_text(json.dumps(self.design_json(), indent=2))
            self.status.configure(text=f"saved {path}", style="Ok.TLabel")

    def _open(self):
        path = filedialog.askopenfilename(filetypes=[("BladeGen design", "*.json")])
        if not path:
            return
        try:
            d = json.loads(Path(path).read_text())
            if d.get("format") != DESIGN_FORMAT:
                raise ValueError(f"not a BladeGen design file (format {d.get('format')!r})")
            st = [Station.from_json(s) for s in d["stations"]]
            if len(st) != len(self.stations):
                raise ValueError("a design has exactly the hub, mid and tip stations")
        except Exception as e:
            messagebox.showerror("BladeGen", f"cannot open {path}:\n{e}")
            return
        self.stations[:] = st
        self.overlay_data.clear()
        self._overlay_keys.clear()
        self._load_station()

    def _export_coords(self):
        res = self.last
        if not res or "q" not in res:
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")],
                                            initialfile=f"{self.st.name.lower()}_section.csv")
        if not path:
            return
        q = res["q"]
        with open(path, "w", encoding="utf-8") as f:
            f.write(f"# BladeGen {bg.__version__}: {self.st.name} station, {self._section_name(res)}\n")
            f.write("# chord 1, LE (0,0), TE (1,0), suction side +y; side LE -> TE\nside,x,y,kappa\n")
            for side, tag in (("ss", "SS"), ("ps", "PS")):
                s = q[side]
                for (x, y), k in zip(s["xy"], s["kappa"]):
                    f.write(f"{tag},{x:.12g},{y:.12g},{k:.10g}\n")
        self.status.configure(text=f"wrote {path}", style="Ok.TLabel")

    def turbolenza_specs(self) -> list[dict]:
        """TURBOlenza CascadeSpec fields (turbolenza.b2b.geometry.CascadeSpec), one per NACA 65 station."""
        out = []
        for s in self.stations:
            if s.method != "NACA 65":
                continue
            out.append({"station": s.name, "span": s.span, "blade_type": "NACA65", "solidity": s.solidity,
                        "stagger_deg": s.stagger_deg, "cl0": s.cl0, "t_max_over_c": s.t_max_over_c,
                        "meanline_type": "A10", "profile_variant": s.variant, "surface_continuity": "C2"})
        return out

    def _export_turbolenza(self):
        specs = self.turbolenza_specs()
        skipped = [s.name for s in self.stations if s.method != "NACA 65"]
        if not specs:
            messagebox.showinfo("BladeGen", "TURBOlenza's CascadeSpec takes NACA 65 sections; no station uses the NACA 65 "
                                            "method. Export the coordinates of custom sections instead.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json")],
                                            initialfile="turbolenza_cascades.json")
        if not path:
            return
        Path(path).write_text(json.dumps({"format": "turbolenza-cascade-spec", "source": f"BladeGen {bg.__version__}",
                                          "note": "fields of turbolenza.b2b.geometry.CascadeSpec; pass each entry "
                                                  "without 'station' / 'span' to CascadeSpec(**entry)",
                                          "cascades": specs}, indent=2))
        msg = f"wrote {len(specs)} cascade(s) to {path}"
        if skipped:
            msg += f"; skipped (not NACA 65): {', '.join(skipped)}"
        self.status.configure(text=msg, style="Ok.TLabel")

    def _save_figure(self):
        path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"),
                                                                                ("SVG", "*.svg")],
                                            initialfile=f"{self.st.name.lower()}_section.png")
        if path:
            self.fig.savefig(path, dpi=150)


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="python -m bladegen.gui", description=__doc__.split("\n")[0])
    ap.add_argument("design", nargs="?", help="a BladeGen design JSON to open")
    ap.add_argument("--screenshot", metavar="PNG", help="render the window's figure to PNG and exit (no interaction)")
    args = ap.parse_args(argv)
    if sys.platform == "win32":  # sharp text and correct sizes on scaled (high-DPI) displays
        try:
            import ctypes

            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass
    stations = None
    if args.design:
        d = json.loads(Path(args.design).read_text())
        stations = [Station.from_json(s) for s in d["stations"]]
    root = tk.Tk()
    root.geometry("1600x950")
    app = DesignerApp(root, stations)
    if args.screenshot:
        def snap():
            if app.last is None:
                root.after(100, snap)
                return
            app.fig.savefig(args.screenshot, dpi=110)
            root.destroy()
        root.after(200, snap)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
