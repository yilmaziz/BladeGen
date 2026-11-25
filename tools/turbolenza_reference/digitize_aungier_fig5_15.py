"""
digitize_aungier_fig5_15.py
===========================
Validation data for TURBOlenza's blade-surface boundary-layer analysis: digitizes Aungier (2003) Fig. 5-15,
"Boundary Layer Analysis Result" (Axial-Flow Compressors, ASME Press, Ch. 5, printed p. 113): the shape factor
H = delta*/theta (Eq. 5-129) on the suction and pressure surfaces of the NACA 65(12A10)-10 cascade at alpha = 12 deg,
beta_in = 45 deg, M_in = 0.65, against x/c.

    python tools/turbolenza_reference/digitize_aungier_fig5_15.py <aungier_axial_2003_ch05.pdf> [out_dir]

Needs PyMuPDF and NumPy. The figure is the one raster image of the chapter file's PDF page 37 (1057 x 801 px,
stored as white ink on black). Writes fig5_15.json and the two CSV files into out_dir (default
tests/turbolenza/validation/aungier2003_fig5_15). Neither the book nor the image is copied anywhere.
"""

import collections
import json
import sys
from pathlib import Path

import numpy as np
import pymupdf

PDF = sys.argv[1]
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[2] / "tests" / "turbolenza" / "validation" / "aungier2003_fig5_15"
OUT.mkdir(parents=True, exist_ok=True)

doc = pymupdf.open(PDF)
page = doc[36]                                    # PDF page 37 of the chapter file (printed p. 113)
xref = page.get_images(full=True)[0][0]
pix = pymupdf.Pixmap(doc, xref)
if pix.n != 1:
    pix = pymupdf.Pixmap(pymupdf.csGRAY, pix)
a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
if a.shape != (801, 1057):
    raise SystemExit(f"unexpected image size {a.shape}: not the expected Fig. 5-15 raster")

# connected components (8-neighbour) of the ink inside the plot frame
ink = a > 128
Hh, Ww = ink.shape
lab = np.zeros(ink.shape, int)
inner = np.zeros_like(ink)
inner[20:668, 112:1034] = True
m_ = ink & inner
n_ = 0
for i in range(Hh):
    for j in range(Ww):
        if m_[i, j] and lab[i, j] == 0:
            n_ += 1
            q = collections.deque([(i, j)])
            lab[i, j] = n_
            while q:
                y, x = q.popleft()
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)):
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < Hh and 0 <= xx < Ww and m_[yy, xx] and lab[yy, xx] == 0:
                        lab[yy, xx] = n_
                        q.append((yy, xx))

ink = a > 128
ids, counts = np.unique(lab[lab > 0], return_counts=True)
main = lab == ids[np.argmax(counts)]

# ---- axis calibration: tick centres (pixel) -> values; least squares
xt_px = np.array([143.0, 229, 314, 399, 484, 569, 654, 739, 824, 909, 994])   # 105 + run centres, rows 668-675
xt_val = np.linspace(0.0, 1.0, 11)
yt_px = np.array([110.0, 205, 299, 393, 488, 582])                            # 18 + run centres, cols 106-112
yt_val = np.array([3.4, 3.0, 2.6, 2.2, 1.8, 1.4])
ax, bx = np.polyfit(xt_px, xt_val, 1)
ay, by = np.polyfit(yt_px, yt_val, 1)
X = lambda c: ax * c + bx
Y = lambda r: ay * r + by
cal = {"x_per_px": ax, "x_tick_rms": float(np.sqrt(np.mean((X(xt_px) - xt_val) ** 2))),
       "H_per_px": ay, "H_tick_rms": float(np.sqrt(np.mean((Y(yt_px) - yt_val) ** 2))),
       "frame_check": {"top_row16_H": float(Y(16)), "bottom_row677_H": float(Y(677)),
                       "left_frame_col102_x": float(X(102)), "right_frame_col1037_x": float(X(1037))}}


def runs(v):
    out, s = [], None
    for k, b in enumerate(v):
        if b and s is None:
            s = k
        if not b and s is not None:
            out.append((s, k - 1)); s = None
    if s is not None:
        out.append((s, len(v) - 1))
    return out


def centre_points(mask, cols, pick):
    """(col, row centre) of the top / bottom run of each column."""
    pts = []
    for c in cols:
        r = runs(mask[:, c])
        if not r:
            continue
        s, e = r[0] if pick == "top" else r[-1]
        pts.append((c, 0.5 * (s + e)))
    return np.array(pts, float)


def fit(pts):
    """Total-least-squares line through points: (point, direction)."""
    p0 = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - p0)
    return p0, vt[0]


def intersect(l1, l2):
    (p, d), (q, e) = l1, l2
    A = np.array([d, -e]).T
    t = np.linalg.solve(A, q - p)
    return p + t[0] * d


def line_rms(pts, line):
    p, d = line
    n = np.array([-d[1], d[0]])
    return float(np.sqrt(np.mean(((pts - p) @ n) ** 2)))


# ---- leading-edge segments (windows chosen from the column runs; branches identified in a 3x zoom)
segs = {
    "SS_a": centre_points(main, range(157, 179), "bottom"),   # (0, 2.96) -> local min
    "SS_b": centre_points(main, range(184, 215), "bottom"),   # local min -> junction with the PS line
    "SS_c": centre_points(main, range(222, 242), "bottom"),   # steep drop to the transition dot
    "PS_a": centre_points(main, range(157, 179), "top"),      # (0, 2.87) -> peak
    "PS_b": centre_points(main, list(range(184, 215)) + list(range(221, 251)), "top"),  # peak -> kink (through the junction)
    "PS_c": centre_points(main, range(257, 282), "top"),      # steep drop to the transition dot
}
lines = {k: fit(v) for k, v in segs.items()}
fit_rms = {k: line_rms(segs[k], lines[k]) for k in segs}

# start points: the two runs of the first curve column (142)
r142 = runs(main[:, 142])
ss_start = np.array([142.0, 0.5 * sum(r142[0])])
ps_start = np.array([142.0, 0.5 * sum(r142[-1])])

# transition dots: centroid of the disc core (pixels whose 7x7 neighbourhood is all ink)
def disc_centre(c0, c1, r0, r1):
    core = []
    for r in range(r0, r1):
        for c in range(c0, c1):
            if main[r - 3:r + 4, c - 3:c + 4].all():
                core.append((c, r))
    core = np.array(core, float)
    return core.mean(axis=0), len(core)

ss_dot, n_ss = disc_centre(240, 266, 450, 480)
ps_dot, n_ps = disc_centre(278, 304, 508, 537)

V1 = intersect(lines["SS_a"], lines["SS_b"])
J = intersect(lines["SS_b"], lines["SS_c"])         # SS vertex at the junction
P1 = intersect(lines["PS_a"], lines["PS_b"])
P3 = intersect(lines["PS_b"], lines["PS_c"])         # PS kink
J_on_ps = line_rms(np.array([J]), lines["PS_b"])     # distance of the SS junction vertex from the PS line [px]

# ---- downstream of the transition dots: column scan, SS = top, PS = bottom (+ the dashed SS part)
dash = np.zeros_like(main)
for k, n in zip(ids, counts):
    if n < 200:
        ys, xs = np.nonzero(lab == k)
        if xs.min() >= 735 and xs.max() <= 1000 and ys.min() >= 270 and ys.max() <= 335:
            dash |= lab == k
curve = main | dash


def scan(c0, c1, pick, mask=None):
    mask = curve if mask is None else mask
    pts = []
    for c in range(c0, c1 + 1):
        r = runs(mask[:, c])
        if not r:
            continue
        s, e = r[0] if pick == "top" else r[-1]
        if e - s > 8:            # a dot or a steep line: skip
            continue
        pts.append((c, 0.5 * (s + e)))
    return np.array(pts, float)


def douglas_peucker(p, tol):
    if len(p) < 3:
        return p
    a_, b_ = p[0], p[-1]
    d = b_ - a_
    n = np.array([-d[1], d[0]]) / np.hypot(*d)
    dist = np.abs((p - a_) @ n)
    k = int(np.argmax(dist))
    if dist[k] > tol:
        return np.vstack([douglas_peucker(p[:k + 1], tol)[:-1], douglas_peucker(p[k:], tol)])
    return np.vstack([a_, b_])


dash_c0 = int(np.nonzero(dash.any(axis=0))[0].min())
ss_solid = scan(int(ss_dot[0]) + 12, dash_c0 - 1, "top", main)
jump = np.nonzero(np.abs(np.diff(ss_solid[:, 1])) > 10)[0]   # the solid SS line ends; the PS line is next below
if len(jump):
    ss_solid = ss_solid[:jump[0] + 1]
ss_dash = scan(dash_c0, 995, "top", dash)        # the separated (dashed) part: dash pixels only
ps_tail = scan(int(ps_dot[0]) + 12, 995, "bottom")
TOL = 1.0    # px: above the +-0.5 px staircase of a rasterized near-flat line
ss_v = douglas_peucker(ss_solid, TOL)
ss_d = np.vstack([ss_dash[0], ss_dash[-1]])       # a straight dashed line: its two ends
dash_fit = fit(ss_dash)
dash_rms = line_rms(ss_dash, dash_fit)
ps_v = douglas_peucker(ps_tail, TOL)

ss_poly = np.vstack([ss_start, V1, J, ss_dot, ss_v])
ss_sep = np.vstack([ss_v[-1], ss_d])          # separated part: from the end of the solid line along the dashes
ps_poly = np.vstack([ps_start, P1, P3, ps_dot, ps_v])
sep_col = float(ss_v[-1][0])   # separation: the end of the solid suction-surface line

def snap_ends(poly):
    """The 3 px line caps put the first and last centre point about 1.5 px outside x/c = 0 and 1: move the end points
    along their own end segments onto x/c = 0 and 1 (drops an end vertex left inside that cap)."""
    xs = [X(c) for c, _ in poly]
    pts = [list(q) for q in poly]
    while len(pts) > 2 and X(pts[-2][0]) >= 1.0 - 2 * ax:   # a cap vertex next to the end
        pts.pop(-2)
    for end, nxt, target in ((0, 1, 0.0), (-1, -2, 1.0)):
        (c0, r0), (c1, r1) = pts[end], pts[nxt]
        if abs(X(c0) - target) < 3 * ax:
            t = ((target - bx) / ax - c0) / (c1 - c0)
            pts[end] = [c0 + t * (c1 - c0), r0 + t * (r1 - r0)]
    return np.array(pts)


# the solid SS line ends at the separation point (no TE end); only its start is snapped to x/c = 0.  The dashed
# separated part is one straight line: its far end is put on x/c = 1 along the line fitted through the dashes.
ss_poly_snapped = np.vstack([snap_ends(ss_poly[:2])[:1], ss_poly[1:]])
(p_d, d_d) = dash_fit
t_end = ((1.0 - bx) / ax - p_d[0]) / d_d[0]
ss_sep_snapped = np.vstack([ss_sep[:-1], p_d + t_end * d_d])


def to_xy(poly):
    return [[round(float(X(c)), 4), round(float(Y(r)), 4)] for c, r in poly]

out = {
    "source": "Aungier, R. H., Axial-Flow Compressors (ASME Press, 2003), Ch. 5, Fig. 5-15 'Boundary Layer Analysis "
              "Result', printed p. 113: NACA 65(12A10)-10, alpha = 12 deg, beta_in = 45 deg, M_in = 0.65; "
              "boundary-layer shape factor H = delta*/theta (Eq. 5-129) on the suction and pressure surfaces vs x/c, "
              "edge conditions from the potential-flow solution of Fig. 5-6",
    "method": "digitized from the PDF's embedded raster (1057 x 801 px) by tools/turbolenza_reference digitization: axes "
              "calibrated on the 11 x and 6 y tick marks (least squares); the plotted curves are polylines, so their "
              "vertices are taken: near the LE as intersections of straight-line fits through the line centres, the "
              "transition points as the centres of the transition dots, downstream of transition by Douglas-Peucker "
              f"simplification ({TOL} px) of the column-scanned line centres",
    "uncertainty": f"about 0.5 px: x/c +- {0.5 * ax:.4f}, H +- {0.5 * abs(ay):.4f} (one pixel {ax:.5f} in x/c, {abs(ay):.5f} in H)",
    "calibration": cal,
    "le_fit_rms_px": fit_rms,
    "junction_distance_from_ps_line_px": J_on_ps,
    "transition": {"suction": to_xy([ss_dot])[0], "pressure": to_xy([ps_dot])[0], "disc_core_pixels": [n_ss, n_ps]},
    "separation_start_suction_x_c": round(float(X(sep_col)), 4),
    "suction": to_xy(ss_poly_snapped),
    "suction_separated": to_xy(ss_sep_snapped),
    "suction_separated_line_rms_px": dash_rms,
    "pressure": to_xy(snap_ends(ps_poly)),
    "note": "H on the suction surface after separation_start is the dashed part (separated, H_k limited, Aungier p. 111); "
            "the LE segments before transition are laminar (Gruschwitz), after transition turbulent (Head / Green)",
}
(OUT / "fig5_15.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
def regime(x, x_tr, separated=False):
    """laminar before the transition dot, transition at it, turbulent after; the dashed part separated."""
    if separated:
        return "separated"
    return "laminar" if x < x_tr - 1e-9 else ("transition" if abs(x - x_tr) <= 1e-9 else "turbulent")


x_tr_ss, x_tr_ps = out["transition"]["suction"][0], out["transition"]["pressure"][0]
rows_ss = [(x, h, regime(x, x_tr_ss)) for x, h in out["suction"]]
rows_ss += [(x, h, regime(x, x_tr_ss, separated=True)) for x, h in out["suction_separated"][1:]]
rows_ps = [(x, h, regime(x, x_tr_ps)) for x, h in out["pressure"]]
NL = chr(10)
for name, rows in (("H_suction.csv", rows_ss), ("H_pressure.csv", rows_ps)):
    with open(OUT / name, "w", encoding="utf-8", newline="") as f:
        f.write("x_c,H,regime" + NL)
        for x, h, part in rows:
            f.write(f"{x:.4f},{h:.4f},{part}" + NL)
print(json.dumps({k: out[k] for k in ("calibration", "le_fit_rms_px", "junction_distance_from_ps_line_px", "transition",
                                      "separation_start_suction_x_c", "uncertainty")}, indent=1))
print("suction", len(out["suction"]), out["suction"])
print("suction separated", out["suction_separated"], "rms", dash_rms)
print("pressure", len(out["pressure"]), out["pressure"])
