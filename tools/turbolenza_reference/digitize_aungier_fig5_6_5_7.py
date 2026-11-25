"""Digitize the potential-flow curves (solid lines) of Aungier (2003) Figs. 5-6 and 5-7.

Usage:  python digitize_aungier_fig5_6_5_7.py <aungier_axial_2003_ch05.pdf> <output dir>

The chapter PDF is local only (never committed).  The figures are the embedded rasters of PDF pages 15
(Fig. 5-6, M_in 0.65) and 16 (Fig. 5-7, M_in 0.757); the white-on-black curves are tracked column by
column (a line predicted from the last points, the nearest ink run accepted within 3 px; for the Fig. 5-6
pressure side the lowest), in segments
where the dotted experimental curves cross them.  Axis calibration from the frame and the tick marks:
  Fig. 5-6: x/c 0 at column 99.5, 1 at 1025.6; M = (698 - row) / 567.5 (0 .. 1.2)
  Fig. 5-7: x/c 0 at column 100, 1 at 1020;   M = 1.4 - (row - 13) / 574.1 (0.2 .. 1.4)
Fig. 5-6 is fitted piecewise-linearly with knots at k/23 (Aungier's 24 stations, see
tests/turbolenza/validation/aungier2003_fig5_15/README.md); Fig. 5-7 is resampled every 0.02 by the median
of the tracked points within +-0.006, gaps (crossings near the LE and at x/c 0.36-0.38) linear.
"""
import csv
import os
import sys

import numpy as np
import pymupdf
from PIL import Image


def page_image(doc, page_index):
    page = doc[page_index]
    xref = max(page.get_images(full=True), key=lambda im: im[2] * im[3])[0]
    pix = pymupdf.Pixmap(doc, xref)
    return Image.frombytes("L" if pix.n == 1 else "RGB", (pix.width, pix.height), pix.samples).convert("L")


def tracker(a):
    ink = a > 128

    def runs(c, lo, hi):
        col = ink[lo:hi, c]
        out, r = [], 0
        while r < len(col):
            if col[r]:
                s = r
                while r < len(col) and col[r]:
                    r += 1
                seg = a[lo + s:lo + r, c]
                out.append(((np.arange(lo + s, lo + r) * seg).sum() / seg.sum(), r - s))
            else:
                r += 1
        return out

    def track(c0, r0, c1, lo, hi, tol=3.0, hist=8, lowest=False):
        pts, last = [], r0
        for c in range(c0, c1):
            if len(pts) >= hist:
                q = np.array(pts[-hist:])
                pred = np.polyval(np.polyfit(q[:, 0], q[:, 1], 1), c)
            else:
                pred = last
            cand = [x for x in runs(c, lo, hi) if x[1] <= 9 and abs(x[0] - pred) <= tol]
            if cand:
                b = max(cand, key=lambda x: x[0]) if lowest else min(cand, key=lambda x: abs(x[0] - pred))
                pts.append((c, b[0]))
                last = b[0]
        return np.array(pts)

    return track


def hat_fit(x, y, knots):
    def basis(x):
        B = np.zeros((len(x), len(knots)))
        for j in range(len(knots)):
            if j > 0:
                m = (x >= knots[j - 1]) & (x <= knots[j])
                B[m, j] = (x[m] - knots[j - 1]) / (knots[j] - knots[j - 1])
            if j < len(knots) - 1:
                m = (x >= knots[j]) & (x <= knots[j + 1])
                B[m, j] = (knots[j + 1] - x[m]) / (knots[j + 1] - knots[j])
        return B

    w = np.ones_like(x)
    for _ in range(10):  # robust (iteratively reweighted) least squares
        B = basis(x)
        sol = np.linalg.lstsq(B * w[:, None], y * w, rcond=None)[0]
        res = y - B @ sol
        s = 1.4826 * np.median(abs(res)) + 1e-4
        w = 1 / np.maximum(1, abs(res) / (2.5 * s))
    return sol


def main(pdf, out):
    doc = pymupdf.open(pdf)
    os.makedirs(out, exist_ok=True)
    # Fig. 5-6
    a = np.asarray(page_image(doc, 14)).astype(float)
    track = tracker(a)
    X = lambda c: (c - 99.5) / (1025.6 - 99.5)
    Y = lambda r: (698 - r) / 567.5
    R = lambda y: 698 - 567.5 * y
    knots = np.arange(24) / 23
    rows = []
    for y0 in (0.7709, 0.6476):
        p = track(104, R(y0), 1024, 60, 690, hist=6 if y0 > 0.7 else 12, lowest=y0 < 0.7)
        rows.append(hat_fit(X(p[:, 0]), Y(p[:, 1]), knots))
    with open(os.path.join(out, "fig5_6_potential.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "M_suction", "M_pressure"])
        for k, x in enumerate(knots):
            w.writerow([f"{x:.6f}", f"{rows[0][k]:.4f}", f"{rows[1][k]:.4f}"])
    # Fig. 5-7
    a = np.asarray(page_image(doc, 15)).astype(float)
    track = tracker(a)
    X = lambda c: (c - 100.0) / 920.0
    Y = lambda r: 1.4 - (r - 13.0) / 574.1
    R = lambda y: 13.0 + (1.4 - y) * 574.1
    segments = {"SS": [(102, 0.905, 128), (130, 1.03, 410), (420, 1.03, 1019)], "PS": [(102, 0.80, 134), (136, 0.672, 1019)]}
    xs = np.round(np.arange(0, 1.0001, 0.02), 3)
    res = {}
    for name, segs in segments.items():
        p = np.vstack([track(c0, R(y0), c1, 16, 699, hist=6 if c1 - c0 < 60 else 12) for c0, y0, c1 in segs])
        x, y = X(p[:, 0]), Y(p[:, 1])
        v = np.array([np.median(y[abs(x - xi) <= 0.006]) if (abs(x - xi) <= 0.006).any() else np.nan for xi in xs])
        res[name] = v
    with open(os.path.join(out, "fig5_7_potential.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "M_suction", "M_pressure", "interpolated"])
        gap = np.isnan(res["SS"]) | np.isnan(res["PS"])
        for k in ("SS", "PS"):
            m = np.isnan(res[k])
            res[k][m] = np.interp(xs[m], xs[~m], res[k][~m])
        for i, x in enumerate(xs):
            w.writerow([f"{x:.2f}", f"{res['SS'][i]:.4f}", f"{res['PS'][i]:.4f}", int(gap[i])])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
