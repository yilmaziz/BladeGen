# Aungier (2003) Fig. 5-15: boundary-layer shape factor, NACA 65(12A10)-10

Validation data for TURBOlenza's blade-surface boundary-layer analysis (Aungier 2003, Sec. 5.6).

- **Source.** R. H. Aungier, *Axial-Flow Compressors* (ASME Press, 2003), Ch. 5, Fig. 5-15 "Boundary Layer
  Analysis Result", printed p. 113.
- **Case.** The NACA 65(12A10)-10 cascade at α = 12°, β_in = 45°, M_in = 0.65. The same case as Fig. 5-6 and
  TURBOlenza's reference run.
- **Quantity.** The shape factor H = δ*/θ (Eq. 5-129) on the suction and pressure surfaces against x/c.
- **How Aungier computed it.** A 2-D boundary-layer analysis whose edge conditions are the inviscid
  potential-flow solution:
  - laminar: Gruschwitz, Eqs. 5-110 to 5-125;
  - transition: Re_θ > 250, Eq. 5-126;
  - turbulent: Head's entrainment method with Green's compressibility corrections, Eqs. 5-127 to 5-135;
  - separation: H_k ≥ 2.4.

  So this figure is a prediction, not a measurement. It verifies an implementation of the same method.

## Files

| File | Content |
|---|---|
| `H_suction.csv`, `H_pressure.csv` | `x_c, H, regime`: the vertices of the plotted polylines. `regime` is `laminar` before the transition dot, `transition` at it, then `turbulent`, and `separated` on the suction surface's dashed part. |
| `fig5_15.json` | The same polylines, plus the transition points, the start of separation, the axis calibration and the fit residuals. |

Key points:

- **Transition** (the dot centres): suction surface x/c = 0.130, H = 1.897; pressure surface x/c = 0.173, H = 1.655.
- **Suction-surface separation:** the solid line ends at x/c = 0.746 (H = 2.66). The dashed separated part runs
  straight to H = 2.51 at the TE.

## How it was digitized

The data comes from the figure's own raster in the PDF (1057 × 801 px), using
`tools/turbolenza_reference/digitize_aungier_fig5_15.py <chapter-5 pdf>`. That command reproduces these files
exactly.

1. **Axes.** Calibrated by least squares on the 11 x and 6 H tick marks. The residuals are 0.0003 in x/c and
   0.0011 in H. The frame lands on H = 3.80 / 1.00 within 0.003.
2. **Polylines.** The plotted curves are polylines, so their vertices are taken:
   - near the LE: the intersections of straight-line fits through the line centres (rms 0.0–0.2 px);
   - the transition points: the centres of the transition dots;
   - downstream of transition: a Douglas–Peucker simplification (1 px) of the column-scanned line centres;
   - the curve ends: snapped to x/c = 0 and 1 along their end segments, which removes the 3 px line caps.
3. **Which curve is which near the LE.** The two curves meet at x/c ≈ 0.086, H ≈ 3.0. A 3× zoom shows which
   continues where:
   - **suction surface:** 2.95 → 2.88 → 3.00, then the straight drop to its transition dot;
   - **pressure surface:** 2.88 → peak 3.21 at x/c 0.043 → 2.80 at x/c 0.130, then its drop.
4. **Uncertainty.** About ±0.5 px: x/c ± 0.0006 and H ± 0.002. This is the reading of the figure. The figure's
   own accuracy is unknown. On the near-flat turbulent plateaus some vertices are pixel steps (±0.004 in H), not
   computed stations.

The book and the figure image are not in the repository: only these numbers, with the source cited.
