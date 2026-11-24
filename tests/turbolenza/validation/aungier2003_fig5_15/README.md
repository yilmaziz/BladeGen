# Aungier (2003) Fig. 5-15: boundary-layer shape factor, NACA 65(12A10)-10

Validation data for TURBOlenza's blade-surface boundary-layer analysis (Aungier 2003, Sec. 5.6).

- **Source.** R. H. Aungier, *Axial-Flow Compressors* (ASME Press, 2003), Ch. 5, Fig. 5-15 "Boundary Layer
  Analysis Result", printed p. 113.
- **Case.** The NACA 65(12A10)-10 cascade at α = 12°, β_in = 45°, M_in = 0.65. The same case as Fig. 5-6 and
  TURBOlenza's reference run.
- **Quantity.** The shape factor H = δ*/θ (Eq. 5-129) on the suction and pressure surfaces against the axis
  labelled x/c. See "Which abscissa" below: we read it as Aungier's blade-station coordinate, not as the
  chordwise position of each surface point.
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

## Which abscissa (our reading, 2026-11-24)

The `x_c` column holds the figure's own abscissa. We read it as the **blade-station fraction**
f = (m − m_first) / (m_last − m_first) of Aungier's grid, which has 24 blade stations (23 intervals), not as
the chordwise x/c of each wall node. Aungier does not state this. The evidence:

- The polyline vertices on both surfaces sit at k/23 (within 0.003). Examples: 0.0426/0.0429, 0.0863 and 0.1297/0.1296 near
  the LE; on the pressure surface 0.2191, 0.4789, 0.8705, 0.914, 0.9587 (= 5, 11, 20, 21, 22 /23). Fig. 5-6
  (p. 91) shows the same 1/23 vertex spacing.
- The suction- and pressure-surface vertices **coincide** (within 0.0004). For a 33°-staggered blade, the
  chordwise x/c of the two wall nodes of one station differ by about 0.01 to 0.07 (median 0.049) on
  TURBOlenza's medium grid. Plotting against x/c would separate them by about 8 to 60 px (median 40).
- Both curves start at 0 and end at 1, although the book says the first and last blade nodes lie inside the
  passage (2003 p. 91).
- At the first station H is the Λ = 0 value of the station's Mach number in Fig. 5-6: 0.76 gives 2.965 and
  0.66 gives 2.864, against the figure's 2.955 and 2.880. This is the θ = 0 start.

Read this way, the transition dots are station 3 (0.130) and station 4 (0.174), and the suction-surface
separation starts at station 17 (0.739; the dashed line begins at 0.746, past the 3 px line cap). The H at
separation, 2.66, is H_k = 2.4 converted by Eq. 5-131 at M ≈ 0.62. At the TE, H = 2.51 is H_k = 2.4 at
M ≈ 0.40. Fig. 5-6's "Predicted Separation Point" sits at station 18 (0.78).

`test_bl_validation.cpp` compares on both abscissae. The comparison that counts is the one on f.

The book and the figure image are not in the repository: only these numbers, with the source cited.
