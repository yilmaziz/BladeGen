# Aungier (2003) Figs. 5-6 and 5-7: potential-flow blade loading, digitized

The solid "Potential Flow Analysis" curves (blade-surface Mach number) of Aungier, Axial-Flow Compressors
(2003), Fig. 5-6 (Dunavant cascade NACA 65-(12A10)10, alpha 12 deg, beta_in 45 deg, M_in 0.65; beta_out 24.4
deg; with the predicted separation point, the filled triangle at x 0.778, M 0.587) and Fig. 5-7 (the same
cascade at M_in 0.757, transonic; beta_out 25.0 deg). Digitized by
`tools/turbolenza_reference/digitize_aungier_fig5_6_5_7.py` from the chapter PDF (local only, never committed),
which reproduces both files exactly.

- `fig5_6_potential.csv`: 24 values per surface at x = k/23 (a robust piecewise-linear fit with knots at
  Aungier's 24 stations; Fig. 5-15 study: the figure's "x/c" axis is read as the blade-station fraction).
- `fig5_7_potential.csv`: every 0.02 in x; `interpolated` = 1 where the dotted experimental curves cross the
  solid ones (near the LE and at 0.36-0.38) and the value is linear between tracked points.
- Accuracy: about +-0.005 in M (pixel size 0.0018 in M; line width ~3 px).

The values are compiled into `apps/designer/turbolenza/TurbolenzaPlots.cpp` (`kFig56*`, `kFig57*`) and must be
kept equal to these files. The TURBOlenza window draws them on the Mach-number loading plot when the
"Aungier (2003)" options are checked.
