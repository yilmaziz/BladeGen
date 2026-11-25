# Dunavant (1955): measured blade loading

`loading.csv`: the measured blade-surface loading of the Dunavant cascade (NACA 65-(12A10)10), as supplied by
the user on 2026-11-25 from Dunavant (1955): per surface, x/c, S = (Pt1 - p)/q1, Cp and the surface Mach number.
The bibliographic details of the report and the test conditions of these points were not checked here; they are
taken to be the Dunavant case of `../aungier2003_fig5_15/README.md`.

The TURBOlenza window overlays these points on the blade-loading plot (Mach number and Dunavant S) when
"Dunavant (1955), measured" is checked; the values are compiled into
`apps/designer/turbolenza/TurbolenzaPlots.cpp` (`kDunavant1955`) and must be kept equal to this file.
