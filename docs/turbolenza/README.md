# TURBOlenza in BladeGen

TURBOlenza is a blade-to-blade compressible potential-flow solver (Aungier 2003, Ch. 5) for a planar
cascade of one BladeGen section. It opens as its own window from BladeGen Designer:
**TURBOlenza → TURBOlenza...**. The numerics are a C++ port of the verified TURBOlenza v1 Python
solver (`turbolenza.b2b`, commit 67599a6). The port uses the same equations, discretization, grid
generation, relaxations and stopping rules; see [PORT.md](PORT.md).

## Using it

Work through the tree on the left, top to bottom:

1. **Blade.** Pick a layer (Hub, 25 %, Mid, 75 %, Tip) and press **Take from layer**. This copies the
   layer's profile, stagger and 2-D solidity. Edit the copy with BladeGen's own knobs (NACA 65 /
   Custom / NACA 65 + Custom). The 3-D design does not change.
2. **Cascade.** Chord [m], solidity, stagger, inlet flow angle and Mach number, Pt, Tt, γ, R. The page
   shows the pitch, the blade angles κ_LE / κ_TE, the incidence and the throat.
3. **Mesh.** A standard level (very_coarse … very_fine) or an explicit n_passage × n_eta grid.
   - **Generate mesh** shows the grid and its numbers.
   - The grid lines at the far inlet and exit stations lie along the flow angle (Aungier p. 86).
   - The first and last blade stations lie inside the passage, with one node on each surface
     (Aungier Fig. 5-5). An explicit grid that breaks this is refused, and the message names the
     nearest admissible n_passage.
4. **Solver.** Every input of v1's `scripts/run_full_solver.py`, with its current values as defaults.
   - Cores: 1 to 10.
   - Linear solver: Eigen SparseLU, or Intel MKL PARDISO when BladeGen was built with MKL.
   - **Run** starts the solve and **Stop** ends it after the current pass.
5. **Run.** Every pass prints one line to the program's terminal and to the log pane. The residual
   plot is live.
6. **Post-processing.** After a run these open as separate windows:
   - **Residuals:** the same history, drawn once.
   - **Contours:** pressure or Mach, filled, with iso-lines and the periodic copies.
   - **Blade loading:** Mach, p or Dunavant's S against x/c.

   **Save results...** writes CSV/JSON files and the plots as PNG to a folder you choose.

## Results and verification

On v1's own frozen NACA 65-(12A10)10 blade, the port reproduces v1 to round-off. Measured on
2026-11-24, MSVC 19.51, NumPy 2 / SciPy reference:

- **Grids, all five standard levels:** nodes within 3e-15. The first and last blade stations and
  the wall identity are exact.
- **Each pass, from v1's inputs:** metrics, A..E, the assembled matrix and the flow are bit for bit
  the same. ψ differs by about 1e-14 (Eigen SparseLU against SciPy's SuperLU).
- **The user's case** (M1 0.65, β1 45°, medium level, λ_ρ = λ_K = 0.3, tolerances 1e-5): converged in
  254 passes, as in v1. β_out is 23.5712° in both, with a relative difference of 2.7e-14. It took
  6.3 s on one core with SparseLU.
- **PARDISO on 4 cores:** the same fixed point, to 3e-14.
- **Cores:** with SparseLU, 1 and 4 cores give bit-identical results.

A BladeGen NACA 65 layer is not v1's blade. BladeGen builds its own LE and TE fairs. The body agrees
with v1 to 5e-7 c, the TE region differs by up to 1.9e-4 c, and so the results differ slightly too.

## Building with PARDISO

PARDISO is used when CMake finds Intel MKL:

```
conda install -n bladegen -c conda-forge mkl-devel
```

MKL runs on its TBB threading layer. The per-node loops use the compiler's OpenMP.

## Reference data for the tests

```
python tools\turbolenza_reference\dump_v1.py <path-to-turbolenza-v1>
python tools\turbolenza_reference\dump_v1.py <path-to-turbolenza-v1> --medium
```

These commands write `tests/turbolenza/data/` (git-ignored). Without that folder, the
v1-comparison tests SKIP.
