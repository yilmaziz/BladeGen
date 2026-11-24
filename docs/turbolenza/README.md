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
   shows the pitch, the blade angles κ_LE / κ_TE, the incidence and the throat. **Set P0 from a
   reference Reynolds number** solves Re = ρ1 W1 c / μ(T1) at a reference Mach number for Pt (Tt held,
   Sutherland at the inlet static temperature). Dunavant's case: Re 550,000 at M 0.30, Tt 335 K,
   c 0.127 m gives Pt = 79,043 Pa and Re = 1.034e6 at M1 0.65.
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
   - **Boundary layer (Aungier Sec. 5.6):** off by default. One way, or two way by transpiration;
     the laminar c_f form, the wake model, ω_BL, the tolerance, the BL iteration limit, the
     transition Re_θ and the separation H_k. See "Boundary layer" below.
5. **Run.** Every pass prints one line to the program's terminal and to the log pane. The residual
   plot is live. A viscous run prints "[BL iteration k]" at the start of each inviscid solve.
6. **Post-processing.** After a run these open as separate windows:
   - **Residuals:** the same history, drawn once.
   - **Contours:** pressure or Mach, filled, with iso-lines and the periodic copies.
   - **Blade loading:** Mach, p or Dunavant's S against x/c. A two-way run shows the coupled
     loading, with the inviscid one (BL iteration 0) dashed beside it (a check box).
   - **Boundary layer:** H, θ/c, δ*/c, c_f or Re_θ of both surfaces, with transition (o) and
     separation (x). For Aungier's case (β1 45°, M1 0.65, σ 1, stagger 33°) his Fig. 5-15 is
     overlaid by default; the check box switches it for any case. Below: the coupling history
     (change of m and β_out per BL iteration).

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

## Boundary layer (Aungier 2003 Sec. 5.6)

Not in v1. `src/turbolenza/BoundaryLayer.hpp` has the method and `Viscous.hpp` has the coupling.
Both headers list every equation with its number and mark every choice that is ours.

- **The layer.** It is marched along each blade surface from the first blade station (θ = 0 there)
  to the last one:
  - laminar: Gruschwitz (Eqs. 5-110 to 5-125). The laminar c_f form is a setting: 2 + Λ/6
    (Centrifugal 2000 Eq. 15-20, the profile's wall slope, the default) or 1 + Λ/6 (Axial 2003
    Eq. 5-123 as printed);
  - transition: Re_θ > 250 (Eq. 5-126);
  - turbulent: Head / Green (Eqs. 5-127 to 5-135), H_k limited to 2.4;
  - loss: Eqs. 5-139 to 5-141;
  - μ: Sutherland's law at the edge temperature;
  - x: the arc length of each surface curve from the LE point.
- **One way** is Aungier's analysis: the layer from the inviscid edge state.
- **Two way:** the equivalent inviscid flow (Lighthill 1958) on the same grid. The walls carry the
  mass deficit ρ_e W_e δ* as their stream function, and the exit periodic stations carry the
  wake's. The wake is either constant (the TE sum) or marched with τ_w = 0. The BL iterations
  relax m by ω_BL and restart the inviscid solve from the previous density and exit angle.
- **Viscous off:** the solve is v1's, bit for bit. The 144 v1 comparisons of the tests measure the
  same values as before the boundary layer was added.

Validation, measured on 2026-11-24 (hidden test `[.validation]`, a report and not a pass/fail
test). The case is Dunavant's NACA 65-(12A10)10 on v1's frozen blade: β1 45°, M1 0.65,
Re = 1.034e6, medium level, λ_ρ = λ_K = 0.3, tolerances 1e-5. The reference is Aungier's
Fig. 5-15, digitized:

| run | BL iterations | β_out [°] | ω | SS transition x/c | PS transition x/c | SS separation x/c |
|---|---|---|---|---|---|---|
| Fig. 5-15 (one way) | – | – | – | 0.130 | 0.173 | 0.746 |
| one way, c_f 2 + Λ/6 | 1 | 23.571 | 0.0394 | 0.213 | 0.180 | 0.939 |
| one way, c_f 1 + Λ/6 | 1 | 23.571 | 0.0310 | 0.293 | 0.327 | 0.959 |
| two way, constant wake, ω_BL 0.5 | diverged in k = 2 | (30.97 at k = 1) | – | – | – | – |
| two way, marched wake, ω_BL 0.5 | diverged in k = 2 | (31.88 at k = 1) | – | – | – | – |
| two way, constant wake, ω_BL 0.2 | 36 (converged) | 25.868 | 0.0179 | 0.213 | 0.180 | 0.920 |
| two way, marched wake, ω_BL 0.2 | 35 (converged) | 26.465 | 0.0186 | 0.213 | 0.180 | 0.920 |

H against Fig. 5-15 (one way, 2 + Λ/6), at the figure's vertices downstream of both transitions:

- pressure surface: max |ΔH| 0.115, mean 0.059;
- suction surface: max |ΔH| 0.60, mean 0.25. Our SS layer separates later: at x/c 0.7, H is 2.07
  against about 2.6 in the figure;
- TE: H = 2.511 (SS) and 1.564 (PS), against Aungier's 2.512 and 1.577.

With the default ω_BL = 0.5 the two-way iteration diverges on this case. What was measured:

- The inviscid speed of the last SS station falls by 16 % in its last step, so θ grows 2.7-fold in
  that step and δ* reaches 0.055 c.
- The transpiration velocity this asks for at that wall node, d(ρ_e W_e δ*)/dx / ρ, is about 1.5 times the
  edge speed there (222 against 146 m/s, by the backward difference of m).
- The BL iteration k = 1 moves β_out by 7.4°, and the inviscid solve of k = 2 goes sonic at that
  node (TransonicFlowError, refused).

With ω_BL = 0.2 the iteration contracts by about 0.77 per iteration and converges. The constant
wake is continuous at the TE by construction (mismatch 0). The marched wake's first exit jump
differs from the TE sum by 3.0e-3 mdot'.

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
