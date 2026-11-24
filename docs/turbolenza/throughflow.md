# The meridional through-flow analysis (Aungier 2003, Ch. 7)

`src/turbolenza/Throughflow*.hpp/.cpp` implement the meridional through-flow analysis of R. H. Aungier,
*Axial-Flow Compressors* (ASME Press, 2003), Chapter 7 (printed pp. 153–173). It is the
streamline-curvature method on straight quasi-normals with blade rows between them, in the four model
levels of Sec. 7.6. It runs on a BladeGen design's meridional channel and blade row
(`ThroughflowBladeGen.hpp`). It gives the stream-surface geometry r(m) and the stream-sheet thickness b(m)
to the blade-to-blade solver, and it is the hub-to-shroud half of the quasi-3D analysis of Chapter 12
(not in this version).

Every equation was read from the page images and is cited by number where it is implemented. A choice
the chapter does not make is marked "ours" in the code.

| File | Content |
|---|---|
| `Throughflow.hpp` | conventions, the case, `BladeRowModel`, `RowLossModel`, `EndwallModel`, the solution structs, `solve_throughflow`, `novak_damping`, `BstarControl` |
| `Throughflow.cpp` | the solver (Figs. 7-3, 7-4), the gas, `SpecifiedBladeRow`, `ConstantEndwall` |
| `ThroughflowNumerics.hpp/.cpp` | Eqs. 7-67 to 7-73 and the ODE solution, Eqs. 7-34 to 7-36 |
| `ThroughflowBladeGen.hpp/.cpp` | BladeGen design → quasi-normals and blade row (in `turbolenza_bladegen`) |
| `tests/turbolenza/test_throughflow*.cpp` | the tests below |

## The procedure

1. **Initial stream surfaces.** Eqs. 7-54 to 7-56 for uniform flow with ε = 0 (p. 165). r is linear in y,
   so the positions are a closed form.
2. **On every quasi-normal, in order** (the linking loop of Figs. 7-3 and 7-4):
   - **Link to the upstream quasi-normal (Sec. 7.3).** With no row: Eqs. 7-40 to 7-42. With a row:
     Eqs. 7-43 to 7-49 with the row models, then Eqs. 7-50 to 7-53.
   - **The flow field (Sec. 7.2).**
     - Eq. 7-34, with f1 from Eq. 7-16 (C_θ known) or Eq. 7-19 (β' known).
     - f4 = f2 + f3 ρΔA/Δṁ, from Eqs. 7-17/7-18 or 7-20/7-21 and Eqs. 7-31 to 7-33.
     - ∂W_m/∂m in f1 from Eq. 7-29, with the limit of Eq. 7-30.
     - The solution by Eqs. 7-35 and 7-36. W_m(0) is found from Eq. 7-5, and choke is tested by Eq. 7-39.
     - The state from Eqs. 7-37 and 7-38 at constant entropy.
   - **End-wall model (Chapter 8).** Gives K_B for Eq. 7-5.
   - Repeat until the flow field converges. The simplified models reposition their stream surfaces
     inside this loop, undamped (p. 168).
3. **Full model only.**
   - Reposition the stream surfaces (Eq. 7-54 inverted at the mass fractions, Eqs. 7-55 and 7-56).
   - Damp the move with Novak's factor (Eq. 7-57, M_m ≤ 0.95). Lower B* by Eq. 7-59 when the position
     error grows.
   - Recompute m (Eq. 7-58), φ (Eq. 7-2) and κ_m (Eq. 7-9) with Eqs. 7-67 and 7-70, then go back to 2.
4. **Choke.** The full model ends on a choke only after one completed repositioning (Fig. 7-3). The
   simplified models end on it at once (Fig. 7-4). A choke is a reported status, not an exception.

The four models (Sec. 7.6, `ThroughflowModel`). In every model, the mass equations (7-5, 7-32/7-33,
7-54) use the **geometric** φ and ε = φ − λ; the model level sets only what the normal-momentum equation
uses:

| Model | Normal momentum: φ, κ_m, ε | Geometric φ (mass equations) | Entropy term | Loop |
|---|---|---|---|---|
| `FullNormalEquilibrium` | from the stream surfaces | the same | yes | outer loop with Novak damping |
| `Approximate` | φ and κ_m linear in y between the end-wall values | the same | yes | marching |
| `SimpleNonIsentropic` | ε = κ_m = 0 | Eq. 7-2 by Eq. 7-70 on the marching pair | yes | marching |
| `SimpleIsentropic` | ε = κ_m = 0 | as above | T ∂s/∂y dropped | marching |

## The equations as read, and the reading notes

All of Chapter 7 (Eqs. 7-1 to 7-73) was read from the 170-dpi page images. So were the cited Eqs. 2-22,
2-25, 2-26, 2-28, 3-13, 3-21 to 3-30 and 6-83/6-84. An independent transcription of the same pages was
compared with these readings equation by equation. **No equation reading disagrees.**

1. **Eq. 7-58 (p. 165)** is printed m = ∫ [1 + (∂r/∂z)²] dz, without the square root.
   - Arc length is dm = [1 + (∂r/∂z)²]^½ dz (Eq. 3-5). **The square root is implemented**, and the code
     calls the printed form a misprint.
   - The test on a 45° cone gives m = √2 Δz on the walls to 2.2e-16. The printed form would give 2 Δz.
   - ∂r/∂z comes from Eqs. 7-67 to 7-69 and the integral from Eqs. 7-72 and 7-73 (ours). z must increase
     along every stream surface; otherwise the case is refused.
2. **Eq. 3-29 (p. 48)** lacks 1/r on the θ-derivatives. This does not affect Ch. 7, since the flow is
   axisymmetric.
3. **Eq. 7-53** prints C for C₂ (cosmetic).
4. **Slips in the prose:**
   - p. 159 says "f3 includes a term involving ∂W_m/∂m". It is f1 (Eqs. 7-16 and 7-19).
   - p. 159 says "Combining Eqs. (3-27) and (7-21)". It means (7-22).
5. **The simplified models (p. 168).** "Set ε = κ_m = 0 in Eqs. (7-12) through (7-19)".
   - That list is the normal-momentum equation only. ε and κ_m are also zeroed in Eq. 7-29 (ours). There
     the ∂W_m/∂m term is multiplied by sin ε = 0 and drops out.
   - The mass equations keep the geometric ε. p. 172 says Eq. 7-70 "will also be required for all
     quasi-normals when the simplified through-flow models … since data at m = m₃ will not be known".
     So the simplified models compute the stream surfaces' φ from Eq. 7-2 by the two-point Eq. 7-70 on
     the marching pair (i − 1, i).
   - On quasi-normal 0 there is no pair. There φ is linear in y between the end-wall angles (ours, as
     the Approximate model does).
6. **K_B in Eqs. 7-32 and 7-33.** Both omit K_B, while Eq. 7-5 has it. Δṁ in Eq. 7-33 is taken as a
   tube's share of the K_B-free mass function of Eq. 7-54, Δψ ṁ / K_B. ρΔA/Δṁ = 1/W_m holds only then.
7. **Eq. 7-57 (Novak 1973)** was checked only for its units and limits. Novak's paper was not available
   (**not verified**).
8. **Eq. 7-70** reads [∂f/∂m]₁ = [∂f/∂m]₂ = (f₂ − f₁)/(m₂ − m₁). An end point takes the slope of its own
   pair. It is used:
   - at the first and last quasi-normals for φ and κ_m in the full model (p. 171);
   - for φ of the marching pair in the simplified models (p. 172).
9. **Eqs. 7-71 to 7-73.** A segment whose right end is interior is integrated by Eq. 7-72 about that
   end. The last segment uses Eq. 7-73 about its left end.
10. **Fig. 7-3.**
    - "Flow choked? No" goes back to "Is I_QN > 1?": the same quasi-normal is linked again and
      recomputed.
    - "Yes" goes to "Have stream surfaces been repositioned?": Yes ends the run, No repositions.

## The stream tubes (Eqs. 7-31 to 7-33)

The tubes of Eq. 7-31 are the sheets between neighbouring surfaces. Sheet k carries
Δṁ = ṁ/((n_ss − 1) K_B). Its ΔA = 2π r K_W cos ε Δy (Eq. 7-32) uses the mean of its two surfaces. The
sheet's ρΔA/Δṁ is 1/W_m at the sheet. It is carried to the surfaces at second order: linear
interpolation between the neighbouring sheet midpoints inside, linear extrapolation at the walls (ours).

An earlier version (c72c100) used node-centred tubes with half tubes at the walls. Those are a
rectangle rule at the wall, and they were first order. The forced vortex then converged at ratios of
about 1.7 per halving.

`stream_tube_form = false` uses f3/W_m of Eq. 7-15 as it stands. Both forms are kept.

## Our choices (the book is silent)

- **Gas.** Thermally and calorically perfect (Eqs. 2-22, 2-25, 2-28), h = c_p T, with s = 0 at
  (288.15 K, 101325 Pa).
- **Viscosity.** Sutherland's law at the local static temperature (the project's choice;
  `sutherland_viscosity` of BoundaryLayer.hpp). It gives μ₁ and Re_c = ρ₁W₁c/μ(T₁) to the blade-row
  model and μ to the end-wall model. A constant is an option.
- **Inlet profiles.** Given at y/y_s of the first quasi-normal and interpolated linearly. An absolute
  angle α is handled as "β' known" with ω = 0 (p. 162).
- **Stream surfaces.** Equal mass, ψ_j = j/(n_ss − 1).
- **W_m(0) search.**
  - It starts where the smallest W_m is 0 and brackets upward.
  - It is stopped by Eq. 7-39 (or by Ch. 12's Eq. 12-37, `ChokeTest`), then solved by Illinois regula
    falsi to |Δṁ|/ṁ ≤ 1e-12.
  - A target that would need reverse flow is refused, and so is W_m ≤ 0 on a converged quasi-normal.
- **Repositioning (Eq. 7-54, "by interpolation from this function").** The three-point series of
  Eqs. 7-72 and 7-73 is integrated inside the segment, and the position found by bisection.
- **Novak's Δm.** The smaller meridional spacing to the two neighbouring quasi-normals.
- **Eq. 7-59** (`BstarControl`). Applied when the position error exceeds 1.05 times the smallest error
  since the last change of B*. A steady growth of less than 5 % per pass therefore also triggers it.
  Compared with the previous pass only, a slowly rising plateau drove B* to 1.
- **Approximate model.** End-wall φ and κ_m as given in `QuasiNormal`, or computed from the wall points
  (Eqs. 7-2 and 7-9 with Eq. 7-67 inside and Eq. 7-70 at the ends; exact for straight walls).
- **Blade geometry on a stream surface (Exercise 7.1).**
  - The straight segment between the surface's points on the two quasi-normals of the row crosses the
    LE and TE lines. The surface's row span is the mean of the two crossing fractions; the walls are 0
    and 1.
  - The section is interpolated linearly in span, and the pitch is taken at the mean radius.
  - BladeGen's throat o (Eq. 6-53) is the smallest distance from the neighbouring blade's LE to the
    section, in the conformal plane at the LE radius (approximate).
- **Wake blockage.** K_W acts only on the row's discharge quasi-normal.
- **Blade height** = the mean of the LE and TE line lengths. N_row is the row's position in the flow
  direction, and rows must be listed in that order.
- **Row-level losses** (`RowLossModel`).
  - Spread linearly in ψ, zero at the opposite wall, with the mass-weighted mean equal to the row's
    ΔP_t (p. 147).
  - Added to ω̄_j as ΔP_t,j/(P'_t1 − P₁)_j.
  - The inlet wall states are handed over with W_θ in the row's frame. The discharge ones use the
    per-surface losses only.
- **End-wall data.**
  - The wall rotation ω_w belongs to the wall segment between two quasi-normals.
  - ω_c (Eq. 8-61) is the case's compressor speed. The BladeGen adapter sets the hub segment across the
    rotor to ω.
- **Stream-sheet thickness** (Ch. 7 defines none).
  - Per sheet at its y-midpoint: b = Δṁ_sheet/(2π r ρ W_m), with the actual sheet mass flow, and its
    geometric thickness.
  - Per surface: the same two thicknesses for a sheet of Δψ centred on the surface. One is cos ε times
    the sheet widths interpolated to the surface (the area term of Eq. 7-32); the other is Eq. 12-39
    with Ch. 7's blockages.
  - Inside a blade row the consumer interpolates linearly in m (p. 287).

## Signs

- **Aungier's frame throughout.** W_θ = C_θ − ωr, and tan β' = W_θ/W_m signed toward +θ.
- **BladeGen's θ is this θ**, and **ω = −2π·rpm/60**.
  - A positive rpm turns BladeGen's default rotor in the compressor sense. With axial inflow,
    β'₁ = atan(−ωr/C_m) > 0, the sign of its LE metal angle: 56.5° against κ₁ = 52.1° at mid-span.
  - The equations of Chs. 3, 7 and 8 are used unchanged with the signed ω. They are invariant under
    θ → −θ.
- **Chapter 6's correlations** hold for positive camber only (p. 128). Each row carries `s_row` (+1 when
  κ₁ > κ₂). The Chapter 6 module multiplies the signed angles and swirls by it on entry, and its
  returned β'₂ on exit.
- **Losses.** Every loss crossing into Eq. 7-49 is ΔP'_t/(P'_t1 − P₁) at the row inlet, on the same
  surface.

## The interfaces (for the Ch. 6 and Ch. 8 implementations)

All of them are declared in `Throughflow.hpp`, with their conventions at its top.

- **`BladeRowModel::evaluate(const BladeRowInlet&) const -> BladeRowDischarge`**
  - Called once per stream surface in every iteration of the discharge quasi-normal.
  - Input, the row and the surface: row, N_row, surface, ψ, span, ω, Z, s_row, IGV flag, profile
    family, blade height, ṁ; (z, r) at both quasi-normals; pitch and solidity; the section (c, γ, θ, κ₁,
    κ₂, t_b/c, a/c, x of t_max, throat o).
  - Input, the inlet state: β'₁, α₁, W₁, W_m1, W_θ1, C_θ1, M'₁, M₁, Re_c, ρ₁, P₁, T₁, h₁, a₁, μ₁, H'₁,
    P'_t1, T'_t1, s₁, I₁.
  - Input, the discharge estimate: W_m2, ρ₂.
  - Output: β'₂ (book sign), ω̄ (Eq. 7-49 form), K_W.
- **`RowLossModel::evaluate(const RowLossInput&) const -> RowLossResult`**
  - Called once per iteration, after the per-surface calls.
  - Input: the row data (N_row, ω, Z, s_row, shrouded, clearances, ṁ, blade height, the hub and shroud
    sections), and the wall states at both quasi-normals with W_θ in the row frame.
  - Output: clearance and seal ΔP_t [Pa] with their sides, and ṁ_leak.
- **`EndwallModel::begin_pass(pass, n_qn)`** and **`EndwallModel::evaluate(const EndwallStation&) -> EndwallResult`**
  - Called once per quasi-normal in every flow iteration, in order (Fig. 7-3, inside the outer loop and
    before the repositioning). Repeated for the same quasi-normal until its flow converges.
  - Input, the quasi-normal: index, pass, ṁ, ω_c, y_s, λ, the frame ω.
  - Input, both wall edges: z, r, m along the wall, φ (geometric), V_me, C_θ, W_θ, V_θ = C_θ − rω_w,
    ρ, P, T, a, μ, P_t, T_t, ω_w of the segment.
  - Input, the rows bounding it: N_row, ω, Z, shrouded, ṁ_leak; per wall the section, γ̄, s̄, δ_c, r and
    W_θ (row frame) at inlet and discharge, Δm, Δ(rC_θ).
  - Output: K_B_endwall only (Eqs. 8-76, 8-77). No force defect and no loss enters Chapter 7.
- **Stand-ins of this version:** `SpecifiedBladeRow` (β'₂, ω̄ and K_W profiles in span) and
  `ConstantEndwall` (K_B constant, default 1).

## Not in this version (doors kept open for Chapter 12)

- **Quasi-normals inside a blade row** (`QuasiNormal::inside_blade_row`, `BladeRowSpec::qn_outlet`) are
  refused. Eqs. 7-19 to 7-21 are invalid there; Ch. 12's 12-17 to 12-19 apply.
- **The blade-metal blockage** K_B_blade (Eq. 12-7) has its slot inside the integral (1 here).
- **Ch. 12's in-passage damping** (Eq. 12-35) is not implemented.
- **Already in place:** the curvature damping D of Eq. 12-38 (default 0), the Mach cap and B*, and
  Ch. 12's choke test (Eq. 12-37).
- **Also not implemented:** annulus sizing (Sec. 7.7) and a real-gas equation of state.

## Tests and measured values

`tools\turbolenza_build.cmd tl-tf test_turbolenza "[throughflow]"` runs 36 cases. Their numbers are
printed in the test output.

### Numerics

- **Eqs. 7-67 to 7-71 on non-uniform points.** Exact for a quadratic: |f′ − exact| ≤ 1.6e-14 and
  |f″ − 10| ≤ 6.4e-14.
- **Eq. 7-70.** Exact for lines; its error on a quadratic is the Taylor remainder.
- **Eqs. 7-72 and 7-73** integrate quadratics exactly (n = 3, 4, 9).
- **∫ e^{2y} dy.**
  - Graded points: errors 2.6e-3, 7.1e-4, 2.0e-4, 3.3e-5 for n = 6, 11, 21, 41.
  - Equal points: errors 3.3e-3, 7.0e-4, 1.1e-4, 1.5e-5, with ratios tending to 8 (third order).
- **ODE 7-34 by 7-35/7-36 against RK4** (200,000 steps). Max relative errors 2.3e-3, 3.5e-4, 4.6e-5,
  5.8e-6, 7.4e-7 for n = 6, 11, 21, 41, 81.

### Closed forms

The annulus is r = 0.3 to 0.5 m with 7 radial quasi-normals, ṁ = 50 kg/s, and air at 101325 Pa,
288.15 K.

- **Uniform axial flow, all four models.** W_m = 83.7105787 m/s, with a spread ≤ 1.3e-12. ṁ against the
  closed form agrees to 7e-15.
- **The sheet thickness, uniform flow.** |b_mass/b − 1| is 4.6e-3 at the walls and 1.4e-3 inside.
- **Free vortex.** W_m is uniform to 4.2e-7.
- **Forced vortex C_θ = 150 r**, against W_m² = W_mh² − 2Ω²(r² − r_h²):

  | n_ss | Eq. 7-33 (sheets, second-order transfer) | f3/W_m |
  |---|---|---|
  | 6 | 2.68e-3 | 8.86e-4 |
  | 11 | 1.15e-3 | 1.75e-4 |
  | 21 | 3.85e-4 | 2.90e-5 |
  | 41 | 1.13e-4 | 4.25e-6 |

  The ratios per halving are 2.3, 3.0, 3.4 (tube form, approaching 4) and 5.1, 6.0, 6.8 (direct form).
  With Ω = 200 the profile at the shroud falls steeply (W_m ≈ 40 m/s). There the tube form's wall
  extrapolation of 1/W_m converges more slowly (1.06e-2 at n_ss = 11).
- **45° cone with swirl, self-convergence against 81 surfaces** (ψ = 0, 0.5, 1):
  - tube form: 4.1e-5, 1.2e-5, 2.5e-6 for n_ss = 11, 21, 41 (ratios 3.6, 4.6);
  - direct form: 5.0e-5, 1.7e-5, 3.8e-6 (ratios 3.0, 4.3).
- **Loss-free rotor, ω = −300 rad/s, β'₂ = 0.7** (simple model), against W_m = C r^−k −
  2ωt r/((1 + t²)(1 + k)):
  - radial quasi-normals: 9.1e-5, 1.5e-5, 2.1e-6, 2.8e-7 (n_ss = 6 … 41);
  - quasi-normals leaned 20°: 3.0e-7 at n_ss = 41 (the cos λ of Eq. 7-20).
  - The full model departs by 1.6e-2: it curves the stream surfaces behind the rotor, and the closed
    form has none.
- **The four models on a straight annulus.** They agree to 2.5e-7. With an entropy gradient the
  isentropic model departs by 0.72.

### The geometric ε in every model

ṁ is recomputed from the node data by the trapezoid rule with the test's own λ and ε, on 41 surfaces:

| Geometry | Full | Approximate | Simple |
|---|---|---|---|
| 45° cone, radial | 8.9e-6 | 3.5e-6 | 5.4e-6 |
| straight annulus, leaned 20° | 8.5e-13 | 8.5e-13 | 8.5e-13 |
| 45° cone, leaned 18.4° | 4.9e-6 | 2.0e-6 | 2.3e-6 |

- **The four models' mean W_m** spread by 6.1e-3 (radial cone) and 2.2e-3 (leaned cone).
- **Before the fix,** the simple models gave W_m 46.84 against 66.95 on the cone (cos 45°), and
  ρW_mA/ṁ = 0.9397 = cos 20° on the leaned annulus.
- **|φ − 45°|** on the cones: ≤ 1.3e-2 (simple, the two-point chords), 7e-3 (full), 6e-16 (approximate).

### Equation by equation

- **Eqs. 7-43 to 7-49, surface by surface.** A rotor with pre-swirl (rC_θ = 20) on a 10° cone, with
  ω̄ = 0.04 and a 500 Pa clearance loss.
  - |I₂ − (H₁ − ω rC_θ1)| = 0.
  - |s₂ − Eq. 7-49| = 4.3e-14 J/(kg K), with the loss spread linearly to 5e-17. r₂ − r₁ is up to 0.036 m.
- **Eq. 7-19 against Eq. 7-16.** The bump channel at QN 0, |κ_m| up to 0.12 1/m: the inlet by α and by
  C_θ agree to 2.0e-6, and the node α = 0.4 to 6e-17.
- **Radial against leaned quasi-normals** (ε up to 0.62 rad), hub W_m:
  - bump with a free vortex: 2.0e-3 (full), 8.7e-3 (approximate);
  - swirling 45° cone (M_θ = 0.46): 4.1e-4 (full), 6.3e-3 (approximate).
- **Eq. 7-57.**
  - By hand: 1/1.375 and 1/1.04875 (the 0.95 cap).
  - In the solver: F = 0.6740856508 against 0.6740856594 recomputed from the converged nodes.
- **Eq. 7-59.** The sequence 1, 0.9, 0.93, 0.96 drops B* from 8 to 4.5. In the solver, B* = 1e4 is first
  lowered to 5000.5.
- **Eq. 7-30.** With the limit 0.97 the limit is active on the 10° cone, and W_m differs by 3.3e-4 from
  the run with 0.1.
- **Choke.**
  - ṁ = 500 chokes at 121.2603 kg/s (by Eq. 7-39 and by Eq. 12-37).
  - With swirl it chokes at 107.4274 kg/s with mean M_m = 1.000000 (trapezoid 1.000045) and max M 1.064.
  - The full model ends after 1, 2, 3 passes for 0, 1, 2 required repositionings.
- **K_W = 0.8 and K_B = 0.85** (a stator, an entropy gradient).
  - The discharge mass from the nodes agrees to 5.8e-6.
  - The tube form against the direct form: 1.1e-6.
- **Eq. 7-70 at the end quasi-normals:** exact (0).
- **Approximate model on a conical diffuser.**
  - The computed wall φ is exact (3e-16).
  - Approximate against full: 6.7e-3.
- **Eq. 8-4:** V_θ = C_θ − rω_w is handed over exactly. μ is from Sutherland's law.

### Novak damping on a 21-quasi-normal bump

- **B* = 8:** 36 passes, error 9.7e-2 → 7.3e-8, min F 0.337.
- **B* = 1e4:** converges in 47 passes. B* is lowered to 40.1.
- **B* = 1e4 without Eq. 7-59:** does not converge (choked after 30 passes).

### BladeGen's default design

Radial quasi-normals (one station per extension), ṁ = 3 kg/s.

- **Full model:** 734 passes. W_m differs from the linear approximation by 2.6e-2.
  - Before the LE, at hub/shroud: 82.34/81.28 m/s (full) against 83.38/83.38 (approximate).
  - The full model's κ_m there is −2.6 and +3.4 1/m.
  - The approximate model equals the simple one, because every station lies on cylindrical
    extensions.
- **Sheet thickness:** |b_mass/b − 1| is 3.9e-2 at the walls and 1.1e-2 inside. The sheets'
  b/(K_B b_geometric) − 1 is 7.2e-4.
- **Rotor at 15000 rpm** (β'₂ = κ₂ of the layers, ω̄ = 0.05):

  | Model | Passes | β'₁ / β'₂ [deg] | W_m1 → W_m2 [m/s] | max ΔH [J/kg] | max P_t2/P_t1 |
  |---|---|---|---|---|---|
  | full | 732 | 56.3 / 41.0 | 83.9 → 90.1 | 7314 | 1.0878 |
  | approximate | 1 | 56.5 / 41.0 | 83.4 → 90.2 | 7292 | 1.0875 |
  | simple | 1 | 56.5 / 41.0 | 83.4 → 90.6 | 7306 | 1.0877 |
