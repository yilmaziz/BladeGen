# The meridional through-flow analysis (Aungier 2003, Ch. 7)

`src/turbolenza/Throughflow*.hpp/.cpp` implement the meridional through-flow analysis of R. H. Aungier,
*Axial-Flow Compressors* (ASME Press, 2003), Chapter 7 (printed pp. 153–173). It is the
streamline-curvature method on straight quasi-normals with blade rows between them, in the four model
levels of Sec. 7.6. It runs on a BladeGen design's meridional channel and blade row
(`ThroughflowBladeGen.hpp`). It gives the stream-surface geometry r(m) and the stream-tube thickness b(m)
to the blade-to-blade solver, and it is the hub-to-shroud half of the quasi-3D analysis of Chapter 12
(not in this version).

Every equation was read from the page images and is cited by number where it is implemented. A choice
the chapter does not make is marked "ours" in the code. The interfaces follow the quasi-3D data
contract between Chs. 5, 6, 7, 8 and 12 (the "contract").

| File | Content |
|---|---|
| `Throughflow.hpp` | conventions, the case, `BladeRowModel`, `RowLossModel`, `EndwallModel`, the solution structs, `solve_throughflow` |
| `Throughflow.cpp` | the solver (Figs. 7-3, 7-4), the gas, `SpecifiedBladeRow`, `ConstantEndwall` |
| `ThroughflowNumerics.hpp/.cpp` | Eqs. 7-67 to 7-73 and the ODE solution, Eqs. 7-34 to 7-36 |
| `ThroughflowBladeGen.hpp/.cpp` | BladeGen design → quasi-normals and blade row (in `turbolenza_bladegen`) |
| `tests/turbolenza/test_throughflow*.cpp` | the tests below |

## The procedure

1. **Initial stream surfaces.** Eqs. 7-54 to 7-56 for uniform flow with ε = 0 (p. 165). r is linear in y,
   so the positions are a closed form.
2. **On every quasi-normal, in order** (the linking loop of Figs. 7-3 and 7-4):
   - **Link to the upstream quasi-normal (Sec. 7.3).** With no row: Eqs. 7-40 to 7-42. With a row:
     Eqs. 7-43 to 7-49 and the row model, then Eqs. 7-50 to 7-53.
   - **The flow field (Sec. 7.2).** Eq. 7-34 with f1 from Eq. 7-16 (C_θ known) or 7-19 (β' known), and
     f4 = f2 + f3 ρΔA/Δṁ (Eqs. 7-17/7-18 or 7-20/7-21, and 7-31 to 7-33). The ∂W_m/∂m in f1 comes from
     Eq. 7-29 with the limit of Eq. 7-30. The solution is Eqs. 7-35 and 7-36. W_m(0) is found from
     Eq. 7-5. The state comes from Eqs. 7-37 and 7-38 at constant entropy. Choke is tested by Eq. 7-39.
   - **End-wall model (Chapter 8).** Gives K_B for Eq. 7-5.
   - Repeat until the flow field converges.
   - The simplified models reposition their stream surfaces inside this loop, undamped (p. 168).
3. **Full model only.**
   - Reposition the stream surfaces (Eq. 7-54 inverted at the mass fractions, Eqs. 7-55 and 7-56).
   - Damp the move with Novak's factor (Eq. 7-57, M_m ≤ 0.95).
   - Lower B* by Eq. 7-59 when the position error grows.
   - Recompute m (Eq. 7-58), φ (Eq. 7-2) and κ_m (Eq. 7-9) with Eqs. 7-67 and 7-70, then go back to 2.
4. **Choke.** The full model ends on a choke only after one completed repositioning (Fig. 7-3). The
   simplified models end on it at once (Fig. 7-4). A choke is a reported status, not an exception.

The four models (Sec. 7.6, `ThroughflowModel`):

| Model | Stream-surface geometry | Entropy term | Loop |
|---|---|---|---|
| `FullNormalEquilibrium` | φ, κ_m, ε from the stream surfaces | yes | outer loop with Novak damping |
| `Approximate` | φ and κ_m linear in y between the exact end-wall values | yes | marching |
| `SimpleNonIsentropic` | ε = κ_m = 0 | yes | marching |
| `SimpleIsentropic` | ε = κ_m = 0 | T ∂s/∂y dropped | marching |

## The equations as read, and the reading notes

All of Chapter 7 (Eqs. 7-1 to 7-73) was read from the 170-dpi page images. So were the cited Eqs. 2-22,
2-25, 2-26, 2-28, 3-13, 3-21 to 3-30 and 6-83/6-84.

An independent transcription (`q3d/eq_ch07.md`) was compared with these readings equation by equation.
**No equation reading disagrees.** Both find the same misprints and ambiguities, resolved as follows:

1. **Eq. 7-58 (p. 165)** is printed m = ∫ [1 + (∂r/∂z)²] dz, without the square root.
   - Arc length is dm = [1 + (∂r/∂z)²]^½ dz (Eq. 3-5).
   - **The square root is implemented**, and the code says this is a misprint.
   - The test on a 45° cone gives m = √2 Δz on the walls to 2.2e-16. The printed form would give 2 Δz.
   - ∂r/∂z comes from Eqs. 7-67 to 7-69 and the integral from Eqs. 7-72 and 7-73 (ours).
   - z must increase along every stream surface. Otherwise the case is refused.
2. **Eq. 3-29 (p. 48)** lacks 1/r on the θ-derivatives. This does not affect Ch. 7, since the flow is
   axisymmetric.
3. **Eq. 7-53** prints C for C₂ (cosmetic).
4. **Slips in the prose:**
   - p. 159 says "f3 includes a term involving ∂W_m/∂m". It is f1 (Eqs. 7-16 and 7-19); the code
     follows the equations.
   - p. 159 says "Combining Eqs. (3-27) and (7-21)". It means (7-22).
   - p. 168 says "set ε = κ_m = 0 in Eqs. (7-12) through (7-19)". Ours: the simple models also have
     ε = κ_m = 0 in Eq. 7-29. There the ∂W_m/∂m term is multiplied by sin ε = 0 and drops out.
5. **Eqs. 7-32 and 7-33 omit K_B, while Eq. 7-5 has it.** Δṁ in Eq. 7-33 is taken as the tube's share of
   the K_B-free mass function of Eq. 7-54, Δψ ṁ / K_B. ρΔA/Δṁ = 1/W_m holds only then (contract F7).
6. **Eq. 7-57 (Novak 1973)** was checked only for its units and limits. Novak's paper was not available
   (**not verified**).
7. **Eq. 7-70** reads [∂f/∂m]₁ = [∂f/∂m]₂ = (f₂ − f₁)/(m₂ − m₁). Reading: an end point takes the slope
   of its own pair. This is used at the first and last quasi-normals for φ and κ_m (p. 171). The
   simplified models here need no m-derivative: the simple ones set ε = κ_m = 0, and the approximate
   one takes the walls' exact values. So the "all quasi-normals" use of p. 172 does not arise.
8. **Eqs. 7-71 to 7-73.** The integral over a segment whose right end is interior uses Eq. 7-72 about
   that end. The last segment uses Eq. 7-73 about its left end. An integral from 0 to y_s is the sum
   of the segments (p. 172).
9. **Fig. 7-3.** "Flow choked? No" goes back to "Is I_QN > 1?": the same quasi-normal is linked again
   and recomputed. "Yes" goes to "Have stream surfaces been repositioned?" (Yes: End; No:
   reposition). Both readings agree.

## Our choices (the book is silent)

- **Gas.** Thermally and calorically perfect (Eqs. 2-22, 2-25, 2-28), h = c_p T, with s = 0 at
  (288.15 K, 101325 Pa). The viscosity is a constant (1.8e-5 Pa s); it is used only for Re_c and
  handed to the Ch. 6 and Ch. 8 models.
- **Inlet profiles.** Given at y/y_s of the first quasi-normal and interpolated linearly. An absolute
  angle α is handled as "β' known" with ω = 0 (p. 162).
- **Stream tubes.** Equal mass, ψ_j = j/(n_ss − 1). The tube of surface j runs from mid-way to j − 1 to
  mid-way to j + 1, with half tubes at the walls.
- **W_m(0) search.** It starts where the smallest W_m is 0 and brackets upward. It is stopped by
  Eq. 7-39 (or by Ch. 12's Eq. 12-37, `ChokeTest`), then solved by Illinois regula falsi to
  |Δṁ|/ṁ ≤ 1e-12. A target that would need reverse flow is refused. W_m ≤ 0 on a converged quasi-normal
  is refused.
- **Repositioning (Eq. 7-54, "by interpolation from this function").** The three-point series of
  Eqs. 7-72 and 7-73 is integrated inside the segment and the position found by bisection.
- **Novak's Δm.** The smaller meridional spacing to the two neighbouring quasi-normals along the
  surface.
- **Eq. 7-59.** Applied when the position error exceeds 1.05 times the previous pass's
  (`bstar_growth_threshold`). See "Measured" below for why it is not 1.00.
- **Tolerances.** The flow field to max |ΔW_m|/max W_m ≤ 1e-10; positions to max |Δy|/y_s ≤ 1e-7.
- **Blade geometry on a stream surface (Exercise 7.1).** The straight segment between the surface's
  points on the two quasi-normals of the row crosses the LE and TE lines. The surface's row span is the
  mean of the two crossing fractions; the walls are 0 and 1. The section is interpolated linearly in
  span, and the pitch is taken at the mean radius (r₁ + r₂)/2.
- **Wake blockage.** K_W acts only on the row's discharge quasi-normal.
- **Re_c** = ρ₁W₁c/μ.
- **Blade height** = the mean of the LE and TE line lengths.
- **Row-level losses** (clearance, seal: `RowLossModel`).
  - Spread linearly in ψ, zero at the opposite wall, with the mass-weighted mean equal to the row's
    ΔP_t (p. 147).
  - Added to ω̄_j as ΔP_t,j/(P'_t1 − P₁)_j.
  - The discharge wall states handed to the row-level model use the per-surface losses only.
- **First estimate on the first quasi-normal.** ṁ/(ρ_t A), with A from Eq. 7-60, capped at half of
  √(2H). It is an initial iterate only.
- **Stream-tube thickness** (Ch. 7 defines none). Per surface: b = cos ε Δy (the area term of Eq. 7-32)
  and b_mass from Eq. 12-39. Per sheet (j, j + 1) at its y-midpoint: b = Δṁ_sheet/(2π r ρ W_m), with the
  actual sheet mass and the mean of the two surfaces' rρW_m (contract Sec. 2.6), and the geometric
  thickness. Inside a blade row the consumer interpolates linearly in m (p. 287).

## Signs (contract Sec. 0.3)

- The book's frame is used throughout: W_θ = C_θ − ωr, and tan β' = W_θ/W_m signed toward +θ.
- BladeGen's θ is the book's θ, and **ω = −2π·rpm/60**.
- A positive rpm turns BladeGen's default rotor in the compressor sense. With axial inflow,
  β'₁ = atan(−ωr/C_m) > 0, the same sign as its LE metal angle. The default design gives
  β'₁ = 56.5° at mid-span, against κ₁ = 52.1°.
- TURBOlenza's blade-to-blade convention (positive relative inlet angle for rpm > 0) needs no sign
  change.
- The Ch. 6 cascade sign `s_row` is +1 when κ₁ > κ₂ at mid-span. It is carried by the row; the Ch. 6
  module applies it.

## The interfaces (for the Ch. 6 and Ch. 8 implementations)

All of them are declared in `Throughflow.hpp`, with their conventions at its top.

- **`BladeRowModel::evaluate(const BladeRowInlet&) const -> BladeRowDischarge`**
  - Called once per stream surface in every iteration of the discharge quasi-normal.
  - Input: row and surface indices; ψ; span; ω; Z; s_row; IGV flag; profile family; blade height;
    (z, r) at both quasi-normals; pitch and solidity; the interpolated section (c, γ, θ, κ₁, κ₂, t_b/c,
    a/c, the position of t_max).
  - Inlet state: β'₁, α₁, W₁, W_m1, W_θ1, C_θ1, M'₁, M₁, Re_c, ρ₁, P₁, T₁, h₁, a₁, μ₁, H'₁, P'_t1,
    T'_t1, s₁, I₁.
  - Discharge estimate: W_m2, ρ₂.
  - Output: β'₂ (book sign), ω̄ (Eq. 7-49 form), K_W.
- **`RowLossModel::evaluate(const RowLossInput&) const -> RowLossResult`**
  - Called once per iteration, after the per-surface calls.
  - Input: the row data (ω, Z, s_row, shrouded, clearances, ṁ, blade height, the hub and shroud
    sections) and the wall states at both quasi-normals (z, r, m, φ, C_m, C_θ, W_θ, ρ, P, T, a, μ, P_t,
    T_t).
  - Output: ΔP_t of the clearance and of the seal [Pa] with their sides, and ṁ_leak (forwarded to the
    end-wall model).
- **`EndwallModel::begin_pass(pass, n_qn)`** and **`EndwallModel::evaluate(const EndwallStation&) -> EndwallResult`**
  - Called once per quasi-normal in every flow iteration, in order (Fig. 7-3), and repeated for the same
    quasi-normal until its flow converges.
  - Input:
    - the quasi-normal (index, pass, ṁ, y_s, λ, frame ω);
    - both wall edges: z, r, m along the wall, φ, V_me, C_θ, W_θ, V_θ = C_θ − rω_w, ρ, P, T, a, μ, P_t,
      T_t, ω_w;
    - the rows bounding it: ω, Z, shrouded, ṁ_leak, and per wall the section, γ̄, s̄, δ_c, r and W_θ
      (row frame) at inlet and discharge, Δm and Δ(rC_θ).
  - Output: K_B only (Eqs. 8-76, 8-77). No force defect enters Chapter 7.
- `SpecifiedBladeRow` (β'₂, ω̄, K_W profiles in span) and `ConstantEndwall` (K_B constant, default 1)
  are the stand-ins of this version.

## Not in this version (doors kept open for Chapter 12)

- **Quasi-normals inside a blade row** (`QuasiNormal::inside_blade_row`, `BladeRowSpec::qn_outlet`) are
  refused. Eqs. 7-19 to 7-21 are invalid there; Ch. 12's 12-17 to 12-19 apply.
- **Ch. 12's blade-metal blockage** K_B_blade (Eq. 12-7) has its slot inside the integral (1 here).
- **Ch. 12's in-passage damping** (Eq. 12-35) is not implemented.
- **Already in place:** the curvature damping D of Eq. 12-38 (default 0), the Mach cap and B*, and
  Ch. 12's choke test (Eq. 12-37).
- **Annulus sizing** (Sec. 7.7) and a real-gas equation of state.

## Tests and measured values

`tools\turbolenza_build.cmd tl-tf test_turbolenza "[throughflow]"` runs 19 cases. Their numbers are
printed in the test output.

### Numerics

- **Eqs. 7-67 to 7-71 on non-uniform points.** Exact for a quadratic: |f′ − exact| ≤ 1.6e-14 and
  |f″ − 10| ≤ 6.4e-14.
- **Eq. 7-70.** Exact for lines; its error on a quadratic is the Taylor remainder.
- **Eqs. 7-72 and 7-73** integrate quadratics exactly (n = 3, 4, 9).
- **∫ e^{2y} dy.** Errors 2.6e-3, 7.1e-4, 2.0e-4, 3.3e-5 (n = 6, 11, 21, 41, graded points). On equal
  points: 3.3e-3, 7.0e-4, 1.1e-4, 1.5e-5, ratios tending to 8 (third order).
- **ODE 7-34 by 7-35/7-36 against RK4** (200,000 steps). Max relative error 2.3e-3, 3.5e-4, 4.6e-5,
  5.8e-6, 7.4e-7 (n = 6 … 81).

### Closed forms in a straight annulus

The annulus has r = 0.3 to 0.5 m, 7 radial quasi-normals, ṁ = 50 kg/s and air at 101325 Pa, 288.15 K.

- **Uniform axial flow, all four models.** W_m = 83.7105787 m/s. The spread is ≤ 1.3e-12, and ṁ against
  the closed form agrees to 7e-15.
- **Free vortex rC_θ = 32 m²/s.** W_m is uniform to 4.3e-7, which is the linear interpolation of the
  401-point inlet profile.
- **Forced vortex C_θ = 200 r**, against W_m² = W_mh² − 2Ω²(r² − r_h²):

  | n_ss | Eq. 7-33 (Aungier, default) | f3/W_m (Eq. 7-15 as it stands) |
  |---|---|---|
  | 6 | 1.84e-2 | 5.31e-2 |
  | 11 | 1.06e-2 | 7.42e-3 |
  | 21 | 6.44e-3 | 1.30e-3 |
  | 41 | 3.04e-3 | 3.27e-4 |

  Aungier's stream-tube form converges at **first order**, the direct form at second order. The cause
  is the half tubes at the walls, where ρΔA/Δṁ ≠ 1/W_m by O(h): the uniform-flow |b_mass/b − 1| is
  4.3e-2 at the walls against 2.9e-3 inside, at n_ss = 11. This is Exercise 7.4's question.
  `ThroughflowOptions::stream_tube_form = false` selects the direct form.
- **Loss-free rotor, ω = −300 rad/s, β'₂ = 0.7 rad** (simple model; Eqs. 7-19 to 7-21 with f2 ≠ 0),
  against W_m = C r^−k − 2ωt r/((1 + t²)(1 + k)). Errors 9.2e-5, 1.5e-5, 2.1e-6, 2.8e-7
  (n_ss = 6 … 41). |Δs| ≤ 6e-14, and rothalpy is conserved.
  - The full model departs from this by 1.6e-2 at every n_ss. It moves the stream surfaces radially
    behind the rotor, so κ_m and ε ≠ 0 there; the closed form has no curvature.
- **The four models on a straight annulus.** They agree to 1.8e-7 (the position tolerance of the full
  model). With an inlet entropy gradient, the isentropic model departs by 0.71: it drops T ∂s/∂y.
- **Mass on every quasi-normal (Eq. 7-5).** |error| ≤ 1e-12 in every test.
- **A 45° cone (Eq. 7-58).** m = √2 Δz and φ = 45° on the walls, to 2.2e-16 and 5.6e-16.
- **Choke.** ṁ = 500 kg/s is choked on quasi-normal 0 at 121.2603 kg/s, mean M_m = 1.0000. Eq. 12-37
  gives the same.

### Novak damping, on a 21-quasi-normal channel with a wall bump

- **B* = 8.** 36 passes. The position error falls from 9.7e-2 by about 3.4 per pass to 7.5e-8. The
  smallest F is 0.337.
- **B* = 1e4.** Eq. 7-59 lowers B* to 40.1, and the run converges in 54 passes.
- **B* = 1e4 without Eq. 7-59.** The iteration diverges until a quasi-normal has no valid state, and the
  case is refused.

### On BladeGen's default design

The channel is the Bell & Mehta hub 40→43 mm and casing 106→102 mm. The quasi-normals are radial (far
inlet, 1, before LE, after TE, 1, far exit), with ṁ = 3 kg/s.

- **The exact wall curvature** agrees with −Δφ/Δm to 1e-6.
- **Full model.** 734 passes to 1e-7. F = 0.015, since y_s/Δm ≈ 20 on these short extensions.
- **Full model against the linear approximation.** They differ by up to 2.6e-2 in W_m. Before the LE,
  W_m at hub/shroud is 82.34/81.28 m/s (full) against 83.38/83.38 (approximate). The full model's
  κ_m there is −2.6 and +3.4 1/m.
- **Approximate against simple.** On this channel the approximate model equals the simple one, because
  every quasi-normal lies on the cylindrical extensions, where the walls have φ = κ = 0.
- **Quasi-normal count.** With four stations per extension instead of one, the full model needs more
  than 3000 passes: F falls to 0.002.
- **Eq. 7-59 taken literally** ("if that error is increasing") stalls the run. A plateau of a slowly
  decaying mode rises by 1e-5 relative per pass, which drives B* to 1 within ten passes. Hence the 5 %
  growth threshold.
- **Stream-tube thickness.** |b_mass/b − 1| is at most 0.144 on the walls and 0.021 inside. The sheets'
  |b/(K_B b_geometric) − 1| is at most 7.2e-4. For the blade-to-blade solver the sheet values are the
  consistent ones.
- **Rotor at 15000 rpm** (SpecifiedBladeRow: β'₂ = the layers' κ₂, ω̄ = 0.05). Converges in all three
  models tried; the full model takes 745 passes.
  - Simple and approximate models, at mid-span: β'₁ 56.5°, β'₂ 41.0°, W_m 83.4 → 90.2 m/s. Max ΔH
    7.29 kJ/kg, max P_t2/P_t1 1.0875.
  - Full model: β'₁ 56.3°, W_m 83.9 → 90.1 m/s, max ΔH 7.31 kJ/kg, max P_t2/P_t1 1.0877.
  - Euler's ΔH = ωΔ(rC_θ) holds on every surface.
