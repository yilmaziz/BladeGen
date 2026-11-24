# The integration adapters: Chapters 6 and 8 in Chapter 7, Chapter 7 to the blade-to-blade solver

R. H. Aungier, *Axial-Flow Compressors* (ASME Press, 2003) splits a compressor analysis into modules:
the cascade correlations of Ch. 6, the meridional through-flow of Ch. 7, the end-wall boundary layers of
Ch. 8 and the blade-to-blade analysis of Ch. 5, which Ch. 12 (p. 287) feeds from the through-flow in
"the simple one-way approach". TURBOlenza implements each module on its own
(`Correlations.*`, `Throughflow*.*`, `Endwall.*`, `StreamSheet.*` with `Case.*`). The adapters below connect
them. They add no correlation and no equation: they map one module's state into the other's frame and
back. Every mapping the book does not state is marked "ours" in the headers and below.

| File | Content |
|---|---|
| `src/turbolenza/ThroughflowCoupling.hpp/.cpp` | `CorrelationBladeRow` (Ch. 6 as Ch. 7's `BladeRowModel`), `CorrelationRowLoss` (Ch. 6 Secs. 6.9, 6.10 as its `RowLossModel`), `with_correlations`, `ChapterEightEndwall` (Ch. 8 as its `EndwallModel`) |
| `src/turbolenza/ThroughflowSheet.hpp/.cpp` | `throughflow_surface_sheet`, `solver_stream_sheet`, `surface_inlet_state`, `blade_to_blade_on_surface` (Ch. 7 to the blade-to-blade case) |
| `tests/turbolenza/test_throughflow_coupling.cpp` | the Ch. 6 and Ch. 8 adapters, the synthetic stage, BladeGen's default rotor (`[coupling]`) |
| `tests/turbolenza/test_throughflow_sheet.cpp` | the stream sheet, the inlet state, the case settings, one blade-to-blade solve (`[coupling]`; `[.b2bgrid]` hidden) |
| `tests/turbolenza/CouplingFixture.hpp` | the through-flow cases (stage, BladeGen rotor, cone rotor), shared by both test files |

The sheet adapter has a header of its own because it needs the blade-to-blade case (`Case.hpp`), which the
Ch. 6 and Ch. 8 adapters do not.

## Frames (binding)

- **Chapters 7 and 8:** Aungier's frame with a signed ω: W_θ = C_θ − ωr, angles from the meridional direction,
  signed toward +θ. TURBOlenza's and BladeGen's θ is this θ, and ω = −2π·rpm/60. The two modules exchange
  quantities verbatim.
- **Chapter 6:** the positive-camber cascade frame of Fig. 6-1. The row sign s_row maps the book frame onto it
  (cascade angle = s_row × book angle) for the stagger, metal and flow angles, W_θ and ω. `Correlations.hpp`
  applies it (`CascadeGeometry::row_sign = s_row`) and returns β'₂ in the book frame. A guide vane of negative
  camber gets p. 128's own flip on top.
- **Blade-to-blade solver:** its blade has the suction side toward +y. A row with s_row = −1 runs mirrored
  (rpm and the inlet angle sign-changed, the caller's blade the mirrored section); its angles times s_row are
  the book's (ours).

## Chapter 6 → Chapter 7, per stream surface (`CorrelationBladeRow`)

The through-flow calls the model once per stream surface in every iteration of the discharge quasi-normal
(p. 163), with the inlet state at QN1 and its current estimate of the discharge W_m2.

| Ch. 6 input | from the through-flow | note |
|---|---|---|
| r₁, W_m1, W_θ1, ω | `BladeRowInlet` r1, Wm1, Wtheta1 (= C_θ1 − ωr₁, Eq. 7-44), omega | book frame; Ch. 6 multiplies by s_row |
| M'₁ | W₁/a₁ | |
| ρ₁, μ₁ (Re_c) | rho1, mu1 = μ(T₁) by Sutherland's law at the static T₁ | Re_c = ρ₁W₁c/μ₁ (ours; the book names Re_c only) |
| r₂, W_m2 | r2, Wm2 | the iterate (W_m2 = W_m1 on the first call, p. 163) |
| γ, R | R = P₁/(ρ₁T₁), γ = a₁²ρ₁/P₁ | exact for the through-flow's perfect gas (ours) |
| c, θ, t_b/c | the section on the surface | the through-flow's linear interpolation in span of `BladeRowSpec::sections` (a BladeGen layer set or a plain table) |
| a/c | 0.5 for the circular arc and the NACA 65 equivalent arc (Eq. 4-35); the section's value otherwise | ours: BladeGen's largest camber ordinate (e.g. 0.516) is not 0.5 |
| γ (stagger) | MetalAngles (default): γ = κ₁ − χ₁(θ, a/c) in the cascade frame | ours: Ch. 6 measures i and δ from κ₁, κ₂ (6-2, 6-3), so the section's metal angles are kept; `Geometric` passes the section's chord line |
| σ | c/s, s = 2π r_mean/Z, r_mean = (r₁ + r₂)/2 of the quasi-normal points | the through-flow's pitch; ours: the q3d contract takes s at the mean of the LE / TE radii of the cut (0.07244 against 0.07246 m on BladeGen's rotor) |
| o/s | throat / s | ours: o where the section gives it, s at r_mean |
| s̄/h̄ | 2π r̄/Z over the through-flow's blade height, r̄ the mean of the four LE / TE end radii | ours (Eqs. 6-99, 6-100) |
| profile, camberline | `BladeRowSpec::profile_family`: NACA65, C4, DCA, A4K6 | BladeGen's "custom" must be given in the options (refused otherwise) |
| kind | IGV flag, else rotor (ω ≠ 0) or stator | |

Output: β'₂ (book frame), ω̄ = `CascadePerformance::loss` (profile, Mach, shock and off-design: already the
Eq. 7-49 form ΔP'_t/(P'_t1 − P₁) at QN1), K_W = 1 − B_wake (Eq. 7-6).

## Chapter 6 → Chapter 7, the row-level losses (`CorrelationRowLoss`)

Called once per iteration after the per-surface calls, with the end-wall states at QN1 and QN2.

- An unshrouded row: the tip gap at the wall that has a clearance (shroud, or hub for a cantilevered stator;
  gaps at both walls are refused: `RowLossResult` carries one clearance loss).
- A shrouded row: the seal on `seal_wall` (ours: the hub, a stator's inner band) with the fins, fin thickness
  and pitch of the options; a clearance at the other wall is refused.
- Eqs. 6-85 to 6-91 by `endwall_leakage`, with the gap-side wall states (absolute C_θ, book sign; Ch. 6 applies
  s_row), the gap radius the mean of the two wall radii (ours, or the seal radius given), the wall section's chord
  and stagger, N_row, and γ, R from the inlet wall state.
- `leakage_mass_flow` = the seal's ṁ_leak signed + from discharge to inlet (Ch. 8's convention); 0 for a tip gap.
- The through-flow spreads the ΔP_t linearly in mass fraction, zero at the opposite wall (p. 147), and adds
  ΔP_t,j/(P'_t1 − P₁)_j to ω̄_j.

## Chapter 8 → Chapter 7 (`ChapterEightEndwall`)

Fig. 7-3 calls the end-wall analysis once per quasi-normal in every flow iteration, in order, and repeats a
quasi-normal until its flow converges.

- QN 0: `endwall_inlet` (Sec. 8.8; B of `EndwallInlet`, default 0.02, split equally, Eq. 8-72).
- QN i > 0: `endwall_step` from the layer that the last call for QN i − 1 left, with the row data when a row lies
  between i − 1 and i (the step spans it, p. 303).
- Then `endwall_blockage`: K_B = 1 − B (Eqs. 8-76, 8-77) of QN i, the only value returned (p. 195; no force
  defect and no loss enter Ch. 7, p. 175).
- A repeated call for QN i starts from the same layer i − 1, so the same input gives the same K_B. A call out
  of order is refused.
- The edge: the end-wall stream surface at the quasi-normal foot (p. 164), velocities absolute, ω_w of the
  segment ending there (`QuasiNormal::hub_wall_omega`, `shroud_wall_omega`), ω_c = `compressor_omega`
  (Eq. 8-61). Ch. 8 forms V_θ = C_θ − rω_w itself.
- Row walls: s̄ and γ̄ as the through-flow gives them, the clearance at that wall. A shrouded row carries its
  band on `seal_wall` (ours: the hub): clearance 0 there (the seal's leakage replaces it, Sec. 8.6) and the
  row-level model's seal leakage.
- No damping of K_B. p. 195 asks for "some limits and numerical damping" and gives no form; none is applied
  (ours). The end-wall `EndwallLayer::theta12` may be NaN at V_θe = 0; nothing here requires it finite.

## Chapter 7 → the blade-to-blade case (`ThroughflowSheet.hpp`)

Ch. 12 p. 287: the stream-surface geometry from the through-flow at the quasi-normals bounding the row,
"assuming a linear variation between those stations".

- `throughflow_surface_sheet`: surface j of the last pass. z, r, m (Eq. 7-58) at every quasi-normal and the
  sheet thickness b at the surface, in metres:
  - `Carried` (default; the q3d contract's, Sec. 2.6 and F7): b = Δψ ṁ/(2π r ρ W_m), the thickness that carries
    the stream tube's actual mass at the through-flow's ρW_m without blockage, so that the blockage-free
    blade-to-blade solver reproduces Ch. 7's W_m;
  - `Geometric`: the stream-tube spacing normal to the streamline (cos ε × the tubes' y-width, the area term of
    Eq. 7-32), about Carried/(K_B K_W);
  - `Mass`: Eq. 12-39 with Ch. 7's blockages, exactly Carried/(K_B K_W K_B,blade).
  Geometric and Mass depart from the contract where K_B K_W varies along the surface; they are kept for
  comparison. The measured argument for the default: on BladeGen's rotor the blade-to-blade W_m on QN2 differs
  from Ch. 7's by +0.19 % (Carried), +0.39 % (Mass), +0.59 % (Geometric). r(m) and b(m) are linear between neighbouring quasi-normals (p. 287 inside the row; ours outside it)
  and constant beyond the first and last (ours). The LE and TE lie where the surface's segment QN1 → QN2 crosses
  the row's edge lines. `span` is the surface's row span, as the through-flow takes it for the row model.
- `solver_stream_sheet`: M = (m − m_LE)/L, r/L and b on a dense table (StreamSheet's spline follows the
  piecewise-linear r away from the two kinks).
- `surface_inlet_state`: the absolute state at the far inlet (`blade_to_blade_on_surface` takes it at the grid's
  own far inlet, which a standard level rounds to whole cells: it builds the case once to find it, so that
  C_θ = (rC_θ)/r is formed at the radius of the case's velocity triangle), linear in m over the blade-free stretch that ends
  at the row's QN1 (rC_θ, H, s conserved there, Eqs. 7-40 to 7-42); QN 0's state before QN 0, QN1's beyond QN1
  (ours); a far inlet inside an upstream row is refused.
- `blade_to_blade_on_surface`: the caller's blade (BladeGen's `stream_surface_blade` of the layer at the
  surface's span, or a section by hand) with its sheet replaced. The section is not re-placed on the
  through-flow surface (ours; `r_le_ratio`, `sigma_te` and `sigma_te_blade` report the mismatch at the LE and TE). A rotor runs as `RowMotion::Rotating` with
  rpm = −60ω/(2π), a stator as `Stationary` with the absolute Mach number and angle, both × s_row.

## Changes to the reviewed modules

- `Correlations.cpp` (the Ch. 6 review's cosmetic item): with `StallRangeModel::DunavantA4K6`, β_1c and β_1s
  are now those of the model used, γ + α* ∓ 10°, and the "Eq. 6-47 evaluated at β_1c = 20°" note (a limit of
  the unused Eq. 6-47) is no longer written; Eq. 6-47 is not evaluated at all with that model.
- `Correlations.hpp/.cpp`: the Ch. 6 inlet state is renamed `CascadeInletState`. `Flow.hpp` declares an
  `InletState` of its own in the same namespace (17 fields against 7), so a program with both translation units
  had two definitions of one class (an ODR violation: the implicit constructors are emitted under one name, and
  the linker keeps one of them). Every use in Correlations, its tests, the adapters and correlations.md is
  renamed; no other duplicate type name exists in `BladeGen::Turbolenza` (a scan of every header; the forward
  declarations of GridData.hpp are the same classes).
- `Throughflow.cpp`: the section handed to the row model now carries the throat o, interpolated in span like the
  other fields. Before, `BladeRowInlet::geometry.throat` was always 0, so no Ch. 6 model could have formed o/s
  (Eq. 6-53).

## Tests and measured values

`tools\turbolenza_build.cmd <dir> test_turbolenza "[coupling]"` (16 cases, about 9 s in Release).

**The Ch. 6 row adapter, by hand** (a rotor mid surface: r 0.30 → 0.302 m, W_m 101 → 97 m/s, ω = −523.6 rad/s,
κ 54° → 34°):
- every Ch. 6 input against the value written out in the test (Eq. 7-44's W_θ1, M'₁, μ(T₁), o/s, s̄/h̄,
  γ = 44° from the metal angles); γ and R recovered exactly from a state of γ = 1.3, R = 300;
- Re_c equal to the through-flow's ρ₁W₁c/μ(T₁) to 1e-14; K₂ of Eqs. 6-99, 6-100 with K_Re at that Re_c;
- β'₂ = κ₂ + δ (Eq. 6-3) and i = β₁ − κ₁ (Eq. 6-2) with the section's own κ (i 0.052°, i* −0.868°, δ 6.180°,
  δ* 5.697°, β'₂ 40.180°, ω̄ 0.03113, K_W 1, M'₁ 0.513, Re_c 6.9e5);
- the same numbers as `evaluate` called with inputs built in the test (≤ 1e-13); the mirrored surface gives
  −β'₂ and the same ω̄ (≤ 1e-15); a guide vane's β'₂ does not move with β'₁, also through p. 128's flip;
- `Geometric` stagger moves Ch. 4's κ₂ by the mismatch; a parabolic arc (a/c 0.4) keeps both metal angles;
- refusals: Z = 0, a zero end radius, an unknown family, a non-physical state, h = 0, θ ≠ κ₁ − κ₂, a NaN
  angle; Ch. 6's own: the wrong row sign, no throat.

**The wake blockage:** with W_m 101 → 55 m/s on the same surface D_eq = 2.922 and K_W = (2/D_eq)^0.9 = 0.711 to
1e-14 (Eqs. 6-84, 7-6). A single stator (22 vanes, α₁ = −40°) with the row model's discharge estimate frozen
(W_m2 = W_m1; the unfrozen loop does not settle, Findings 1) converges with K_W down to 0.886 on its discharge
quasi-normal: every node's K_W is the adapter's (exactly), (2/D_eq)^0.9 or 1 (Eqs. 6-83, 6-84), and 1 on the
other quasi-normals; Eq. 7-5 from the nodes with K_W carries the mass flow to 1.5e-5, without K_W it is 10.4 % off.

**The row-level losses:** the rotor tip gap against Eqs. 6-85 to 6-89 written out, with a rising casing (r 0.350
→ 0.354 m, the gap radius their mean) and a shroud section of its own axial chord: τ 3.269 N m, ΔP 11609 Pa,
U_c 99.56 m/s, ṁ_c 0.0941 kg/s, ΔP_t 49.68 Pa (to 1e-13); the same for the mirrored row; a hub gap and a hub seal
against `endwall_leakage` (1e-13); the seal's sign (+ for P₂ > P₁, − for P₂ < P₁), a negative leakage handed to
Ch. 8 with its sign, a given seal radius; the refusals.

**The end-wall adapter:**
- an annulus (r = 0.25 to 0.35 m, 7 quasi-normals, 22 kg/s, no rows): K_B(QN 0) = 0.98 exactly (Eq. 8-72 with
  B = 0.02); every station's K_B equal to 1 − Σ2π rρ_eV_meδ₁*/ṁ of its own walls (Eq. 8-76, 1e-13); the adapter's
  march identical to `analyze_endwalls` on the converged edges (|ΔK_B| = 0); Eq. 7-5 from the node data with K_B
  to 1e-4; purity and the refusals;
- a 15° cone with inlet swirl (C_θ = 20 m/s): φ handed over = 15° and m = the segment lengths (1e-12); the march
  identical to `analyze_endwalls` built from the contours by hand (0); K_B 0.980, 0.974, 0.968, 0.960, 0.952;
- BladeGen's curved walls: m and φ integrated along the hub and casing splines; φ agrees to 7e-16, m is 5.3e-4 m
  shorter in Ch. 7 (Findings 4), K_B against Ch. 8's march from the contours 3.8e-5.

**The synthetic stage** (IGV, rotor with a 0.5 mm tip gap, shrouded stator with a hub seal; 5000 rpm; every
row Ch. 6, the end walls Ch. 8), approximate and full model (1 and 9 passes), and approximate with Ch. 6's
Reynolds correction:
- mass error 1.1e-11; cp ΔT₀ = ω Δ(rC_θ) on every surface to 1.3e-14; I₂ = I₁ exactly; the stator does no work
  and turns the flow to within 5° of axial;
- the Ch. 6 values of the mid surface of every row rebuilt in the test from the converged nodes (μ by Sutherland
  at the node's T₁) and evaluated directly: β'₂, ω̄ and K_W agree to 1e-7 rad, 1e-8, 1e-8 (the flow tolerance).
  Measured (approximate model): IGV β'₂ −9.67°, ω̄ 0.0192; rotor β'₂ 44.88°, ω̄ 0.0324 (+ 0.0013 tip gap),
  D_eq 1.51; stator β'₂ −0.04°, ω̄ 0.0364 (+ 0.0005 seal), D_eq 1.67; with K_Re: rotor ω̄ 0.0309, stator 0.0357;
- Eq. 7-49 from the nodes to 9.2e-15; the tip-gap loss zero at the hub, the hub seal's zero at the shroud;
- K_B 0.980, 0.978, 0.976, 0.968, 0.960, 0.959, 0.958; the adapter's march equal to `analyze_endwalls` (0);
- K_B's effect: ∫K_W rρW_m cos ε dy with and without the end walls differs by 1/K_B to 1.7e-5;
- the mirrored stage (every angle, swirl, ω and s_row sign-changed): W_m, losses, K_W, K_B identical and the
  angles mirrored, to 0 (bit for bit).

**BladeGen's default rotor** (15000 rpm, 3 kg/s, its custom profile taken as NACA 65, a 0.3 mm tip gap; simple,
approximate and full model: 1, 1 and 799 passes): converged, mass error ≤ 1e-11, Euler work = cp ΔT₀ to
1.7e-14. Mid surface (ψ = 0.5): β'₁ 55.9° → β'₂ 43.3° (κ₂ 35.65°), ω̄ 0.036, P_t2/P_t1 1.044. K_B 0.980 before the
LE and 0.983 after the TE: the end-wall blockage falls across this rotor (its channel contracts and
accelerates the flow).

**The stream sheet** (the synthetic stage's rotor, surface 5): r, b, m equal the through-flow's at every
quasi-normal (exactly), linear between them, constant beyond; the crossings by hand (radial edge lines); the
carried thickness carries the tube's mass, 2π rρW_m b = ṁ/(n_ss − 1), to 2.2e-16, and is down to 0.959 of the
geometric one (K_B K_W); Eq. 12-39's thickness times K_B K_W likewise; the solver sheet's r inside the row equals
the linear rule to round-off and sin φ its slope to 1.4e-13; the inlet state equals the nodes at the
quasi-normals (Mach 1e-13), is linear between them and uses the radius where it is taken; the refusals.

**The case settings:** rpm 5000 recovered to 1e-14; the inlet state at the grid's own far inlet; the case's
relative inlet angle equals the triangle of the through-flow's absolute state there (1e-12); the thickness option
reaches the case's sheet; σ_TE by hand; the mirrored stage's rotor (s_row −1) gives the same case; the stator
(s_row −1) runs mirrored with a positive inlet angle. On a 15° conical inlet (the cone rotor, surface 5): the grid's
far inlet at M = −1.019 chords (nominal −1), r_in 0.2898 m against QN1's 0.3006 m; rC_θ conserved to 1e-13 at the
radius used, and the case's triangle at the same radius (2e-16) with the same C_θ and C_m (1e-9).

**One blade-to-blade solve against Chapter 6** (BladeGen's default rotor, approximate model; the surface whose
row span 0.4945 is nearest the Mid layer's 0.5, surface 4 of 11, ψ = 0.4; the Mid layer's blade): the solve's
station radii and b/b_in equal the sheet's (≤ 1e-6, measured 1.1e-8); σ_LE 1.0084 (ψ = 0.5 would give 1.093),
σ_TE 0.9768 on the through-flow's sheet against 0.9932 on the layer's own.

| | through-flow (Ch. 6) | blade-to-blade |
|---|---|---|
| β'₁ | 53.695° | 53.708° |
| β'₂ on QN2 (0.084 L behind the TE), very coarse / coarse / medium | 37.909° | 40.732° / 40.245° / 39.791° |
| β'₂ at the far exit (1 L behind the TE), very coarse / coarse / medium | – | 40.679° / 40.091° / 39.735° |
| deviation from the metal angle (very coarse) | 2.759° from κ₂ = 35.151° (the section at span 0.4945) | 5.026° from κ_TE = 35.653° (the Mid layer's) |
| W_m on QN2 (very coarse; carried / Mass / geometric b) | 96.03 m/s | 96.20 / 96.40 / 96.60 m/s |
| Euler work (very coarse) | 4105 J/kg | 3211 J/kg |

At matched metal angles the deviations differ by 2.3° (the blade-to-blade blade is the Mid layer's, whose κ_TE is
0.50° above Ch. 6's κ₂ at the surface's span; the default design's κ₂ runs −22.6°, 12.8°, 35.7°, 49.5°, 54.0° at spans 0, 0.25, 0.5, 0.75, 1, so 0.0055 of span is 0.5°). The difference of β'₂ on QN2 falls with the grid (2.8°, 2.3°, 1.9°).
No bound is set on it: it is a comparison of two models, not a test of either.

## Seeded mutations

Each mutation was applied alone, the tests `[coupling]`, `[correlations]` and `[throughflow]` (69 cases) were
built and run, and the file restored (the scripts are not committed). M1 to M23 are the first round's, R1 to R14
the independent review's (seven of them survived 2a5f5f7: R2 to R6, R10, R11), N1 to N4 new for this revision.
All 41 are killed.

| # | mutation | killed by (test cases) |
|---|---|---|
| M1 | s_row ignored on entry | the field-by-field, guide-vane, K_W, stage and sheet cases (8) |
| M2 | ω sign-changed into Ch. 6 | the field-by-field case, the stage's hand ω̄, the case settings |
| M3 | Re_c with μ at the relative total T | the field-by-field case, the stage with K_Re |
| M4 | o/c for o/s | the field-by-field case, BladeGen's rotor and its blade-to-blade cases (5) |
| M5 | K_B of the previous quasi-normal | the end-wall annulus and BladeGen-wall cases, the stage, the sheet cases (5) |
| M6 | loss = the minimum loss | the field-by-field case, the stage's hand ω̄ |
| M7 | the metal-angle placement sign | 11 cases |
| M8 | Eq. 12-39's b for the geometric b | the sheet case |
| M9 | the solver sheet's b divided by r | the sheet, case-settings and BladeGen blade-to-blade cases |
| M10 | γ fixed at 1.4 | the field-by-field case (γ = 1.3) |
| M11 | s/h with the local pitch | the field-by-field case, the stage |
| M12 | the tip gap with the hub section | the row-level loss case |
| M13 | the seal leakage unsigned | the row-level loss case |
| M14 | the shrouded wall keeps its clearance in Ch. 8 | the stage and the sheet cases (6) |
| M15 | the inlet C_θ with QN1's radius | the case settings, the conical inlet |
| M16 | the blade-to-blade rpm not mirrored | the case settings (the mirrored rotor) |
| M17 | the wall frame ω_w dropped in Ch. 8 | the stage (against `analyze_endwalls`), the case settings |
| M18 | the LE crossing's m from z | the sheet case |
| M19 | the far inlet across an upstream row | the sheet, inlet-state and case-settings cases |
| M20 | s_row applied again on exit | the field-by-field, guide-vane, mirrored-stage, stage and sheet cases |
| M21 | the Dunavant β_1c from Eq. 6-47's start | `test_correlations.cpp`'s guide-vane case |
| M22 | the throat not handed over by the through-flow | the stage, K_W and sheet cases (9) |
| M23 | the Dunavant model keeps Eq. 6-47's note | `test_correlations.cpp`'s guide-vane case |
| R1 | the sheet's β'₂ = the absolute α on QN2 | the sheet case |
| R2 | the through-flow's μ₁ at the relative total T | the stage with K_Re |
| R3 | the inlet state at QN1's radius | the case settings, the conical inlet |
| R4 | the thickness option ignored | the case settings, BladeGen's blade-to-blade case |
| R5 | K_W dropped (1) | the K_W case |
| R6 | the seal leakage into Ch. 8 unsigned | the Ch. 8 row-wall case (a negative leakage) |
| R7 | the hub section's throat everywhere | BladeGen's rotor, its blade-to-blade case, the curved-wall case |
| R8 | the stagger placement with a/c = 0.5 | the field-by-field case (parabolic arc) |
| R9 | no p. 128 flip in the stagger placement | the guide-vane case |
| R10 | the gap radius at QN1 only | the row-level loss case (a rising casing) |
| R11 | the Ch. 8 edge φ sign-changed | the cone case |
| R12 | the Ch. 8 edge m = z | the annulus, cone and BladeGen-wall cases |
| R13 | the hub step with the shroud row data | the stage, the case settings |
| R14 | the far inlet half a chord upstream | the case settings |
| N1 | the carried b without K_B | the sheet, case-settings and BladeGen blade-to-blade cases |
| N2 | the inlet at the nominal far inlet only | the case settings, the conical inlet, BladeGen's blade-to-blade case |
| N3 | σ_TE taken at the LE | the case settings |
| N4 | the Ch. 6 viscosity doubled | the field-by-field case, the stage with K_Re |

## Validation data

Chapters 7, 8 and 12 were searched for a published case with data. The book's results are Figs. 8-8 to 8-10
(end-wall blockage of the NACA 8-, 10- and 5-stage compressors by meridional station), Fig. 12-5 (blade loading
of the NACA 10-stage compressor's first stator at the hub) and Figs. 12-7, 12-8 (its end-wall layers). None
can be reproduced: the compressors' geometry and operating data (Voit 1953, Johnsen 1952, Budinger and Thomson
1952, Kovach and Sandercock 1961) are not in the book, and Ch. 7 has no example at all (Figs. 7-1 to 7-4 are
definitions and flow charts). The one number that can be compared, B = 0.02 at the first station of Fig. 8-9,
is the inlet value the march imposes (Eq. 8-72), reproduced exactly. Nothing was digitized.

## Findings and open items

1. **Open item for Ch. 12 (and Ch. 7): the K_W loop.** Where D_eq > 2, the plain substitution of the discharge
   quasi-normal (Ch. 7's flow-field loop, Fig. 7-3) does not settle. Eq. 6-84 gives K_W = (2/D_eq)^0.9; D_eq
   ≈ (W_max/W₁)(W₁/W₂) with W₂ ∝ W_m2, and Eq. 7-5 makes W_m2 ∝ 1/K_W, so D_eq ∝ K_W and the map
   K_W → (2/(cK_W))^0.9 has the slope −0.9 in ln K_W: at best a factor 0.9 per iteration, alternating.
   Measured on a single stator (22 vanes, α₁ = −40°, K_W 0.90 to 0.97 on the mid surface): a 2-cycle whose
   amplitude falls by 0.3 % per iteration at first and then holds (not converged in 40 000 iterations, the kink of
   Eqs. 6-83, 6-84 at D_eq = 2 reached at the outer surfaces). In the first stage design (44 vanes) the same loop
   failed with and without the end walls depending on the design. The test stage is designed below D_eq = 2, and
   K_W's way into Eq. 7-5 is tested with the discharge frozen (test wrapper `FrozenDischarge`). An
   under-relaxation of K_W (or of W_m2 in the row model's input) would be ours; p. 195's "some limits and
   numerical damping are normally required" gives no form.
2. The Ch. 6 type `InletState` clashed with Flow.hpp's (an ODR violation): renamed `CascadeInletState` (above).
3. The through-flow did not hand the throat to the row model (fixed, above).
4. The wall arc length m that Ch. 7 hands to Ch. 8 is Eq. 7-58 integrated over the quasi-normal points only. On
   BladeGen's default rotor, whose hub rises inside the row between two quasi-normals, it is 5.3e-4 m (1.3 %)
   shorter than the contour's own arc; K_B moves by 3.8e-5 (measured against Ch. 8's march from the contours).
5. Ch. 8's row step (`endwall_step`) does not converge on the 15° cone rotor of the tests ("the blade force defect
   did not converge in 200 iterations"); the cone test of the adapter therefore runs without the rotor. Not
   investigated (Endwall.cpp is not this branch's module).
6. The blade-to-blade and Ch. 6 discharge angles differ by about 2.8° on the mid surface (below): by the
   independent review, the Euler-work difference follows from it (about 312 J/kg per degree) and is a model
   difference, not an adapter defect.

## Not verified

- The book gives no example of any adapter's mapping. Every mapping marked ours is a reading, not a check.
- The blade-to-blade comparison is on the very coarse to medium grids only and on one surface; the section is
  BladeGen's layer, not re-placed on the through-flow surface (σ_LE 1.0084, σ_TE 0.9768 on the through-flow sheet
  against 0.9932 on the layer's own).
- The Ch. 8 adapter is compared with Ch. 8's own march on cylinders, a 15° cone (exact m and φ) and BladeGen's
  curved walls (m and φ integrated along the contour splines); not with any published boundary layer.
