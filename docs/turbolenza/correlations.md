# TURBOlenza: the empirical cascade performance models (Aungier 2003, Ch. 6)

`src/turbolenza/Correlations.hpp/.cpp` implements the empirical cascade models of R. H. Aungier,
*Axial-Flow Compressors* (ASME Press, 2003), Ch. 6, printed pp. 117-152. The meridional through-flow
(Ch. 7, Sec. 7.3) uses them for each blade row on each stream surface: they give the discharge relative
angle beta'_2, the loss coefficient and the wake blockage K_W. The quasi-3D analysis (Ch. 12) uses them
later. This is not a v1 port. Rules 2 to 6 of `PORT.md` apply. In place of rule 1, every equation cites
its number as read from the page images, and every choice the book does not make is marked "ours".

## Sources

- Ch. 6, all equations 6-1 to 6-102, read from the page images. A second, independent reading of the
  same pages agrees with every one of them.
- Ch. 4, pp. 66-76: the blade angles of the circular arc (4-7, 4-35), the parabolic arc (4-17 to 4-20)
  and the NACA C_l0-to-camber relation (4-29 = 6-7). The parabolic arc is also used for A4K6 (p. 66).
- Ch. 7, Sec. 7.3, pp. 162-164: how the through-flow uses the models (7-43 to 7-50). Eq. 7-6, p. 157:
  K_W = 1 - B_wake.
- Aungier, *Centrifugal Compressors* (2000), Ch. 4, pp. 59-64: the seal coefficients (Eqs. 4-34 to
  4-38), the 0.816 of Eq. 4-17, and the statement that rho and T of Egli's model are taken on the
  higher-pressure side (p. 63).

## Conventions

These are the interface to the through-flow (Ch. 7) and the quasi-3D analysis (Ch. 12). Full text: the
file comment of `Correlations.hpp`.

**Input frame.** Inputs are in the book's Ch. 7 frame:
- Angles are measured from the meridional direction and signed with +theta.
- W_theta = C_theta - omega r, with omega signed.

**Row sign.** Each row carries s_row = `CascadeGeometry::row_sign`.
- Every signed input is multiplied by s_row. This gives the positive-camber cascade frame of Fig. 6-1.
- A row whose camber is still negative after that is refused, unless it is an inlet guide vane. A guide
  vane gets Aungier's own sign change of p. 128 on top.
- `frame_sign` reports the total factor used.
- `beta2` is returned in the book frame. All `_deg` fields are in the cascade frame.

**Units.**
- Structs use radians.
- Functions and fields named `_deg` use degrees, because the correlations are dimensional in degrees.
  The book states this on pp. 135 and 144 and implies it elsewhere.
- Everything else is SI.

**Loss coefficient.** omega_bar = Delta P't / (P't1 - P1) at QN1 (Eq. 7-49).
- The shock Delta P't is converted with the perfect-gas P't1 - P1 at M'_1 (ours).
- Row-level leakage losses are pressures. `leakage_loss_at` distributes them and
  `loss_coefficient_from_pressure` converts them (ours: the book gives no conversion).

**Reynolds number.** Re_c = rho1 W1 c / mu1 (ours: the book gives no formula). It is needed only with the
Reynolds correction, which Aungier does "not normally" use (p. 150).

**Gas.** A perfect gas (gamma, R) is used for choke (6-53), W* (6-59, 6-60), the shock (6-70 to 6-75) and
the seal (6-91). Aungier uses the general equation of state of Ch. 2. A real-gas version is not
implemented.

## API

`evaluate(CascadeGeometry, InletState, DischargeEstimate, Gas, CorrelationOptions) -> CascadePerformance`.
This is the per-stream-surface call of Sec. 7.3.
- Inputs:
  - QN1: r1, W_m1, W_theta1, omega, M'_1, and rho1, mu1 for Re_c.
  - QN2: r2, W_m2, the through-flow's current iterate.
- Outputs:
  - beta'_2, omega_bar (design, minimum, incidence and shock parts), B_wake and K_W.
  - i, i*, delta, delta*, alpha*, beta1*, beta2*.
  - D, D*, D_eq, D_eq*, W_max/W1.
  - R_c, R_s, beta_1c, beta_1s, i_c, i_s, i_m, beta_1,choke, the stall margins, xi and M'_c.
  - Flags (`choke_limit_active`, `design_diffusion_limit_exceeded`, `extrapolated`) and `notes`.
- Options:
  - Design incidence: alpha* (6-8, the default, Aungier's choice on p. 125) or Lieblein (6-12 to 6-18).
  - Profile loss: D_eq (6-46, the default, p. 134) or D (6-45).
  - K1, K2 calibration: 2-D cascade (6-34 / 6-37), constants, or Sec. 6.11 (K1 = 0.0073 with K2 from
    6-99, the default).
  - Reynolds correction (6-100 to 6-102), off by default.
  - alpha of 6-38 for DCA.

`endwall_leakage(RowWallData, WallState qn1, WallState qn2, mdot, Gas) -> RowLeakage` is the row-level
call:
- Tip clearance: 6-85 to 6-89.
- Shroud seal: 6-90 to 6-96. For N > 12 the carryover coefficient uses Aungier (2000) Eqs. 4-36 and
  4-37 (p. 64), where the 2003 text refers (p. 149), instead of refusing.
- `leakage_loss_at(DeltaP_t, psi, gap_side)` gives the linear distribution of p. 147.

Every single correlation is also a public function (degrees, cascade frame), named for its equation in
the header.

## Equations, as read

| Eqs. | p. | Implemented as | Reading notes |
|---|---|---|---|
| 6-1 .. 6-6 | 119-120 | definitions in `evaluate` | |
| 6-7 (= 4-29) | 120 | `naca_camber_from_cl0_deg` | 0.05515 confirmed in Ch. 4, p. 71 |
| 6-8, 6-9 | 121 | `design_attack_angle_deg` | the exponent 0.25 applies to a/c only; sigma^e multiplies the bracket |
| 6-10, 6-11 | 122 | `incidence_thickness_factor` | |
| 6-12 .. 6-15 | 122-123 | `lieblein_zero_camber_incidence_deg`, `lieblein_incidence_slope`, `lieblein_design_incidence` | fixed point on beta1 (ours: substitution to 1e-12 deg) |
| 6-16 .. 6-18 | 124 | `lieblein_design_incidence` | used for every camberline (for a/c = 0.5 and gamma >= 0 it is 6-12 itself) |
| 6-19 .. 6-25 | 125-126 | `lieblein_*deviation*`, `deviation_thickness_factor` | m_1.0 by camberline: 6-22 NACA 65, 6-23 circular arc |
| 6-26 | 128 | `howell_design_deviation_deg` | used for IGVs, the parabolic arc and A4K6 |
| 6-27, 6-29, 6-30 | 128-130 | `wake_loss_coefficient`, `*_from_velocities` | 6-28 is a proportionality and is not implemented |
| 6-31 | 130 | `diffusion_factor` | |
| 6-32, 6-33 | 130 | not implemented | correlation forms only |
| 6-34 .. 6-37, 6-45, 6-46 | 130-133 | `design_loss_*`, `lieblein_wmax_ratio_design`, `equivalent_diffusion` | 6-34 / 6-37 are 6-45 / 6-46 with K1 = 0.0035 / 0.004 and K2 = 1 |
| 6-38, 6-39 | 132 | `wmax_ratio` | see the decisions on 6-39 below |
| 6-40 .. 6-44 | 133 | `howell_*` | Howell's background; not used by `evaluate` (as Aungier) |
| 6-47 .. 6-50 | 135-136 | `*_stall_attack_offset_deg`, `stall_ranges` | 6-48 is linear in beta_1s, so it is solved in closed form |
| 6-51 .. 6-60 | 136-138 | `mach_*`, `choke_inlet_angle_deg`, `minimum_loss_incidence_deg`, `critical_mach_number`, `minimum_loss` | |
| 6-56, 6-57 | 138 | equal to 6-80 for abs(xi) <= 1 | |
| 6-61 .. 6-75 | 139-141 | `shock_expansion_geometry`, `prandtl_meyer_*`, `normal_shock_total_pressure_ratio`, `shock_loss` | 6-61 and 6-63 are misprinted (below); the perfect-gas closed forms are used for 6-70 and 6-72 .. 6-75 |
| 6-76 .. 6-84 | 142-145 | `deviation_incidence_slope`, `off_design_deviation_deg`, `normalized_incidence`, `off_design_loss_factor`, `wake_blockage` | |
| 6-85 .. 6-89 | 146-147 | `tip_clearance_loss`, `endwall_leakage` | |
| 6-90 .. 6-96 | 148-149 | `seal_*`, `shroud_seal_loss`, `endwall_leakage` | 6-93 is misprinted (below); N > 12 uses Aungier (2000) Eqs. 4-36 and 4-37 |
| 6-97, 6-98 | 150 | not implemented | the argument behind 6-99 |
| 6-99 .. 6-102 | 150 | `aspect_ratio_k2`, `reynolds_factor` | |

## Two misprints, and how they are used

- **Eqs. 6-61 and 6-63 (pp. 139-140).** They are printed "2R/c = sin(theta/2)".
  - The book's own Eq. 4-7 (p. 66) reads c/2 = R sin(theta/2). That is the circular arc through both
    blade ends, with R >= c/2.
  - The code uses c/(2R_u) = sin(theta_u/2).
  - A test checks R_u >= c/2 and shows that the printed form would give R_u < c/2.
- **Eq. 6-93 (p. 149).** It is printed C_t = (2.143 ln N - 1.464)/(N - 4.322) [1 - P_R]^(0.375 P_R).
  - In that form C_t is negative for N = 2 to 4 (-0.67 at N = 3) and has a pole at N = 4.322.
  - The source, Aungier (2000) Eq. 4-35 (p. 64), reads 2.143 [ln(N) - 1.464]/(N - 4.322). That form
    gives 0.944, 0.711, 0.592 and 0.517 at N = 1 to 4, which matches Egli's curves in Fig. 4-8 there.
  - The code uses the 2000 form.
  - Its numerator vanishes at e^1.464 = 4.3230, 0.001 from the pole. This near-cancellation never
    matters for integer N.

## Silences and ambiguities: what is done

- **6-39, the C_theta sign.** 6-39 is written with absolute C_theta, as printed. It is evaluated in
  the cascade frame, with C_theta = W_m tan(beta) + omega r and omega sign-changed with the row. This
  makes r1 C_theta1 - r2 C_theta2 positive both for a stator and for a work-adding rotor. Rotation
  enters only through (U1/W_m1)(1 - (r2/r1)^2). A test checks that term.
- **(delta*_0)_10 and (i*_0)_10 for beta1 < 0** (Eqs. 6-20, 6-13). After the p. 128 sign change a
  guide vane usually has a slightly negative beta1*, and Eq. 6-26 (used for ALL guide vanes, p. 128)
  needs (delta*_0)_10 there. Both functions are defined for beta1 < 0 by their odd extension
  f(beta1) = -f(-beta1) (ours). The basis: they describe a cascade of zero camber as a function of beta1
  and sigma only. Mirroring a zero-camber cascade gives a zero-camber cascade with beta1, kappa1, kappa2
  and beta2 sign-changed, so its deviation and incidence change sign. 6-20 is 0 at beta1 = 0, so its
  extension is continuous, with the same slope. 6-13's second term leaves a jump of about
  5e-9 sigma^3 deg. The argument does not hold for the camber slopes n and m (6-15, 6-21), which multiply
  a positive camber, nor for 6-76. Those still refuse beta1 < 0.
- **6-38 and 6-39 for i < i*.** The book states alpha (i - i*)^1.43 "for i >= i*". For i < i* the term
  is omitted (ours), and a note is written. `wmax_ratio` refuses a negative excess.
- **6-38 and 6-39 for DCA.** The book gives alpha only for NACA (0.0117) and C.4 (0.007). For DCA,
  `evaluate` refuses at i > i* unless `CorrelationOptions::lieblein_alpha` is set.
- **The design state.** It uses the current W_m2/W_m1, r2/r1 and U1/W_m1, with W_max from 6-39 at
  i = i* (ours, after p. 133: "for both design and off-design operation"). For a planar stationary
  cascade this is 6-35.
- **beta1 in 6-76.** It is read as beta1*, because the slope is "at the design incidence angle" (ours).
- **beta2 in 6-99 and 6-100.** It is taken as beta2*, because K2 belongs to the design loss 6-46 (ours).
- **The Lieblein route.** The pseudo-blade route 6-16 to 6-18 is given for A4K6 and parabolic arcs
  (p. 123). It is used for every camberline and stagger sign (ours). For a/c = 0.5 and gamma >= 0 it is
  6-12 itself, and |gamma| keeps it defined for gamma < 0.
- **A4K6 stall ranges.** p. 136 quotes Dunavant's nearly constant A4K6 ranges, alpha_s - alpha* = 10 deg
  and alpha_c - alpha* = -10 deg. Aungier computes the stall incidence of all blade sections with
  6-47 and 6-48, and attributes the constancy to the guide vanes' constant inlet angle. So 6-47 and 6-48
  are the default, with a note for A4K6. `StallRangeModel::DunavantA4K6` gives R_c = R_s = 10 deg
  (A4K6 only), and its stall inlet angles beta_1c = gamma + alpha* - 10 deg and beta_1s = gamma + alpha* + 10 deg;
  6-47 and 6-48 are not evaluated then, so 6-47's beta_1c >= 20 deg limit and its note do not apply.
- **6-52.** K_sh <= 1, as stated, and reported. In 6-60, K_sh is used as defined: the constraint is
  stated for 6-52 only (ours).
- **6-60.** M' is M'_1 limited to 1 (p. 138). W*/W1 is taken at the actual M'_1.
- **6-58.** R_s is the low-speed range (6-50). Without the choke limit, i_m >= i* follows from 6-51,
  6-52 and 6-55 because K_sh <= 1 there.
- **6-54.** The equation's +1 deg is used. The prose's "2 % from choke" is not separately modelled.
  If i_c rises to i_s or beyond, the input is refused: no low-loss range is left.
- **6-71.** The geometric mean is used, as the equation gives it. The prose says "average".
- **6-47.** beta_1c >= 20 deg is Aungier's own limit (p. 135). The equation is evaluated at 20 deg and
  a note is written. beta_1c is found by substitution. It converges because the map is a contraction:
  its slope is largest at beta_1c = 20 deg, where it is 0.0070 theta (0.63 at theta = 90 deg). Camber is
  restricted to [0, 90) deg.
- **The bow shock.** For M'_1 <= 1 there is no shock loss, and the Fig. 6-17 geometry is not built.
- **The seal at P1 = P2.** C_t = 0 at P_R = 1, so the leakage and the loss are 0, and the direction is 0.
- **6-85.** |tau| is used. A stator, which removes swirl, gives tau < 0, and 6-87 would then be
  imaginary.
- **rho_bar in 6-87 and 6-88.** (rho1 + rho2)/2 in `endwall_leakage` (ours).
- **6-90 and 6-91.** For the seal: Delta P = |P2 - P1|, P_R = P_low/P_high, and rho, T from the
  higher-pressure side (Aungier 2000 p. 63).
- **The p. 147 distribution.** Linear in the mass fraction, zero at the opposite wall, with a
  mass-weighted mean equal to the row value (ours: the book says only "a linear distribution across
  the annulus").
- **H in 6-74.** The total enthalpy, implicit in the perfect-gas normal-shock relations.
- **The Fig. 6-6 legend.** "t_b = 0.1" means t_b/c = 0.1.
- **Data ranges.** Ranges beyond the figures' data are reported, not refused: beta1* outside 0 to
  70 deg, sigma outside 0.4 to 2.0, t_b/c > 0.12.
- **The leading-edge radius.** Not needed. Sec. 6.7 assumes a sharp leading edge (p. 141).
- **The throat o/s.** Required (6-53). It is the caller's, for example from BladeGen's
  `make_cascade`.

## Figures

No correlation of the chapter exists only as a figure. Each data figure (Figs. 6-2 to 6-9, 6-11, 6-12,
6-14, 6-18) shows the source curve with Aungier's equation plotted over it as circles. The tests check
the equations against those circles, read by eye from the page images to within the stated margins
(about 0.02 of the plotted unit, 0.3 deg in Fig. 6-2). Nothing was digitized.

The remaining figures carry no equation:
- Fig. 6-19 (Mellor) is qualitative. Aungier abandoned it in favour of 6-77.
- Fig. 6-23 is experimental data behind 6-101 and 6-102.

## Tests

`tests/turbolenza/test_correlations.cpp` covers every correlation:
- against hand-computed values written out in the test;
- against the book's figure circles;
- in the limits the book states (zero camber; K_t,i = K_t,delta = 1 at t_b/c = 0.1; i = i* gives
  delta* and omega*; omega = 2 omega_m at i_c and i_s; K_Re = 0 at 2.5e5);
- for the monotonic trends;
- for the misprints and the refusals.

Further checks:
- Mirroring the row (s_row = -1) gives the mirrored answer.
- An IGV's deviation does not depend on incidence.
- Eq. 6-70 integrated numerically equals the closed-form Prandtl-Meyer angle.

`tests/turbolenza/test_correlations_dunavant.cpp` runs Dunavant's NACA 65-(12A10)10 cascade as Aungier
uses it in Ch. 5 (Fig. 5-6, p. 91): beta1 = 45 deg, alpha = 12 deg, stagger 33 deg, sigma = 1, at
M = 0.65 and M = 0.3. The throat o/s = 0.7215 is BladeGen's geometric throat of the section. The planar
calibration (6-37) is used.

| | Ch. 6 | the book |
|---|---|---|
| alpha* (6-8) | 12.557 deg (Lieblein: 13.890) | alpha = 12 deg is "within 0.5 deg of the optimum" (p. 91): 0.557 deg |
| i*, delta* | -2.522 deg, 7.781 deg | - |
| beta2 at beta1 = 45 deg | 25.654 deg | measured 25.0 deg (Fig. 5-6) |
| D*, D_eq*, D, D_eq | 0.411, 1.649, 0.399, 1.630 | - |
| omega* = omega_m; omega at M 0.65 / 0.3 | 0.01242; 0.01249 / 0.01248 | - |
| i_c, i_m, i_s at M 0.65 | -9.81, -2.52, 4.81 deg (choke beta1 34.98 deg, not active) | - |

AVDR is taken as 1: Dunavant's axial-velocity density ratio is not in the book.
