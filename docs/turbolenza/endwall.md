# End-wall boundary layers (Aungier 2003, Ch. 8)

`src/turbolenza/Endwall.hpp/.cpp` computes the hub and shroud (casing) boundary layers of an
axial compressor along the end walls. It follows R. H. Aungier, *Axial-Flow Compressors* (ASME Press,
2003), Ch. 8, printed pp. 175-198. It returns the end-wall blockage factor K_B of Eq. 7-5 to the
meridional through-flow (Ch. 7).

This module is not a v1 port. As for the blade-surface boundary layer, every equation cites the
number read from the page images, and every choice the book does not make is marked "ours" in the
header (`Endwall.hpp`). The header lists every equation as implemented. This page documents the
interface.

## What goes back to the through-flow

The module returns only **K_B = 1 - B per station** (Eqs. 8-76 and 8-77). B is the sum over hub and
shroud of 2 pi r rho_e V_me delta1* / mdot, and K_B multiplies the whole integral of Eq. 7-5.

Nothing else is returned:
- The blade force defects nu stay inside the end-wall momentum integrals 8-6 and 8-7.
- End-wall losses and work come from Ch. 6, not from this module (p. 175).
- The full layer is reported per station and wall, for reports and for the alternative of p. 164
  (walls displaced by delta1*), which Aungier does not use.
- The book prescribes no damping of K_B in the coupled iteration ("some limits and numerical damping
  are normally required", p. 195, with no formula). The module applies none.

## Calling it (Fig. 7-3: once per quasi-normal, repeatable)

```cpp
// station 0, each wall
EndwallLayer h = endwall_inlet(hub_wall[0], inlet, mdot, settings);
// station k >= 1, each wall; row = nullptr between rows, &row_wall when k-1 -> k spans a blade row
EndwallStepResult s = endwall_step(h_prev, hub_wall[k-1], hub_wall[k], row, omega_c, settings);
// K_B of station k
EndwallBlockage b = endwall_blockage(s_hub.layer, hub_wall[k], s_shroud.layer, shroud_wall[k], mdot);
// b.K_B goes into Eq. 7-5 at station k
```

`endwall_step` is pure. The same arguments always give the same result, bit for bit, so the
flow-field loop may call it again at the same station from the converged layer at k-1.
`analyze_endwalls(EndwallInput)` is the convenience march over all stations; the tests use it.

## Conventions

- **SI units, radians, rad/s.**
- **Book signs.** omega > 0 turns toward +theta, and W_theta = C_theta - omega r (8-67). TURBOlenza
  maps to the book with C_theta,book = C_theta,T and **omega_book = -2 pi rpm / 60**; the q3d data
  contract, Sec. 0.3, is binding. Every equation is then used verbatim.
- **Velocities enter in the absolute frame** (V_m and C_theta). The module forms the wall frame
  itself, V_theta = C_theta - r omega_w (8-4), and the blade-row frame, W_theta = C_theta - r omega.
- **m is the arc length along each wall**, so the hub and shroud have their own m. It must strictly
  increase with the station index.
- **wall_omega at station k** is the rotation speed of the wall segment from k-1 to k. At station 0 it
  is the frame of the inlet data. A change of wall frame therefore always falls at a station, and the
  layer at station k is reported in the frame of the segment that ends there.
- **A blade row in the Ch. 7 through-flow is one step**, from its leading-edge station k-1 to its
  trailing-edge station k (p. 303). `analyze_endwalls` also accepts stations inside a row. With the
  Ch. 8 model the row then keeps one nu over all its intervals (ours).

## Input, per station and wall: `EndwallWallStation`

| field | meaning | unit | notes |
|---|---|---|---|
| `z` | axial coordinate of the wall point | m | for the record only |
| `r` | radius of the wall point | m | > 0 |
| `phi` | wall contour angle with the axial direction | rad | sin(phi) = dr/dm |
| `m` | arc length along this wall | m | strictly increasing in k |
| `wall_omega` | omega_w of the wall segment from k-1 to k (station 0: inlet frame) | rad/s | book sign; drum = rotor speed, casing = 0 |
| `Vm` | edge meridional velocity C_m = W_m = V_m (8-3) | m/s | > 0 |
| `Ctheta` | edge absolute tangential velocity | m/s | book sign |
| `rho` | edge static density | kg/m^3 | > 0 |
| `P` | edge static pressure | Pa | used in 8-8, inside rows only |
| `mu` | viscosity at the edge | Pa s | > 0; e.g. `sutherland_viscosity(T_e)` |

"Edge" means the inviscid end-wall stream surface (j = 0 or N) at the quasi-normal foot. Static
temperature is not needed, because mu is passed directly.

For the march, `EndwallStation` holds `{hub, shroud, quasi_normal_length}`. The quasi-normal length
(m) is used only by the diagnostic `K_B_geometric`.

## Input, per blade row and wall: `EndwallRowWall`

| field | meaning | unit | notes |
|---|---|---|---|
| `omega` | blade rotation speed | rad/s | book sign; 0 for a stator; the same on both walls |
| `blade_count` | Z | - | not used by Ch. 8 (kept for Ch. 12's 12-48) |
| `pitch` | s bar at this wall: the mean of inlet and discharge | m | > 0 (8-2) |
| `stagger` | gamma bar at this wall | rad | abs < pi/2 (8-2) |
| `clearance` | delta_c: tip or hub gap at this wall | m | >= 0; 0 for an attached blade; must be 0 on a shrouded wall (refused otherwise) |
| `shrouded` | shroud band with a seal at this wall | - | across the row the wall must rotate at `omega` (p. 193) |
| `seal_leakage` | mdot_leak from Ch. 6 | kg/s | positive from row discharge to row inlet; 0 unless shrouded |

For the march, `EndwallBladeRow` holds `{station_in, station_out, hub, shroud}`.

## Input, global

- `mass_flow`: the compressor mdot in kg/s, used in 8-72 and 8-76.
- `compressor_omega`: omega_c in rad/s, book sign, used in 8-61.
- `EndwallInlet` (Sec. 8.8):
  - `blockage` B. Default 0.02, the value for an analysis started just upstream of the first row. A
    very small B starts the layer at transition through 8-73.
  - `H1`, default 1.4.
  - `m`, default 1/7.
- `EndwallSettings`:

| setting | default | meaning |
|---|---|---|
| entrainment | `FreeStreamH` | 8-43; `MeridionalH1` selects 8-42 |
| force-defect model | `RowChapter8` | Sec. 8.5; `None` gives nu = 0 |
| `force_defect_correction` | on | 8-57 |
| shape-factor limit | 2.4 | text, pp. 184 and 186; Fig. 8-4 draws 2.2 |
| Re_theta limit | 250 | pp. 186 and 194 |
| theta12 precaution | 0.99 | p. 195 |
| `extended_m_min` | -0.85 | ours: lower limit of m (H2 >= 0.40, the end of Fig. 8-4's curve); keeps the layer continuous through V_theta_e = 0 |
| tolerance | 1e-12 | ours |
| `max_iterations` | 200 | ours |

## Output, per station and wall: `EndwallLayer`

The layer is reported in the frame of the segment ending at that station.

| field | meaning | unit |
|---|---|---|
| `Q`, `theta11`, `q12` | the marched, frame-invariant state: delta - delta1*, theta11, V_theta_e theta12 (8-68) | m, m, m^2/s |
| `wall_omega`, `Vtheta_e`, `V_e`, `beta_e` | the frame, and the edge relative to the wall | rad/s, m/s, m/s, rad |
| `delta`, `delta1`, `delta2` | delta, meridional delta1*, tangential delta2* | m |
| `theta12`, `theta22` | tangential momentum thicknesses in this frame | m |
| `n`, `m`, `H1`, `H2` | profile exponents and shape factors, after the 2.4 limits | - |
| `delta_star`, `theta`, `H` | streamwise (8-40, 8-41), H limited | m, m, - |
| `Re_theta`, `Re_theta11` | rho V_e theta / mu and rho V_e theta11 / mu | - |
| `cf`, `E`, `tau_mw`, `tau_thetaw` | 8-36, 8-42 or 8-43, 8-34, 8-35 | -, -, Pa, Pa |
| `nu`, `dnu` | force defect of the interval ending here (0 outside rows) | m |
| flags | `H1_limited`, `H2_limited`, `theta12_limited`, `H_limited`, `re_theta_limited` (8-36), `re_theta11_limited` (8-73) | - |
| `iterations` | substitution passes at this station | - |

## Output, per station: `EndwallBlockage`

- `B_hub` and `B_shroud`: each wall's share.
- `B`: Eq. 8-76.
- `K_B`: Eq. 8-77, **the value for Eq. 7-5**.
- `K_B_geometric` (march only): 1 - (delta1*_hub + delta1*_shroud) / quasi-normal length. This is a
  diagnostic (ours).

## Output, per row and wall: `EndwallRowWallResult`

| field | meaning |
|---|---|
| `g` | 8-2 |
| `theta11_mean` | theta11 bar |
| `nu0`, `dnu`, `nu` | 8-44, 8-57, 8-56 |
| `F` | 8-55 |
| `dnu_bound`, `dnu_limited` | the bound abs(dnu) <= [0.05 + 0.2 (H1 - 1)] delta at the discharge, Aungier's Eq. 12-50, applied to the Ch. 8 row correction (ours), and whether it acted |
| `xi` | 8-53 |
| `collateral_defect` | V_theta_e theta12 - W_theta_e theta11 at the discharge (8-45) |
| `iterations` | row iterations |
| `leakage` | Delta+, Delta-, U_leak at the inlet and discharge, and the states before and after each correction (8-58..8-66) |

## Reading of the book (decisions)

- **8-30 is used as printed**, theta22/delta = 0.95 m - 1.684 m^2. It is a fit, not the integral of
  profile 8-28. The integral is -4.39e-4 + 0.9502 m - 1.5126 m^2, so 8-30 is off by 3.5 % at m = -0.3
  and 6.5 % at m = -0.85 (measured in the tests). 8-29 and 8-31 are exact.
- **Separation limit 2.4**, the text value (pp. 184 and 186). Fig. 8-4 draws its limit line at
  H2 = 2.2.
- **Two Reynolds limits.**
  - The limit on Re_theta of the streamwise theta (p. 186) is applied inside c_f (8-36).
  - The limit on Re_theta11 (8-73) is applied on the marched state. A station below it is scaled up,
    keeping H1 and m.
- **The procedure's m < 0.05 branch** is used. Steps 4 and 5 name only 8-21 and 8-22..8-26, but Sec.
  8.3 says 8-31, 8-29, 8-30 and 8-32 apply below m = 0.05.
- **8-6 and 8-7 are used as printed.**
  - Ch. 3's differential forms on p. 55 (3-46..3-48) are misprinted: r is missing inside d/dm, and the
    sign of d tau/dy is wrong.
  - The integral forms, which this module uses, are consistent.
  - 8-6 omits a centrifugal density-defect term, which is zero for the constant density that Sec. 8.3
    assumes.
  - 8-69..8-71 hold for constant density only.
- **The symbol B in 8-48 and 8-50** is an integration term, not the blockage. 8-47, 8-48 and 8-50 are
  not coded: they serve only to derive xi (8-53).
- **8-45, 8-62 and 8-66** print V_theta and W_theta without the subscript e. Edge values are used.
- **E (8-43) uses the limited H** (ours; p. 186 states the 2.4 limit only for 8-36).
- **The step-1 hold of H2 is not used.** The march carries q12, so 8-75 already fixes theta12 and
  with it m. The start does not change the converged station.
- **The 8-71 bracket is not carried across a frame change.** It is rebuilt from the profile in the
  new frame. The profile family misses 8-69 by -2.9 % and 8-71 by -3.0 % at a
  stationary-to-rotating hub (measured).
- **V_theta_e -> 0 with q12 != 0.** theta12 = q12 / V_theta_e goes to +-infinity there. The H2 limit
  (m <= 0.7) and our lower limit (m >= -0.85) keep the profile finite, so the 8-6 bracket and every
  tangential term go to 0 from both sides. The layer is then continuous through V_theta_e = 0
  (tested). Without the lower limit, 8-31 gives m -> -infinity on one side and the bracket jumps.
- **The force defect as the turning vanishes.** 8-57 alone gives dnu ~ 1/xi, a force defect with
  no blade force. Bounded by 12-50, the defect terms vanish with the blade force (tested).

## Not done

- Ch. 12's in-passage force-defect model (12-44..12-50) is not implemented. The settings enum has room
  for it, and the march keeps one nu per interval. Eq. 12-46 is misprinted (its brackets do not
  balance).
- The typical results of Sec. 8.9 (Figs. 8-8 to 8-10) cannot be reproduced. They need the NACA
  compressors' through-flow solutions, which the book does not give. The only comparable number is
  the inlet B = 0.02 that Fig. 8-9 starts from, and the march reproduces it exactly.
