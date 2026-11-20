"""Reference data from the original blade2d Python engine (not part of the build).

Builds the NACA 65-(12)10 fit with blade2d's profile_c2_connect.build_profile_c2_connection
(the engine BladeGen was ported from) and writes, to tests/data/blade2d_naca65_reference.csv:
  kind,i,x,y   with kind in {ss_ctrl, ps_ctrl, ss, ps, le, te}
The C++ test (tests/profile/test_blade2d_reference.cpp) rebuilds the same profile in
ConnectorMode::G2Legacy and compares.

Run with a Python that has numpy + scipy, from anywhere:
    python tools/reference/make_blade2d_reference.py <path-to-blade2d>
"""
import sys
from pathlib import Path

import numpy as np

blade2d = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\ASUS\Documents\blade2d")
sys.path.insert(0, str(blade2d))
from geometry.characteristicpoints import EdgeParams  # noqa: E402
from geometry.ssbspline import SuctionSideSplineParams  # noqa: E402
from geometry.pscoupling import PressureSideCouplingParams, ThicknessDistributionParams  # noqa: E402
from geometry.profile_c2_connect import ProfileC2ConnectParams, build_profile_c2_connection  # noqa: E402

le = EdgeParams(radius=0.012073, stretch_ab=1.198788, asym=0.003422, beta_deg=15.561637,
                d_beta_ss_deg=22.417183, d_beta_ps_deg=15.638344)
te = EdgeParams(radius=0.000255, stretch_ab=4.0, asym=-0.000526, beta_deg=24.659155,
                d_beta_ss_deg=6.386629, d_beta_ps_deg=0.0)
N = 2000  # dense SS sampling: blade2d's arc length / normals come from samples
ss = SuctionSideSplineParams(le_params=le, te_params=te, handle_le_frac=0.092724, ctrl2_mprime=0.52168,
                             ctrl2_theta=0.177553, handle_te_frac=0.041191, n_samples=N)
ps = PressureSideCouplingParams(ss_params=ss,
                                thickness_params=ThicknessDistributionParams(0.11027, 0.405561, 0.035922,
                                                                             0.002126, 0.685326, 1.145254),
                                handle_le_frac=0.006057, handle_te_frac=0.05123, n_samples=N)
g = build_profile_c2_connection(ProfileC2ConnectParams(ss_params=ss, ps_params=ps, n_edge_samples=201,
                                                       n_ps_samples=N))
out = Path(__file__).resolve().parents[2] / "tests" / "data" / "blade2d_naca65_reference.csv"
rows = []
for kind, arr in (("ss_ctrl", g.ss_geom.ctrl), ("ps_ctrl", g.ps_geom.ctrl),
                  ("ss", np.column_stack([g.ss_geom.curve_x, g.ss_geom.curve_y])[::20]),
                  ("ps", np.column_stack([g.ps_curve_x, g.ps_curve_y])[::20]),
                  ("le", np.column_stack([g.le_x, g.le_y])), ("te", np.column_stack([g.te_x, g.te_y]))):
    rows += [f"{kind},{i},{x:.15g},{y:.15g}" for i, (x, y) in enumerate(arr)]
out.write_text("# blade2d profile_c2_connect reference (NACA 65-(12)10 fit), make_blade2d_reference.py\n"
               "kind,i,x,y\n" + "\n".join(rows) + "\n")
print(f"wrote {out} ({len(rows)} rows)")
