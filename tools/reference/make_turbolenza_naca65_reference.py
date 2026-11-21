"""Reference NACA 65 sections from TURBOlenza (not part of the build).

Builds the frozen C2 geometry with turbolenza.b2b.geometry.make_cascade (stagger 0, so the
cascade frame is the blade frame) and writes, to tests/data/turbolenza_naca65_reference.csv:
  case,cl0,t,variant,side,x,y        side: SS / PS, LE -> TE, sampled on the C2 curves
plus per case the LE / TE radius and the maximum thickness. The C++ test
(tests/naca/test_naca65.cpp) compares BladeGen's NACA 65 sections with it.

    python tools/reference/make_turbolenza_naca65_reference.py <path-to-turbolenza>
"""
import sys
from pathlib import Path

import numpy as np

root = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\ASUS\Documents\turbolenza")
sys.path.insert(0, str(root / "src"))
from turbolenza.b2b.geometry import CascadeSpec, make_cascade  # noqa: E402

CASES = [(1.2, 0.10, "kovach_sandercock"), (1.2, 0.10, "emery_base"), (0.4, 0.06, "kovach_sandercock"),
         (1.8, 0.12, "kovach_sandercock"), (0.0, 0.10, "emery_base")]
out = Path(__file__).resolve().parents[2] / "tests" / "data" / "turbolenza_naca65_reference.csv"
rows, meta = [], []
for k, (cl0, t, var) in enumerate(CASES):
    c = make_cascade(CascadeSpec(blade_type="NACA65", solidity=1.0, stagger_deg=0.0, cl0=cl0, t_max_over_c=t,
                                 profile_variant=var, surface_continuity="C2"))
    n = c.normalized
    meta.append(f"# case {k}: cl0={cl0} t={t} {var} R_LE={n.leading_edge_radius:.10g} R_TE={n.trailing_edge_radius:.10g}"
                f" max_t={c.metadata['max_thickness_over_c']:.10g}")
    for side, curve in (("SS", c.normalized_curves.suction), ("PS", c.normalized_curves.pressure)):
        a, b = curve.domain
        pts = np.asarray(curve.evaluate(np.linspace(a, b, 1201)))
        rows += [f"{k},{cl0},{t},{var},{side},{x:.15g},{y:.15g}" for x, y in pts]
out.write_text("# TURBOlenza frozen C2 NACA 65 sections (make_turbolenza_naca65_reference.py)\n" + "\n".join(meta) +
               "\ncase,cl0,t,variant,side,x,y\n" + "\n".join(rows) + "\n")
print("wrote", out, len(rows), "rows")
