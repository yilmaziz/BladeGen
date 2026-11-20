"""Tests of the profile bindings: every value comes from the C++ library."""
import math

import numpy as np
import pytest

import bladegen as bg


def test_defaults_and_param_table():
    p = bg.ProfileParams.naca65_1210()
    assert p.connector == bg.ConnectorMode.Fair
    table = bg.param_table()
    names = [t["name"] for t in table]
    assert "le.radius" in names and "thickness.t_max" in names
    assert {t["family"] for t in table} >= {"Leading edge", "Trailing edge", "Suction side",
                                            "Thickness distribution", "Pressure side"}
    d = p.to_dict()
    q = bg.ProfileParams.from_dict(d)
    assert q.get("le.radius") == pytest.approx(p.le.radius)
    p.set("le.radius", 0.02)
    assert p.le.radius == 0.02


def test_profile_passes_through_edge_points_and_is_g2():
    prof = bg.Profile2D(bg.ProfileParams.naca65_1210())
    cp = prof.characteristic_points()
    assert np.allclose(cp["P_LE"], [0, 0]) and np.allclose(cp["P_TE"], [1, 0])
    ss = prof.sample_side(bg.Side.Suction, 40, 120)
    assert np.allclose(ss["xy"][0], [0, 0], atol=1e-13)
    assert np.allclose(ss["xy"][-1], [1, 0], atol=1e-13)
    assert ss["xy"].shape == (40 + 119 + 39, 2)
    for j in prof.junctions():
        assert j["gap"] < 1e-12 and j["tangent_angle"] < 1e-9


def test_quality_flags_but_still_builds():
    prm = bg.ProfileParams.naca65_1210()
    prm.le.tension = 3.0
    prof = bg.Profile2D(prm)  # bad shape: no exception
    q = bg.check_quality(prof)
    assert not q["ok"] and q["errors"] >= 1
    assert all(i["severity"] in ("ERROR", "WARNING", "INFO") for i in q["issues"])
    assert q["ss"]["kappa"].shape == (q["ss"]["xy"].shape[0],)


def test_invalid_input_raises_value_error():
    prm = bg.ProfileParams.naca65_1210()
    prm.le.radius = -1.0
    with pytest.raises(ValueError):
        bg.Profile2D(prm)
    with pytest.raises(ValueError):
        prm.set("no.such", 1.0)


def test_thickness_law_closed_form():
    t = bg.ThicknessParams()
    t.t_max, t.t_le, t.t_te = 0.1, 0.02, 0.004
    assert bg.thickness_law(0.0, t, 0.4) == pytest.approx(0.02)
    assert bg.thickness_law(0.4, t, 0.4) == pytest.approx(0.1)
    assert bg.thickness_law(1.0, t, 0.4) == pytest.approx(0.004)


def test_blade_sections_and_json_roundtrip(tmp_path):
    d = bg.BladeDefinition.demo_rotor()
    b = bg.Blade(d)
    s = b.section(0.5, 20, 60)
    assert s["xyz"].shape == (4 * 19 + 2 * 59 + 1, 3)
    assert np.allclose(s["xyz"][0], s["xyz"][-1])
    p = tmp_path / "rotor.json"
    bg.write_blade_json(d, str(p))
    d2 = bg.read_blade_json(str(p))
    s2 = bg.Blade(d2).section(0.5, 20, 60)
    assert np.allclose(s["xyz"], s2["xyz"], atol=1e-14)


def test_fit_recovers_its_own_profile():
    prm = bg.ProfileParams.naca65_1210()
    prof = bg.Profile2D(prm)
    ss = prof.sample_side(bg.Side.Suction, 30, 80)["xy"]
    ps = prof.sample_side(bg.Side.Pressure, 30, 80)["xy"]
    start = prm.copy()
    start.thickness.t_max *= 1.1
    opt = bg.FitOptions()
    opt.params = [bg.FitParam("thickness.t_max", start.thickness.t_max, 0.05, 0.2)]
    r = bg.fit_profile(ss, ps, start, opt)
    assert r["rms"] < 1e-6
    assert r["params"].thickness.t_max == pytest.approx(prm.thickness.t_max, rel=1e-4)


def test_occt_skin_and_solid(tmp_path):
    b = bg.Blade(bg.BladeDefinition.demo_rotor())
    sh = bg.build_blade(b, 0.0, 16, 48)
    assert sh.blade_valid and sh.blade.n_solids == 1
    assert sh.skin_error < 1e-9
    surf = sh.surface
    assert surf.u_periodic and surf.u_degree == 3
    gap, d1, d2 = surf.seam_jumps(0.5)
    assert gap < 1e-9 and d1 < 1e-8 and d2 < 1e-6  # C2 across the periodic seam
    v, f, fid = sh.blade.tessellate(5e-4)
    assert v.shape[1] == 3 and f.shape[1] == 3 and len(fid) == len(f)
    out = tmp_path / "blade.step"
    sh.blade.write(str(out))
    assert bg.read_step(str(out)).volume == pytest.approx(sh.blade_volume, rel=1e-6)


def test_rotor67_reader():
    import os
    path = os.path.join(os.environ.get("BLADEGEN_SOURCE_DIR", str(bg.SOURCE_DIR)), "data", "rotor67",
                        "Rotor_14_Sections_and_Flowpath_EN.txt")
    d = bg.read_rotor67(path)
    assert len(d.sections) == 14 and d.n_blades == 22
    loop = d.sections[0].profile_conformal()
    assert loop.shape == (161 * 2 + 31 * 2 + 1, 2)
