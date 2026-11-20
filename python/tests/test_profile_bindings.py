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
