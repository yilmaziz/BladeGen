"""Tests of the NACA 65 bindings and NACA 65 + Custom."""
import math

import numpy as np
import pytest

import bladegen as bg


def test_naca65_section_metadata_and_samples():
    p = bg.Naca65Params(1.2, 0.10, "kovach_sandercock")
    assert p.designation() == "NACA 65-(12A10)10"
    s = bg.Naca65Section(p)
    m = s.metadata
    assert m["le_radius"] == pytest.approx(0.00687)
    assert m["te_radius"] == pytest.approx(0.008)
    assert math.tan(math.radians(m["camber_angle_deg"]) / 4) == pytest.approx(0.1103 * 1.2)
    ss = s.sample_side(bg.Side.Suction, 40, 120)
    assert np.allclose(ss["xy"][0], [0, 0], atol=1e-14) and np.allclose(ss["xy"][-1], [1, 0], atol=1e-14)
    assert ss["kappa"][0] == pytest.approx(1 / m["le_radius"], rel=1e-10)
    q = bg.check_quality(s)
    assert q["ok"], q["summary"]
    e = bg.Naca65Section(bg.Naca65Params(1.2, 0.10, "emery_base"))
    assert e.sharp_te and any(i["code"] == "sharp_te" for i in bg.check_quality(e)["issues"])


def test_naca65_invalid_input_raises():
    with pytest.raises(ValueError):
        bg.Naca65Params(1.2, 0.10, "dca")
    with pytest.raises(ValueError):
        bg.Naca65Section(bg.Naca65Params(3.0, 0.10, "emery_base"))
    with pytest.raises(ValueError):
        bg.Naca65Section(bg.Naca65Params(1.0, 0.30, "emery_base"))


def test_naca65_plus_custom_fit():
    s = bg.Naca65Section(bg.Naca65Params(1.2, 0.10, "kovach_sandercock"))
    r = bg.fit_naca65(s)
    assert r["rms"] < 3e-4
    q = bg.check_quality(bg.Profile2D(r["params"]))
    assert q["ok"], q["summary"]



def test_s_shaped_edge_is_an_error():
    d = bg.ProfileParams.naca65_1210().to_dict()
    d["le.tension"] = 4.5
    q = bg.check_quality(bg.Profile2D(bg.ProfileParams.from_dict(d)))
    assert not q["ok"] and q["errors"] >= 1
