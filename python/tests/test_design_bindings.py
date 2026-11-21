"""Python side of the v1.0.2 designer core: meridional flow path, five layers, mesh, fields."""
import json
import math

import numpy as np
import pytest

import bladegen as bg


def test_bell_mehta_contour():
    assert bg.bell_mehta(0.0) == 0.0
    assert bg.bell_mehta(0.5) == pytest.approx(0.5)
    assert bg.bell_mehta(1.0) == pytest.approx(1.0)


def test_flowpath_cylindrical_wrap_and_le_line():
    cfg = bg.MeridionalConfig()
    cfg.hub = bg.WallPoints(-15.0, 50.0, 15.0, 50.0)
    cfg.hub_law = bg.WallLaw.Linear
    fp = bg.FlowPath(cfg)
    s = fp.surface(0.0)
    assert s.cylindrical
    assert s.mprime_te == pytest.approx(30.0 / 50.0)
    assert s.xr(0.1)[0] == pytest.approx(-15.0 + 50.0 * 0.1)
    np.testing.assert_allclose(fp.le_point(0.5), [-15.0, 65.0])
    bad = bg.MeridionalConfig()
    bad.casing = bg.WallPoints(-15.0, 30.0, 15.0, 80.0)
    with pytest.raises(ValueError):
        bg.FlowPath(bad)


def test_default_design_evaluates():
    ev = bg.DesignEvaluator()
    r = ev.evaluate(bg.design_defaults())
    assert r["ok"], r["errors"]
    assert len(r["layers"]) == 5
    nj, ni, three = r["blade"].shape
    assert three == 3 and nj == 25 and ni == len(r["layers"][0]["xyz"])
    for lay in r["layers"]:
        assert lay["mfrac"][0] == pytest.approx(0.0, abs=1e-12)
        assert lay["mfrac"][-1] == pytest.approx(1.0, abs=1e-9)
        assert lay["beta_le"] > lay["beta_te"]
    assert r["layers"][4]["cylindrical"]          # the default casing is a cylinder (R = 80)
    assert not r["layers"][0]["cylindrical"]      # the hub rises 40 -> 60
    f = r["beta_field"]
    assert f["value"].shape == f["xr"].shape[:2]


def test_design_json_v1_import_and_errors():
    v1 = {"format": "bladegen-design-1", "stations": [
        {"name": n, "span": s, "method": "NACA 65", "cl0": c, "t_max_over_c": 0.1, "variant": "kovach_sandercock",
         "stagger_deg": g, "solidity": 1.2} for n, s, c, g in (("Hub", 0, 1.5, 25), ("Mid", 0.5, 1.2, 37), ("Tip", 1, 0.8, 48))]}
    d = json.loads(bg.design_normalize(json.dumps(v1)))
    assert d["format"] == "bladegen-design-2"
    assert [lay["span"] for lay in d["layers"]] == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert d["layers"][1]["naca"]["cl0"] == pytest.approx(1.35)
    bad = json.loads(bg.design_defaults())
    bad["layers"][0]["method"] = "DCA"
    with pytest.raises(ValueError):
        bg.design_normalize(json.dumps(bad))
    bad = json.loads(bg.design_defaults())
    bad["meridional"]["hub"]["z_te"] = -20.0
    with pytest.raises(ValueError):
        bg.DesignEvaluator().evaluate(json.dumps(bad))
