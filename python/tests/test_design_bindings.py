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
    cfg.casing = bg.WallPoints(-15.0, 80.0, 15.0, 80.0)  # the default casing is not a cylinder
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
    assert not r["layers"][4]["cylindrical"]      # the default casing falls 106 -> 102
    assert not r["layers"][0]["cylindrical"]      # the hub rises 40 -> 43
    d = json.loads(bg.design_defaults())
    assert d["n_blades"] == 18
    pts = lambda w: [w[k] for k in ("z_le", "r_le", "z_te", "r_te")]
    assert pts(d["meridional"]["casing"]) == [-11.0, 106.0, 10.0, 102.0]
    assert pts(d["meridional"]["hub"]) == [-15.0, 40.0, 15.0, 43.0]
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


def test_rotor67_default_layers():
    lay = bg.rotor67_layer(0.0)
    assert lay["stagger_deg"] == pytest.approx(11.6406, abs=1e-3)
    assert bg.rotor67_layer(1.0)["stagger_deg"] > bg.rotor67_layer(0.5)["stagger_deg"]
    d = json.loads(bg.design_defaults())
    assert all(layer["method"] == "Custom" for layer in d["layers"])
    n = json.loads(bg.design_naca65_defaults())
    assert all(layer["method"] == "NACA 65" for layer in n["layers"])
    with pytest.raises(ValueError):
        bg.rotor67_layer(1.5)


def test_centres_of_gravity():
    r = bg.DesignEvaluator().evaluate(bg.design_defaults())
    assert r["blade_volume"] > 0.0
    x, y, z = r["blade_cg"]
    rad = math.hypot(y, z)
    assert r["layers"][0]["r_le"] < rad < r["layers"][-1]["r_le"]
    for lay in r["layers"]:
        assert lay["section_area"] > 0.0
        cx, cy, cz = lay["cg_xyz"]
        assert math.hypot(cy, cz) == pytest.approx(lay["cg_r"], rel=1e-12)
        assert cx == pytest.approx(lay["cg_z"], abs=1e-12)
        assert math.degrees(math.atan2(cy, cz)) == pytest.approx(lay["cg_theta_deg"], abs=1e-9)
    # the section centres of gravity rise from hub to tip
    assert all(a["cg_r"] < b["cg_r"] for a, b in zip(r["layers"], r["layers"][1:]))


def test_restagger_whole_blade_about_the_hub_le():
    ev = bg.DesignEvaluator()
    d = json.loads(bg.design_defaults())
    assert d["restagger"]["angle_deg"] == 0.0 and d["restagger"]["axis"] == "hub"
    a = ev.evaluate(json.dumps(d))
    d["restagger"]["angle_deg"] = -3.0  # opened by 3 deg
    b = ev.evaluate(json.dumps(d))
    assert a["ok"] and b["ok"], b["errors"]
    for la, lb in zip(a["layers"], b["layers"]):
        np.testing.assert_allclose(lb["beta_deg"], la["beta_deg"] - 3.0, atol=1e-10)
        assert lb["stagger_deg"] == pytest.approx(la["stagger_deg"] - 3.0, abs=1e-12)
        assert lb["chord_mprime"] == pytest.approx(la["chord_mprime"], rel=1e-14)
    # the hub LE is the pivot: it stays
    assert b["layers"][0]["mprime_le"] == pytest.approx(0.0, abs=1e-12)
    assert b["layers"][0]["theta_le_deg"] == pytest.approx(a["layers"][0]["theta_le_deg"], abs=1e-12)
    d["restagger"]["axis"] = "shaft"
    with pytest.raises(ValueError):
        ev.evaluate(json.dumps(d))


def test_stacking_on_the_centre_of_gravity_variants_1_and_2():
    ev = bg.DesignEvaluator()
    d = json.loads(bg.design_defaults())
    assert d["stacking"]["mode"] == "edges"
    assert all(lay["chord_scale"] == 1.0 for lay in d["layers"])
    d["stacking"] = {"mode": "reference_point", "variant": "same_z_theta", "reference": "cg", "z_ref": 0.5,
                     "theta_ref_deg": 3.0}
    r1 = ev.evaluate(json.dumps(d))
    assert r1["ok"], r1["errors"]
    for lay in r1["layers"]:  # variant 1: every centre of gravity on the radial line (Z_ref, theta_ref)
        assert lay["cg_z"] == pytest.approx(0.5, abs=1e-10)
        assert lay["cg_theta_deg"] == pytest.approx(3.0, abs=1e-10)
        assert lay["stack_z"] == pytest.approx(lay["cg_z"], abs=1e-12)
    np.testing.assert_allclose(r1["stacking_xr"][:, 0], 0.5, atol=1e-10)
    d["stacking"]["variant"] = "closest_point"
    r2 = ev.evaluate(json.dumps(d))
    assert r2["ok"], r2["errors"]
    # variant 2: the first layer is the same, the others follow the stream-surface normals
    assert r2["layers"][0]["cg_z"] == pytest.approx(0.5, abs=1e-10)
    assert abs(r2["layers"][-1]["cg_z"] - 0.5) > 0.1
    d["stacking"]["variant"] = "radial"
    with pytest.raises(ValueError):
        ev.evaluate(json.dumps(d))


def test_hub_reference_point_is_idempotent_with_a_restagger():
    # pressing "Take Z_ref, theta_ref from the current hub section" again must not move the blade
    ev = bg.DesignEvaluator()
    d = json.loads(bg.design_defaults())
    d["restagger"]["angle_deg"] = 4.0
    z, th = ev.hub_reference_point(json.dumps(d), "cg")
    for variant in ("same_z_theta", "closest_point"):
        d["stacking"] = {"mode": "reference_point", "variant": variant, "reference": "cg", "z_ref": z,
                         "theta_ref_deg": th}
        for _ in range(5):
            z2, th2 = ev.hub_reference_point(json.dumps(d), "cg")
            assert z2 == pytest.approx(z, abs=1e-10) and th2 == pytest.approx(th, abs=1e-10)
            d["stacking"]["z_ref"], d["stacking"]["theta_ref_deg"] = z2, th2
