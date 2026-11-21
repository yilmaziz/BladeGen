"""Tests of the NACA 65 bindings, NACA 65 + Custom, and the designer GUI (bladegen.gui)."""
import json
import math
import os

import numpy as np
import pytest

import bladegen as bg
from bladegen import gui


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


def test_station_json_round_trip_and_methods():
    st = gui.default_stations()
    assert [s.name for s in st] == ["Hub", "Mid", "Tip"]
    s = st[1]
    s.method = "NACA 65 + Custom"
    s.fitted_to = s.naca_key()
    back = gui.Station.from_json(json.loads(json.dumps(s.to_json())))
    assert back.to_json() == s.to_json()
    bad = s.to_json()
    bad["method"] = "DCA"
    with pytest.raises(ValueError):
        gui.Station.from_json(bad)


@pytest.mark.parametrize("method", gui.METHODS)
def test_evaluate_each_method(method):
    s = gui.default_stations()[1]
    s.method = method
    r = gui.evaluate(s)
    assert r["q"]["ok"], r["q"]["summary"]
    if method == "NACA 65 + Custom":
        assert r.get("unfitted")  # nothing fitted yet: the NACA section is shown
        s.fitted_to = s.naca_key()
        r = gui.evaluate(s)
        assert r["deviation"]["rms"] > 0  # custom (default parameters) vs the NACA target
    if method == "Custom":
        assert "control" in r


def test_evaluate_reports_s_shaped_edge_as_error():
    s = gui.default_stations()[2]
    s.method = "Custom"
    s.custom["le.tension"] = 4.5
    q = gui.evaluate(s)["q"]
    assert not q["ok"] and q["errors"] >= 1


def _tk_root():
    import tkinter as tk

    try:
        return tk.Tk()
    except tk.TclError:
        pytest.skip("no display for Tk")


def test_designer_window_runs_and_exports(tmp_path):
    if os.environ.get("BLADEGEN_SKIP_GUI"):
        pytest.skip("BLADEGEN_SKIP_GUI set")
    import time

    root = _tk_root()
    try:
        app = gui.DesignerApp(root)

        def settle():
            t = time.time()
            while time.time() - t < 120:
                root.update()
                time.sleep(0.02)
                if app.idle():
                    return
            raise AssertionError("GUI did not settle")

        settle()
        assert app.last is not None and app.last["q"]["ok"]
        app.station_var.set(2)
        app._switch_station()
        settle()
        app.method_var.set("Custom")
        app._method_changed()
        app.st.custom["le.tension"] = 4.5
        app._load_custom()
        app._request()
        settle()
        assert app.last["q"]["errors"] >= 1
        specs = app.turbolenza_specs()
        assert [s["station"] for s in specs] == ["Hub", "Mid"]
        assert set(specs[0]) >= {"blade_type", "solidity", "stagger_deg", "cl0", "t_max_over_c", "profile_variant",
                                 "surface_continuity"} and specs[0]["blade_type"] == "NACA65"
        d = app.design_json()
        assert d["format"] == gui.DESIGN_FORMAT and len(d["stations"]) == 3
        app.fig.savefig(tmp_path / "gui.png")
    finally:
        root.destroy()
