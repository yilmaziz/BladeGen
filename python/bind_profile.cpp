// Bindings of the profile layer: ProfileParams, Profile2D, quality checks.
#include "bindings_common.hpp"

#include "profile/Profile2D.hpp"
#include "profile/ProfileParams.hpp"
#include "profile/Quality.hpp"

namespace bladegen_py {

using namespace BladeGen::Profile;

namespace {

py::dict sample_dict(const CurveSample& c)
{
    py::dict d;
    d["xy"] = to_np<2>(c.xy);
    d["t"] = to_np(c.t);
    d["s"] = to_np(c.s);
    d["kappa"] = to_np(c.kappa);
    d["piece"] = to_np(c.piece);
    return d;
}

py::dict issue_dict(const Issue& i)
{
    py::dict d;
    d["severity"] = to_string(i.severity);
    d["code"] = i.code;
    d["message"] = i.message;
    d["where"] = vec_np(i.where);
    d["side"] = i.side;
    d["s"] = i.s;
    return d;
}

}  // namespace

void bind_profile(py::module_& m)
{
    py::enum_<ConnectorMode>(m, "ConnectorMode", "How the LE / TE connectors meet the two sides")
        .value("Fair", ConnectorMode::Fair)
        .value("C2Parametric", ConnectorMode::C2Parametric)
        .value("G2Legacy", ConnectorMode::G2Legacy);

    py::class_<EdgeParams>(m, "EdgeParams", "Leading- or trailing-edge family")
        .def(py::init<>())
        .def_readwrite("radius", &EdgeParams::radius)
        .def_readwrite("stretch_ab", &EdgeParams::stretch_ab)
        .def_readwrite("asym", &EdgeParams::asym)
        .def_readwrite("beta_deg", &EdgeParams::beta_deg)
        .def_readwrite("d_beta_ss_deg", &EdgeParams::d_beta_ss_deg)
        .def_readwrite("d_beta_ps_deg", &EdgeParams::d_beta_ps_deg)
        .def_readwrite("tension", &EdgeParams::tension);
    py::class_<SuctionSideParams>(m, "SuctionSideParams")
        .def(py::init<>())
        .def_readwrite("handle_le_frac", &SuctionSideParams::handle_le_frac)
        .def_readwrite("ctrl2_mprime", &SuctionSideParams::ctrl2_mprime)
        .def_readwrite("ctrl2_theta", &SuctionSideParams::ctrl2_theta)
        .def_readwrite("handle_te_frac", &SuctionSideParams::handle_te_frac);
    py::class_<ThicknessParams>(m, "ThicknessParams")
        .def(py::init<>())
        .def_readwrite("t_max", &ThicknessParams::t_max)
        .def_readwrite("mprime_max", &ThicknessParams::mprime_max)
        .def_readwrite("t_le", &ThicknessParams::t_le)
        .def_readwrite("t_te", &ThicknessParams::t_te)
        .def_readwrite("shape_le", &ThicknessParams::shape_le)
        .def_readwrite("shape_te", &ThicknessParams::shape_te)
        .def_readwrite("enforce_c1", &ThicknessParams::enforce_c1);
    py::class_<PressureSideParams>(m, "PressureSideParams")
        .def(py::init<>())
        .def_readwrite("handle_le_frac", &PressureSideParams::handle_le_frac)
        .def_readwrite("handle_te_frac", &PressureSideParams::handle_te_frac)
        .def_readwrite("mprime_ctrl", &PressureSideParams::mprime_ctrl);

    py::class_<ProfileParams>(m, "ProfileParams", "All parameters of one blade profile (unstaggered, chord = 1)")
        .def(py::init<>())
        .def_readwrite("le", &ProfileParams::le)
        .def_readwrite("te", &ProfileParams::te)
        .def_readwrite("ss", &ProfileParams::ss)
        .def_readwrite("thickness", &ProfileParams::thickness)
        .def_readwrite("ps", &ProfileParams::ps)
        .def_readwrite("connector", &ProfileParams::connector)
        .def_static("naca65_1210", &ProfileParams::naca65_1210)
        .def_static("naca65_1210_blade2d", &ProfileParams::naca65_1210_blade2d)
        .def_static("cda_start", &ProfileParams::cda_start)
        .def("get", [](const ProfileParams& p, const std::string& n) { return get_param(p, n); }, py::arg("name"))
        .def("set", [](ProfileParams& p, const std::string& n, double v) { set_param(p, n, v); },
             py::arg("name"), py::arg("value"))
        .def("to_dict", [](const ProfileParams& p) {
            py::dict d;
            for (const auto& i : param_table()) d[py::str(i.name)] = get_param(p, i.name);
            d["connector"] = to_string(p.connector);
            d["ps.mprime_ctrl"] = p.ps.mprime_ctrl;
            d["thickness.enforce_c1"] = p.thickness.enforce_c1;
            return d;
        })
        .def_static("from_dict", [](const py::dict& d) {
            ProfileParams p;
            for (auto item : d) {
                const auto k = py::cast<std::string>(item.first);
                if (k == "connector") p.connector = connector_mode_from_string(py::cast<std::string>(item.second));
                else if (k == "ps.mprime_ctrl") p.ps.mprime_ctrl = py::cast<std::vector<double>>(item.second);
                else if (k == "thickness.enforce_c1") p.thickness.enforce_c1 = py::cast<bool>(item.second);
                else set_param(p, k, py::cast<double>(item.second));
            }
            return p;
        })
        .def("copy", [](const ProfileParams& p) { return p; })
        .def("__copy__", [](const ProfileParams& p) { return p; })
        .def("__deepcopy__", [](const ProfileParams& p, const py::dict&) { return p; });

    m.def("param_table", [] {
        py::list out;
        for (const auto& i : param_table()) {
            py::dict d;
            d["name"] = i.name;
            d["family"] = i.family;
            d["label"] = i.label;
            d["lo"] = i.lo;
            d["hi"] = i.hi;
            out.append(d);
        }
        return out;
    }, "Scalar profile parameters with their family, label and suggested range");
    m.def("blend", &blend, py::arg("a"), py::arg("b"), py::arg("w"));
    m.def("thickness_law", &thickness_law, py::arg("xi"), py::arg("params"), py::arg("xi_max"));
    m.def("validate", &validate, py::arg("params"));

    py::enum_<Side>(m, "Side").value("Suction", Side::Suction).value("Pressure", Side::Pressure);

    py::class_<Profile2D>(m, "Profile2D", "Closed C2 profile in the unstaggered (m', theta) plane")
        .def(py::init<ProfileParams>(), py::arg("params"))
        .def_property_readonly("params", &Profile2D::params)
        .def_property_readonly("period", &Profile2D::period)
        .def_property_readonly("breaks", [](const Profile2D& p) {
            const auto& b = p.breaks();
            return std::vector<double>(b.begin(), b.end());
        })
        .def_property_readonly("t_le_point", &Profile2D::t_le_point)
        .def_property_readonly("t_te_point", &Profile2D::t_te_point)
        .def("characteristic_points", [](const Profile2D& p) {
            py::dict d;
            d["P_LE"] = vec_np(p.p_le());
            d["P_TE"] = vec_np(p.p_te());
            d["P_LE_SS"] = vec_np(p.p_le_ss());
            d["P_LE_PS"] = vec_np(p.p_le_ps());
            d["P_TE_SS"] = vec_np(p.p_te_ss());
            d["P_TE_PS"] = vec_np(p.p_te_ps());
            for (const auto* g : {&p.le_geometry(), &p.te_geometry()}) {
                const std::string k = (g == &p.le_geometry()) ? "LE" : "TE";
                d[py::str("asym_" + k)] = vec_np(g->p_asym);
                d[py::str("mid_" + k)] = vec_np(g->p_mid);
                d[py::str("V_" + k)] = vec_np(g->v_axis);
            }
            return d;
        })
        .def("control_polygons", [](const Profile2D& p) {
            py::dict d;
            d["ss"] = to_np<2>(p.suction_side().ctrl());
            d["ps"] = to_np<2>(p.pressure_side().ctrl());
            d["le"] = to_np<2>(p.le_connector().ctrl());
            d["te"] = to_np<2>(p.te_connector().ctrl());
            d["ps_free"] = to_np<2>(p.ps_free_points());
            return d;
        })
        .def("pieces", [](const Profile2D& p, int n) {
            // every piece in its own parameter, for plotting with one colour each
            py::dict d;
            auto smp = [n](const BladeGen::Math::BSpline2& c) {
                std::vector<Vec2> out;
                for (int i = 0; i < n; ++i) out.push_back(c.eval(static_cast<double>(i) / (n - 1)));
                return to_np<2>(out);
            };
            d["ss"] = smp(p.suction_side());
            d["ps"] = smp(p.pressure_side());
            d["le"] = smp(p.le_connector());
            d["te"] = smp(p.te_connector());
            return d;
        }, py::arg("n") = 200)
        .def("point", &Profile2D::point, py::arg("t"))
        .def("curvature", &Profile2D::curvature, py::arg("t"), py::arg("from_left") = false)
        .def("sample_side", [](const Profile2D& p, Side s, int ne, int nb) { return sample_dict(p.sample_side(s, ne, nb)); },
             py::arg("side"), py::arg("n_edge") = 60, py::arg("n_body") = 240)
        .def("sample_loop", [](const Profile2D& p, int ne, int nb) { return sample_dict(p.sample_loop(ne, nb)); },
             py::arg("n_edge") = 40, py::arg("n_body") = 160)
        .def("thickness", [](const Profile2D& p, int n) {
            const auto t = p.thickness(n);
            py::dict d;
            d["xi"] = to_np(t.xi);
            d["mprime"] = to_np(t.mprime);
            d["measured"] = to_np(t.measured);
            d["law"] = to_np(t.law);
            d["xi_max"] = p.xi_max();
            d["xi_ctrl"] = to_np(p.xi_ctrl());
            return d;
        }, py::arg("n") = 200)
        .def("junctions", [](const Profile2D& p) {
            py::list out;
            for (const auto& j : p.junctions()) {
                py::dict d;
                d["name"] = j.name;
                d["t"] = j.t;
                d["gap"] = j.gap;
                d["tangent_angle"] = j.tangent_angle;
                d["curvature_jump"] = j.curvature_jump;
                d["d1_jump"] = j.d1_jump;
                d["d2_jump"] = j.d2_jump;
                out.append(d);
            }
            return out;
        });

    py::class_<QualityOptions>(m, "QualityOptions")
        .def(py::init<>())
        .def_readwrite("wiggle_rel", &QualityOptions::wiggle_rel)
        .def_readwrite("max_extrema_body", &QualityOptions::max_extrema_body)
        .def_readwrite("n_edge", &QualityOptions::n_edge)
        .def_readwrite("n_body", &QualityOptions::n_body);

    m.def("check_quality", [](const Profile2D& p, const QualityOptions& o) {
        const auto r = check_quality(p, o);
        py::dict d;
        py::list issues;
        for (const auto& i : r.issues) issues.append(issue_dict(i));
        d["issues"] = issues;
        d["ss"] = sample_dict(r.ss);
        d["ps"] = sample_dict(r.ps);
        d["ok"] = r.ok();
        d["errors"] = r.count(Severity::Error);
        d["warnings"] = r.count(Severity::Warning);
        d["summary"] = r.summary();
        d["max_curvature"] = r.max_curvature;
        d["min_thickness"] = r.min_thickness;
        d["self_intersections"] = r.self_intersections;
        d["ss_inflections"] = r.ss_inflections;
        d["ps_inflections"] = r.ps_inflections;
        return d;
    }, py::arg("profile"), py::arg("options") = QualityOptions{},
       "Curvature along both sides (LE point -> TE point) and every geometry check; never throws on a bad shape");
}

}  // namespace bladegen_py
