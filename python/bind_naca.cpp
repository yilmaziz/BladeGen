// Bindings of the NACA 65 layer: parameters, section, metadata, quality.
#include "bindings_common.hpp"

#include "naca/Naca65.hpp"

namespace bladegen_py {

using namespace BladeGen::Naca;
namespace Profile = BladeGen::Profile;

namespace {

py::dict metadata_dict(const Naca65Metadata& m)
{
    py::dict d;
    d["designation"] = m.designation;
    d["le_radius"] = m.le_radius;
    d["le_center"] = vec_np(m.le_center);
    d["te_radius"] = m.te_radius;
    d["te_center"] = vec_np(m.te_center);
    d["camber_angle_deg"] = BladeGen::Math::rad2deg(m.camber_angle);
    d["chi_le_deg"] = BladeGen::Math::rad2deg(m.chi_le);
    d["chi_te_deg"] = BladeGen::Math::rad2deg(m.chi_te);
    d["max_thickness"] = m.max_thickness;
    d["max_thickness_station"] = m.max_thickness_station;
    d["le_fairness"] = py::make_tuple(m.le_fairness[0], m.le_fairness[1]);
    d["le_anchor_x"] = py::make_tuple(m.le_anchor_x[0], m.le_anchor_x[1]);
    d["le_anchor_angle_deg"] = py::make_tuple(m.le_anchor_angle_deg[0], m.le_anchor_angle_deg[1]);
    d["te_start_x"] = py::make_tuple(m.te_start_x[0], m.te_start_x[1]);
    d["te_circle_deviation"] = m.te_circle_deviation;
    d["camber_length"] = m.camber_length;
    return d;
}

}  // namespace

void bind_naca(py::module_& m)
{
    py::enum_<Variant>(m, "Naca65Variant")
        .value("KovachSandercock", Variant::KovachSandercock)
        .value("EmeryBase", Variant::EmeryBase);
    m.def("naca65_variant", &variant_from_string, py::arg("name"), "'kovach_sandercock' or 'emery_base'");

    py::class_<Naca65Params>(m, "Naca65Params")
        .def(py::init<>())
        .def(py::init([](double cl0, double t, const std::string& v) { return Naca65Params{cl0, t, variant_from_string(v)}; }),
             py::arg("cl0") = 1.2, py::arg("t_max_over_c") = 0.10, py::arg("variant") = "kovach_sandercock")
        .def_readwrite("cl0", &Naca65Params::cl0)
        .def_readwrite("t_max_over_c", &Naca65Params::t_max_over_c)
        .def_property("variant", [](const Naca65Params& p) { return to_string(p.variant); },
                      [](Naca65Params& p, const std::string& v) { p.variant = variant_from_string(v); })
        .def("validate", &Naca65Params::validate)
        .def("designation", &Naca65Params::designation)
        .def("__repr__", [](const Naca65Params& p) { return "<Naca65Params " + p.designation() + " " + to_string(p.variant) + ">"; });

    py::class_<Naca65Section>(m, "Naca65Section")
        .def(py::init([](const Naca65Params& p) {
            py::gil_scoped_release nogil;
            return Naca65Section(p);
        }), py::arg("params"))
        .def_property_readonly("params", &Naca65Section::params)
        .def_property_readonly("metadata", [](const Naca65Section& s) { return metadata_dict(s.metadata()); })
        .def_property_readonly("sharp_te", &Naca65Section::sharp_te)
        .def("station", &Naca65Section::station, py::arg("x"))
        .def("half_thickness", &Naca65Section::half_thickness, py::arg("station"))
        .def("camber_point", [](const Naca65Section& s, double x) { return vec_np(s.camber_point(x)); }, py::arg("x"))
        .def("camber_line", [](const Naca65Section& s, int n) {
            std::vector<Vec2> pts;
            for (int i = 0; i < n; ++i) pts.push_back(s.camber_point(0.5 * (1.0 - std::cos(3.14159265358979323846 * i / (n - 1)))));
            return to_np<2>(pts);
        }, py::arg("n") = 200)
        .def("sample_side", [](const Naca65Section& s, Profile::Side side, int ne, int nb) {
            return sample_dict(s.sample_side(side, ne, nb));
        }, py::arg("side"), py::arg("n_edge") = 60, py::arg("n_body") = 240)
        .def("sample_loop", [](const Naca65Section& s, int ne, int nb) { return sample_dict(s.sample_loop(ne, nb)); },
             py::arg("n_edge") = 40, py::arg("n_body") = 160);

    m.def("check_quality", [](const Naca65Section& s, const Profile::QualityOptions& o) {
        return quality_dict(check_quality(s, o));
    }, py::arg("section"), py::arg("options") = Profile::QualityOptions{},
       "Curvature of a NACA 65 section and every geometry check (sharp Emery TE reported as information)");
}

}  // namespace bladegen_py
