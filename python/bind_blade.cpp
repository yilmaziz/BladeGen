// Bindings of the blade layer (meridional channel, placement, 3-D sections) and JSON I/O.
#include "bindings_common.hpp"

#include "blade/Blade.hpp"
#include "blade/Meridional.hpp"
#include "io/BladeJson.hpp"

namespace bladegen_py {

using namespace BladeGen::Blade;
namespace IO = BladeGen::IO;

namespace {

py::dict section_dict(const Section3D& s)
{
    py::dict d;
    d["span"] = s.span;
    d["params"] = s.params;
    d["placement"] = s.placement;
    d["mprime_le"] = s.mprime_le;
    d["chord_mprime"] = s.chord_mprime;
    d["unstaggered"] = to_np<2>(s.unstaggered);
    d["conformal"] = to_np<2>(s.conformal);
    d["xyz"] = to_np<3>(s.xyz);
    d["piece"] = to_np(s.piece);
    std::vector<Vec2> merid;
    for (const auto& q : s.xyz) merid.emplace_back(q.x(), std::hypot(q.y(), q.z()));
    d["xr"] = to_np<2>(merid);
    return d;
}

Vec2 vec2_of(const InArray& a) { return from_np<2>(a.attr("reshape")(1, 2), "point")[0]; }

}  // namespace

void bind_blade(py::module_& m)
{
    py::class_<MeridionalCurve>(m, "MeridionalCurve", "Stream line (x, r) with arc length m and conformal m' = int dm / r")
        .def_static("through_points", [](const InArray& a) { return MeridionalCurve::through_points(from_np<2>(a, "points")); },
                    py::arg("xr"))
        .def_static("line", [](const InArray& a, const InArray& b) { return MeridionalCurve::line(vec2_of(a), vec2_of(b)); },
                    py::arg("a"), py::arg("b"))
        .def_property_readonly("length", &MeridionalCurve::length)
        .def_property_readonly("mprime_length", &MeridionalCurve::mprime_length)
        .def("at_m", [](const MeridionalCurve& c, double mm) { return vec_np(c.at_m(mm)); }, py::arg("m"))
        .def("mprime_of_m", &MeridionalCurve::mprime_of_m, py::arg("m"))
        .def("m_of_mprime", &MeridionalCurve::m_of_mprime, py::arg("mprime"))
        .def("sample", [](const MeridionalCurve& c, int n) { return to_np<2>(c.sample(n)); }, py::arg("n") = 101)
        .def("control_points", [](const MeridionalCurve& c) { return to_np<2>(c.curve().ctrl()); });

    py::class_<MeridionalChannel>(m, "MeridionalChannel", "Hub and shroud; streamline(span) interpolates them")
        .def(py::init<MeridionalCurve, MeridionalCurve, int>(), py::arg("hub"), py::arg("shroud"), py::arg("n_samples") = 161)
        .def_static("straight", &MeridionalChannel::straight, py::arg("x0"), py::arg("x1"), py::arg("r_hub0"),
                    py::arg("r_hub1"), py::arg("r_tip0"), py::arg("r_tip1"))
        .def_property_readonly("hub", &MeridionalChannel::hub)
        .def_property_readonly("shroud", &MeridionalChannel::shroud)
        .def("streamline", &MeridionalChannel::streamline, py::arg("span"));

    py::class_<SectionPlacement>(m, "SectionPlacement", "Where a profile sits on its stream surface")
        .def(py::init([](double m_le, double m_te, double stagger, double theta, bool mirror) {
                 return SectionPlacement{m_le, m_te, stagger, theta, mirror};
             }),
             py::arg("m_le") = 0.25, py::arg("m_te") = 0.75, py::arg("stagger_deg") = 30.0, py::arg("theta_le_deg") = 0.0,
             py::arg("mirror") = false)
        .def_readwrite("m_le", &SectionPlacement::m_le)
        .def_readwrite("m_te", &SectionPlacement::m_te)
        .def_readwrite("stagger_deg", &SectionPlacement::stagger_deg)
        .def_readwrite("theta_le_deg", &SectionPlacement::theta_le_deg)
        .def_readwrite("mirror", &SectionPlacement::mirror)
        .def("copy", [](const SectionPlacement& p) { return p; });

    py::class_<BladeStation>(m, "BladeStation")
        .def(py::init([](double span, const BladeGen::Profile::ProfileParams& p, const SectionPlacement& pl,
                         std::optional<MeridionalCurve> sl) { return BladeStation{span, p, pl, std::move(sl)}; }),
             py::arg("span"), py::arg("profile"), py::arg("placement"), py::arg("streamline") = std::nullopt)
        .def_readwrite("span", &BladeStation::span)
        .def_readwrite("profile", &BladeStation::profile)
        .def_readwrite("placement", &BladeStation::placement)
        .def_readwrite("streamline", &BladeStation::streamline);

    py::enum_<SpanInterpolation>(m, "SpanInterpolation")
        .value("Linear", SpanInterpolation::Linear)
        .value("Cubic", SpanInterpolation::Cubic);

    py::class_<BladeDefinition>(m, "BladeDefinition", "Channel + spanwise stations of a blade")
        .def(py::init<>())
        .def_readwrite("name", &BladeDefinition::name)
        .def_readwrite("channel", &BladeDefinition::channel)
        .def_readwrite("stations", &BladeDefinition::stations)
        .def_readwrite("n_blades", &BladeDefinition::n_blades)
        .def_readwrite("interpolation", &BladeDefinition::interpolation)
        .def_readwrite("span_extension", &BladeDefinition::span_extension)
        .def_readwrite("n_sections", &BladeDefinition::n_sections)
        .def_readwrite("tip_clearance", &BladeDefinition::tip_clearance)
        .def("validate", &BladeDefinition::validate)
        .def_static("demo_rotor", &BladeDefinition::demo_rotor)
        .def("to_json", [](const BladeDefinition& d) { return IO::to_json(d).dump(2); })
        .def_static("from_json", [](const std::string& s) { return IO::blade_from_json(nlohmann::json::parse(s)); })
        .def("copy", [](const BladeDefinition& d) { return d; });

    py::class_<Blade>(m, "Blade", "A blade built from its definition")
        .def(py::init<BladeDefinition>(), py::arg("definition"))
        .def_property_readonly("definition", &Blade::definition)
        .def("profile_at", &Blade::profile_at, py::arg("span"))
        .def("placement_at", &Blade::placement_at, py::arg("span"))
        .def("streamline_at", &Blade::streamline_at, py::arg("span"))
        .def("section", [](const Blade& b, double span, int ne, int nb) { return section_dict(b.section(span, ne, nb)); },
             py::arg("span"), py::arg("n_edge") = 40, py::arg("n_body") = 160)
        .def("sections", [](const Blade& b, int ne, int nb) {
            py::list out;
            for (const auto& s : b.sections(ne, nb)) out.append(section_dict(s));
            return out;
        }, py::arg("n_edge") = 40, py::arg("n_body") = 160);

    m.def("build_section", [](const BladeGen::Profile::ProfileParams& p, const SectionPlacement& pl, const MeridionalCurve& sl,
                              double span, int ne, int nb) { return section_dict(build_section(p, pl, sl, span, ne, nb)); },
          py::arg("profile"), py::arg("placement"), py::arg("streamline"), py::arg("span") = 0.0, py::arg("n_edge") = 40,
          py::arg("n_body") = 160);

    m.def("read_blade_json", &IO::read_blade_json, py::arg("path"));
    m.def("write_blade_json", &IO::write_blade_json, py::arg("definition"), py::arg("path"));
    m.def("read_profile_json", &IO::read_profile_json, py::arg("path"),
          "BladeGen profile JSON, or a blade2d '<name>_blade_params.json' (detected by its keys)");
    m.def("write_profile_json", &IO::write_profile_json, py::arg("params"), py::arg("path"));
}

}  // namespace bladegen_py
