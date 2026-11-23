// Bindings of the meridional flow path and the designer's blade (five layers, mesh, fields).
// A design crosses the boundary as its JSON ("bladegen-design-2"), so Python scripts and the
// C++ designer read and write the same files.
#include "bindings_common.hpp"

#include "design/BladeDesign.hpp"
#include "design/DesignSolid.hpp"
#include "io/BladeJson.hpp"

#include <pybind11/stl.h>

namespace bladegen_py {

using namespace BladeGen;
namespace D = BladeGen::Design;

namespace {

py::array_t<double> grid_np(const D::Grid& g)
{
    py::array_t<double> a({static_cast<py::ssize_t>(g.nj), static_cast<py::ssize_t>(g.ni), static_cast<py::ssize_t>(3)});
    auto r = a.mutable_unchecked<3>();
    for (int j = 0; j < g.nj; ++j)
        for (int i = 0; i < g.ni; ++i)
            for (int k = 0; k < 3; ++k) r(j, i, k) = g.at(i, j)[k];
    return a;
}

py::dict field_dict(const D::MeridionalField& f)
{
    py::dict d;
    py::array_t<double> xr({static_cast<py::ssize_t>(f.nj), static_cast<py::ssize_t>(f.ni), static_cast<py::ssize_t>(2)});
    py::array_t<double> v({static_cast<py::ssize_t>(f.nj), static_cast<py::ssize_t>(f.ni)});
    auto rx = xr.mutable_unchecked<3>();
    auto rv = v.mutable_unchecked<2>();
    for (int j = 0; j < f.nj; ++j)
        for (int i = 0; i < f.ni; ++i) {
            const auto k = static_cast<std::size_t>(i + f.ni * j);
            rx(j, i, 0) = f.xr[k].x();
            rx(j, i, 1) = f.xr[k].y();
            rv(j, i) = f.value[k];
        }
    d["xr"] = xr;
    d["value"] = v;
    d["lo"] = f.lo;
    d["hi"] = f.hi;
    return d;
}

py::dict result_dict(const D::DesignResult& r)
{
    py::dict d;
    d["ok"] = r.ok();
    d["errors"] = r.errors;
    d["seconds"] = r.seconds;
    d["blade_volume"] = r.blade_volume;
    d["blade_cg"] = py::make_tuple(r.blade_cg.x(), r.blade_cg.y(), r.blade_cg.z());
    py::list layers;
    for (const auto& L : r.layers) {
        py::dict l;
        l["error"] = L.error;
        if (L.error.empty() && L.section) {
            l["designation"] = L.section->designation;
            l["quality"] = quality_dict(L.section->quality);
            l["loop"] = sample_dict(L.section->loop);
            l["camber_xy"] = to_np<2>(L.section->camber.point);
            l["cylindrical"] = L.surface.cylindrical();
            l["stagger_deg"] = Math::rad2deg(L.gamma);
            l["chord_mprime"] = L.chord_mprime;
            l["mprime_le"] = L.mprime_le;
            l["section_area"] = L.section_area;
            l["cg_xyz"] = py::make_tuple(L.cg_xyz.x(), L.cg_xyz.y(), L.cg_xyz.z());
            l["cg_z"] = L.cg_xr.x();
            l["cg_r"] = L.cg_xr.y();
            l["cg_theta_deg"] = Math::rad2deg(L.cg_theta);
            l["cg_mprime"] = L.cg_mprime;
            l["stack_mprime"] = L.stack_point.x();
            l["stack_theta_deg"] = Math::rad2deg(L.stack_point.y());
            l["stack_z"] = L.stack_xr.x();
            l["stack_r"] = L.stack_xr.y();
            l["theta_le_deg"] = Math::rad2deg(L.theta_le);
            l["chord"] = L.chord;
            l["r_le"] = L.r_le;
            l["r_te"] = L.r_te;
            l["solidity_3d"] = L.solidity_3d;
            l["conformal"] = to_np<2>(L.conformal);
            l["xyz"] = to_np<3>(L.xyz);
            l["mprime"] = to_np(L.mprime);
            l["m"] = to_np(L.m);
            l["mfrac"] = to_np(L.mfrac);
            l["xc"] = to_np(L.xc);
            l["theta_deg"] = to_np(L.theta_deg);
            l["beta_deg"] = to_np(L.beta_deg);
            l["thickness"] = to_np(L.thickness);
            l["camber_xr"] = to_np<2>(L.camber_xr);
            l["beta_le"] = L.beta_le;
            l["beta_te"] = L.beta_te;
        }
        layers.append(l);
    }
    d["layers"] = layers;
    d["hub"] = to_np<2>(r.hub);
    d["casing"] = to_np<2>(r.casing);
    py::list sl;
    for (const auto& s : r.streamlines) sl.append(to_np<2>(s));
    d["streamlines"] = sl;
    d["le_line"] = to_np<2>(r.le_line);
    d["te_line"] = to_np<2>(r.te_line);
    if (!r.blade.p.empty()) {
        d["blade"] = grid_np(r.blade);
        d["hub_wall"] = grid_np(r.hub_wall);
        d["casing_wall"] = grid_np(r.casing_wall);
        d["row_span"] = to_np(r.row_span);
        d["blade_le"] = to_np<2>(r.blade_le);
        d["blade_te"] = to_np<2>(r.blade_te);
        d["stacking_xr"] = to_np<2>(r.stacking_xr);
        d["stacking_xyz"] = to_np<3>(r.stacking_xyz);
        d["beta_field"] = field_dict(r.beta);
        d["thickness_field"] = field_dict(r.thickness);
    }
    return d;
}

}  // namespace

void bind_design(py::module_& m)
{
    py::enum_<Blade::WallLaw>(m, "WallLaw")
        .value("Linear", Blade::WallLaw::Linear)
        .value("BellMehta", Blade::WallLaw::BellMehta);
    m.def("bell_mehta", &Blade::bell_mehta, py::arg("xi"), "Bell & Mehta (1988) fifth-order contour 6ξ⁵ − 15ξ⁴ + 10ξ³");

    py::class_<Blade::WallPoints>(m, "WallPoints")
        .def(py::init([](double z_le, double r_le, double z_te, double r_te) { return Blade::WallPoints{z_le, r_le, z_te, r_te}; }),
             py::arg("z_le"), py::arg("r_le"), py::arg("z_te"), py::arg("r_te"))
        .def_readwrite("z_le", &Blade::WallPoints::z_le)
        .def_readwrite("r_le", &Blade::WallPoints::r_le)
        .def_readwrite("z_te", &Blade::WallPoints::z_te)
        .def_readwrite("r_te", &Blade::WallPoints::r_te)
        .def("__repr__", [](const Blade::WallPoints& w) {
            return "WallPoints(LE (" + std::to_string(w.z_le) + ", " + std::to_string(w.r_le) + "), TE (" + std::to_string(w.z_te) + ", " +
                   std::to_string(w.r_te) + "))";
        });

    py::class_<Blade::MeridionalConfig>(m, "MeridionalConfig",
                                        "Initial meridional configuration: hub and casing LE / TE points (z, R) and their laws")
        .def(py::init<>())
        .def_readwrite("hub", &Blade::MeridionalConfig::hub)
        .def_readwrite("casing", &Blade::MeridionalConfig::casing)
        .def_readwrite("hub_law", &Blade::MeridionalConfig::hub_law)
        .def_readwrite("casing_law", &Blade::MeridionalConfig::casing_law)
        .def_readwrite("inlet_extension", &Blade::MeridionalConfig::inlet_extension)
        .def_readwrite("outlet_extension", &Blade::MeridionalConfig::outlet_extension)
        .def("validate", &Blade::MeridionalConfig::validate);

    py::class_<Blade::StreamSurface>(m, "StreamSurface", "Conformal map (m', θ) -> (x, r, θ) of one layer; m' from its LE")
        .def_property_readonly("span", &Blade::StreamSurface::span)
        .def_property_readonly("cylindrical", &Blade::StreamSurface::cylindrical)
        .def_property_readonly("mprime_te", &Blade::StreamSurface::mprime_te)
        .def_property_readonly("m_te", &Blade::StreamSurface::m_te)
        .def_property_readonly("mprime_inlet", &Blade::StreamSurface::mprime_inlet)
        .def_property_readonly("mprime_outlet", &Blade::StreamSurface::mprime_outlet)
        .def("xr", [](const Blade::StreamSurface& s, double mp) { return vec_np(s.xr(mp)); }, py::arg("mprime"))
        .def("m_of_mprime", &Blade::StreamSurface::m_of_mprime, py::arg("mprime"))
        .def("cartesian", [](const Blade::StreamSurface& s, double mp, double th) { return vec_np(s.cartesian(mp, th)); },
             py::arg("mprime"), py::arg("theta"))
        .def("curvature", &Blade::StreamSurface::curvature, py::arg("mprime"))
        .def("polyline", [](const Blade::StreamSurface& s, int n) { return to_np<2>(s.polyline(n)); }, py::arg("n_per_piece") = 40);

    py::class_<Blade::FlowPath>(m, "FlowPath", "Hub, casing and layer streamlines of a meridional configuration")
        .def(py::init<Blade::MeridionalConfig>(), py::arg("config"))
        .def("surface", &Blade::FlowPath::surface, py::arg("span"))
        .def("le_point", [](const Blade::FlowPath& f, double s) { return vec_np(f.le_point(s)); }, py::arg("span"))
        .def("te_point", [](const Blade::FlowPath& f, double s) { return vec_np(f.te_point(s)); }, py::arg("span"))
        .def("streamline_points", [](const Blade::FlowPath& f, double s, int n) { return to_np<2>(f.surface(s).polyline(n)); },
             py::arg("span"), py::arg("n_per_piece") = 40);

    py::class_<D::Evaluator>(m, "DesignEvaluator", "Evaluates designer blades (cached 2-D sections, parallel layers)")
        .def(py::init<>())
        .def("evaluate", [](D::Evaluator& ev, const std::string& design_json) {
            const auto d = D::design_from_json(nlohmann::json::parse(design_json));
            D::DesignResult r;
            {
                py::gil_scoped_release nogil;
                r = ev.evaluate(d);
            }
            return result_dict(r);
        }, py::arg("design_json"), "Evaluate a design given as its JSON text; returns a dict of numpy arrays")
        .def("hub_reference_point", [](D::Evaluator& ev, const std::string& design_json, const std::string& reference) {
            const auto d = D::design_from_json(nlohmann::json::parse(design_json));
            const Blade::FlowPath fp(d.meridional);
            const auto p = D::hub_reference_point(d, fp, *ev.section(d.layers.front(), d.mesh), D::stacking_reference_from_string(reference));
            return py::make_tuple(p.x(), Math::rad2deg(p.y()));
        }, py::arg("design_json"), py::arg("reference") = "cg",
           "(Z_ref [mm], theta_ref [deg]) of the hub section's reference point before any re-stagger")
        .def("write_step", [](D::Evaluator& ev, const std::string& design_json, const std::string& path) {
            const auto d = D::design_from_json(nlohmann::json::parse(design_json));
            py::gil_scoped_release nogil;
            const auto s = D::build_design_solid(d, ev);
            Occt::write_shape(s.blade, path, Occt::Unit::Millimetre);
            return s.blade_valid;
        }, py::arg("design_json"), py::arg("path"), "OpenCASCADE blade solid (design in mm) to STEP / IGES / BREP / STL");
    m.def("design_defaults", [] { return D::to_json(D::BladeDesign::defaults()).dump(2); },
          "JSON of the default five-layer axial compressor design");
    m.def("design_naca65_defaults", [] { return D::to_json(D::BladeDesign::naca65_defaults()).dump(2); },
          "JSON of the five-layer NACA 65 design");
    m.def("rotor67_layer", [](double span) {
        const auto l = D::rotor67_layer(span);
        py::dict d;
        d["span"] = l.span;
        d["stagger_deg"] = l.stagger_deg;
        d["theta_le_deg"] = l.theta_le_deg;
        d["profile"] = py::module_::import("json").attr("loads")(IO::to_json(l.profile).dump());
        return d;
    }, py::arg("span"), "NASA Rotor 67 custom profile, stagger and LE angle at a span fraction (blended fitted sections)");
    m.def("design_normalize", [](const std::string& j) { return D::to_json(D::design_from_json(nlohmann::json::parse(j))).dump(2); },
          py::arg("design_json"), "Validate a design (v2 or the v1.0.1 format) and return its v2 JSON");
}

}  // namespace bladegen_py
