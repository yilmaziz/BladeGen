// Bindings of the OpenCASCADE layer: NURBS skin, solids, fillet, export, tessellation.
#include "bindings_common.hpp"

#include "occt/Fillet.hpp"
#include "occt/OcctBlade.hpp"

#include <GeomLProp_SLProps.hxx>
#include <Geom_BSplineSurface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <cmath>

namespace bladegen_py {

namespace Occt = BladeGen::Occt;

namespace {

/// A TopoDS_Shape handle for Python.
struct Shape
{
    TopoDS_Shape s;
};

/// A B-spline surface handle for Python.
struct Surface
{
    Handle(Geom_BSplineSurface) s;
};

std::vector<std::vector<Vec3>> loops_from(const py::list& loops)
{
    std::vector<std::vector<Vec3>> out;
    for (const auto& l : loops) out.push_back(from_np<3>(py::cast<InArray>(l), "loop"));
    return out;
}

Occt::Unit unit_of(const std::string& u)
{
    if (u == "mm") return Occt::Unit::Millimetre;
    if (u == "m") return Occt::Unit::Metre;
    throw std::invalid_argument("unit must be 'mm' or 'm'");
}

}  // namespace

void bind_occt(py::module_& m)
{
    py::class_<Shape>(m, "Shape", "An OpenCASCADE shape (solid, shell, face, ...)")
        .def_property_readonly("is_null", [](const Shape& s) { return s.s.IsNull(); })
        .def_property_readonly("volume", [](const Shape& s) { return Occt::volume(s.s); })
        .def_property_readonly("area", [](const Shape& s) { return Occt::area(s.s); })
        .def_property_readonly("is_valid", [](const Shape& s) { return Occt::is_valid(s.s); })
        .def_property_readonly("n_faces", [](const Shape& s) { return Occt::count_faces(s.s); })
        .def_property_readonly("n_solids", [](const Shape& s) { return Occt::count_solids(s.s); })
        .def("write", [](const Shape& s, const std::string& path, const std::string& unit) {
            Occt::write_shape(s.s, path, unit_of(unit));
        }, py::arg("path"), py::arg("unit") = "mm", "STEP / IGES / STL / BREP by extension")
        .def("tessellate", [](const Shape& s, double lin, double ang) {
            const auto mesh = Occt::tessellate(s.s, lin, ang);
            py::array_t<int> tri({static_cast<py::ssize_t>(mesh.triangles.size()), static_cast<py::ssize_t>(3)});
            auto r = tri.mutable_unchecked<2>();
            for (std::size_t i = 0; i < mesh.triangles.size(); ++i)
                for (int k = 0; k < 3; ++k) r(static_cast<py::ssize_t>(i), k) = mesh.triangles[i][static_cast<std::size_t>(k)];
            return py::make_tuple(to_np<3>(mesh.vertices), tri, to_np(mesh.face_of_triangle));
        }, py::arg("linear_deflection") = 2e-4, py::arg("angular_deflection") = 0.3,
           "(vertices (N,3), triangles (M,3), face id per triangle)")
        .def("distances", [](const Shape& s, const InArray& pts) { return to_np(Occt::distances(from_np<3>(pts, "points"), s.s)); },
             py::arg("points"));

    py::class_<Surface>(m, "Surface", "A B-spline surface (the blade or fillet skin)")
        .def_property_readonly("u_degree", [](const Surface& w) { return w.s->UDegree(); })
        .def_property_readonly("v_degree", [](const Surface& w) { return w.s->VDegree(); })
        .def_property_readonly("n_poles", [](const Surface& w) { return py::make_tuple(w.s->NbUPoles(), w.s->NbVPoles()); })
        .def_property_readonly("u_periodic", [](const Surface& w) { return static_cast<bool>(w.s->IsUPeriodic()); })
        .def_property_readonly("continuity", [](const Surface& w) {
            static const char* names[] = {"C0", "G1", "C1", "G2", "C2", "C3", "CN"};
            return std::string(names[static_cast<int>(w.s->Continuity())]);
        })
        .def("poles", [](const Surface& w) {
            py::array_t<double> a({w.s->NbUPoles(), w.s->NbVPoles(), 3});
            auto r = a.mutable_unchecked<3>();
            for (int i = 1; i <= w.s->NbUPoles(); ++i)
                for (int j = 1; j <= w.s->NbVPoles(); ++j) {
                    const gp_Pnt p = w.s->Pole(i, j);
                    r(i - 1, j - 1, 0) = p.X();
                    r(i - 1, j - 1, 1) = p.Y();
                    r(i - 1, j - 1, 2) = p.Z();
                }
            return a;
        })
        .def("grid", [](const Surface& w, int nu, int nv) {
            py::array_t<double> a({nu, nv, 3});
            auto r = a.mutable_unchecked<3>();
            double u1, u2, v1, v2;
            w.s->Bounds(u1, u2, v1, v2);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    const gp_Pnt p = w.s->Value(u1 + (u2 - u1) * i / (nu - 1), v1 + (v2 - v1) * j / (nv - 1));
                    r(i, j, 0) = p.X();
                    r(i, j, 1) = p.Y();
                    r(i, j, 2) = p.Z();
                }
            return a;
        }, py::arg("nu") = 120, py::arg("nv") = 40, "points on a uniform (u, v) grid, shape (nu, nv, 3)")
        .def("point", [](const Surface& w, double u, double v) {
            const gp_Pnt p = w.s->Value(u, v);
            return vec_np(Vec3(p.X(), p.Y(), p.Z()));
        }, py::arg("u"), py::arg("v"))
        .def("curvatures", [](const Surface& w, int nu, int nv) {
            py::array_t<double> K({nu, nv}), H({nu, nv});
            auto k = K.mutable_unchecked<2>();
            auto h = H.mutable_unchecked<2>();
            double u1, u2, v1, v2;
            w.s->Bounds(u1, u2, v1, v2);
            GeomLProp_SLProps pr(w.s, 2, 1e-12);
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j) {
                    pr.SetParameters(u1 + (u2 - u1) * i / (nu - 1), v1 + (v2 - v1) * j / (nv - 1));
                    const bool ok = pr.IsCurvatureDefined();
                    k(i, j) = ok ? pr.GaussianCurvature() : std::nan("");
                    h(i, j) = ok ? pr.MeanCurvature() : std::nan("");
                }
            return py::make_tuple(K, H);
        }, py::arg("nu") = 120, py::arg("nv") = 40, "(Gaussian, mean) curvature on a uniform (u, v) grid")
        .def("seam_jumps", [](const Surface& w, double v) {
            gp_Pnt a, b;
            gp_Vec au, av, auu, avv, auv, bu, bv, buu, bvv, buv;
            double u1, u2, v1, v2;
            w.s->Bounds(u1, u2, v1, v2);
            w.s->D2(u1 + 1e-12, v, a, au, av, auu, avv, auv);
            w.s->D2(u2 - 1e-12, v, b, bu, bv, buu, bvv, buv);
            return py::make_tuple(a.Distance(b), (au - bu).Magnitude() / au.Magnitude(), (auu - buu).Magnitude() / auu.Magnitude());
        }, py::arg("v") = 0.5, "(gap, relative d1 jump, relative d2 jump) across the periodic seam: the C2 check");

    py::class_<Occt::BladeShapes>(m, "BladeShapes", "Everything the OCCT layer builds for one blade")
        .def_property_readonly("blade", [](const Occt::BladeShapes& b) { return Shape{b.blade}; })
        .def_property_readonly("hub", [](const Occt::BladeShapes& b) { return Shape{b.hub}; })
        .def_property_readonly("blade_face", [](const Occt::BladeShapes& b) { return Shape{b.blade_face}; })
        .def_property_readonly("blade_on_hub", [](const Occt::BladeShapes& b) { return Shape{b.blade_on_hub}; })
        .def_property_readonly("blade_with_fillet", [](const Occt::BladeShapes& b) { return Shape{b.blade_with_fillet}; })
        .def_property_readonly("fillet_face", [](const Occt::BladeShapes& b) { return Shape{b.fillet_face}; })
        .def_property_readonly("surface", [](const Occt::BladeShapes& b) { return Surface{b.skin.surface}; })
        .def_property_readonly("skin_error", [](const Occt::BladeShapes& b) { return b.skin.max_error; })
        .def_readonly("fillet_ok", &Occt::BladeShapes::fillet_ok)
        .def_readonly("fillet_message", &Occt::BladeShapes::fillet_message)
        .def_readonly("fillet_g1_deg", &Occt::BladeShapes::fillet_g1_deg)
        .def_readonly("fillet_ball_error", &Occt::BladeShapes::fillet_ball_error)
        .def_readonly("blade_volume", &Occt::BladeShapes::blade_volume)
        .def_readonly("blade_area", &Occt::BladeShapes::blade_area)
        .def_readonly("blade_valid", &Occt::BladeShapes::blade_valid);

    m.def("build_blade", [](const BladeGen::Blade::Blade& b, double fillet, int ne, int nb) {
        py::gil_scoped_release nogil;
        return Occt::build_blade(b, fillet, ne, nb);
    }, py::arg("blade"), py::arg("fillet_radius") = 0.0, py::arg("n_edge") = 40, py::arg("n_body") = 160,
       "C2 NURBS skin, blade solid trimmed on hub and casing, optional rolling-ball hub fillet");
    m.def("build_shapes", [](const py::list& loops, const BladeGen::Blade::MeridionalCurve& hub,
                             const BladeGen::Blade::MeridionalCurve& casing, const InArray& inside, double tip_clearance,
                             double fillet, double extend) {
        auto l = loops_from(loops);
        Occt::SolidOptions o;
        o.hub = hub;
        o.casing = casing;
        o.tip_clearance = tip_clearance;
        o.fillet_radius = fillet;
        const Vec3 p = from_np<3>(inside.attr("reshape")(1, 3), "inside")[0];
        py::gil_scoped_release nogil;
        if (extend > 0.0) l = Occt::extend_loops(l, extend, extend);
        return Occt::build_shapes(l, o, p);
    }, py::arg("loops"), py::arg("hub"), py::arg("casing"), py::arg("inside_point"), py::arg("tip_clearance") = 0.0,
       py::arg("fillet_radius") = 0.0, py::arg("extend") = 0.0,
       "Solid from closed 3-D section loops (same point count), trimmed on hub / casing (x, r) curves");
    m.def("skin_loops", [](const py::list& loops) {
        const auto r = Occt::skin_loops(loops_from(loops));
        return py::make_tuple(Surface{r.surface}, r.max_error);
    }, py::arg("loops"), "C2 periodic NURBS surface through closed loops; returns (Surface, max interpolation error)");
    m.def("read_step", [](const std::string& p) { return Shape{Occt::read_step(p)}; }, py::arg("path"));
}

}  // namespace bladegen_py
