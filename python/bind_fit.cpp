// Bindings of the inverse-design layer and of the Rotor 67 data.
#include "bindings_common.hpp"

#include "fit/InverseFit.hpp"
#include "fit/NacaFit.hpp"
#include "fit/SectionFit.hpp"
#include "io/Rotor67.hpp"

namespace bladegen_py {

using namespace BladeGen::Fit;
namespace IO = BladeGen::IO;

namespace {

py::dict fit_dict(const FitResult& r)
{
    py::dict d;
    d["params"] = r.params;
    d["rms"] = r.rms;
    d["max"] = r.max;
    d["rms_ss"] = r.rms_ss;
    d["rms_ps"] = r.rms_ps;
    d["d_ss"] = to_np(r.d_ss);
    d["d_ps"] = to_np(r.d_ps);
    d["at_bound"] = r.at_bound;
    d["evaluations"] = r.evaluations;
    d["seconds"] = r.seconds;
    d["status"] = r.status;
    return d;
}

py::dict frame_dict(const Normalized& n)
{
    py::dict d;
    d["ss"] = to_np<2>(n.ss);
    d["ps"] = to_np<2>(n.ps);
    d["le"] = vec_np(n.le);
    d["te"] = vec_np(n.te);
    d["chord"] = n.chord;
    d["stagger"] = n.stagger;
    d["mirrored"] = n.mirrored;
    return d;
}

}  // namespace

void bind_fit(py::module_& m)
{
    py::class_<FitParam>(m, "FitParam", "A fitted scalar: profile parameter name, start, bounds")
        .def(py::init([](std::string n, double x0, double lo, double hi) { return FitParam{std::move(n), x0, lo, hi}; }),
             py::arg("name"), py::arg("x0"), py::arg("lo"), py::arg("hi"))
        .def_readwrite("name", &FitParam::name)
        .def_readwrite("x0", &FitParam::x0)
        .def_readwrite("lo", &FitParam::lo)
        .def_readwrite("hi", &FitParam::hi)
        .def("__repr__", [](const FitParam& p) {
            return "FitParam(" + p.name + ", x0=" + std::to_string(p.x0) + ", [" + std::to_string(p.lo) + ", " +
                   std::to_string(p.hi) + "])";
        });
    py::class_<FitOptions>(m, "FitOptions")
        .def(py::init<>())
        .def_readwrite("params", &FitOptions::params)
        .def_readwrite("stages", &FitOptions::stages)
        .def_readwrite("max_iterations", &FitOptions::max_iterations)
        .def_readwrite("fairness_weight", &FitOptions::fairness_weight)
        .def_readwrite("n_edge", &FitOptions::n_edge)
        .def_readwrite("n_body", &FitOptions::n_body);

    m.def("default_fit_params", &default_fit_params, py::arg("start"));
    m.def("naca65_start", &naca65_start, py::arg("naca"),
          "Start parameters for fitting custom parameters to a NACA 65 section");
    m.def("naca65_fit_options", &naca65_fit_options);
    m.def("fit_naca65", [](const BladeGen::Naca::Naca65Section& s, const FitOptions& opt, int n_target) {
        FitResult r;
        {
            py::gil_scoped_release nogil;
            r = fit_naca65(s, opt, n_target);
        }
        return fit_dict(r);
    }, py::arg("section"), py::arg("options") = naca65_fit_options(), py::arg("n_target") = 100,
       "NACA 65 + Custom: custom profile parameters fitted to a NACA 65 section (3-5 s)");
    m.def("fit_profile", [](const InArray& ss, const InArray& ps, const BladeGen::Profile::ProfileParams& start,
                            const FitOptions& opt) {
        const auto a = from_np<2>(ss, "ss"), b = from_np<2>(ps, "ps");
        FitResult r;
        {
            py::gil_scoped_release nogil;
            r = fit_profile(a, b, start, opt);
        }
        return fit_dict(r);
    }, py::arg("ss"), py::arg("ps"), py::arg("start"), py::arg("options") = FitOptions{},
       "Inverse design: profile parameters reproducing suction / pressure side points (unstaggered frame)");
    m.def("side_distances", [](const BladeGen::Profile::Profile2D& p, const InArray& ss, const InArray& ps) {
        std::vector<double> a, b;
        side_distances(p, from_np<2>(ss, "ss"), from_np<2>(ps, "ps"), a, b, 400, 1600);
        return py::make_tuple(to_np(a), to_np(b));
    }, py::arg("profile"), py::arg("ss"), py::arg("ps"));
    m.def("normalize_contour", [](const InArray& pts) { return frame_dict(normalize_contour(from_np<2>(pts, "points"))); },
          py::arg("points"));
    m.def("fit_section", [](const InArray& loop, const BladeGen::Blade::MeridionalCurve& sl, double span,
                            const BladeGen::Profile::ProfileParams& start, const FitOptions& opt) {
        const auto l = from_np<2>(loop, "loop");
        SectionFitResult r;
        {
            py::gil_scoped_release nogil;
            r = fit_section(l, sl, span, start, opt);
        }
        py::dict d;
        d["frame"] = frame_dict(r.frame);
        d["fit"] = fit_dict(r.fit);
        d["station"] = r.station;
        d["published_xyz"] = to_np<3>(r.published_xyz);
        d["fitted_xyz"] = to_np<3>(r.fitted_xyz);
        d["rms_3d"] = r.rms_3d;
        d["max_3d"] = r.max_3d;
        return d;
    }, py::arg("conformal_loop"), py::arg("streamline"), py::arg("span"), py::arg("start"), py::arg("options") = FitOptions{});

    py::class_<IO::EdgeCircle>(m, "EdgeCircle")
        .def_readonly("rho", &IO::EdgeCircle::rho)
        .def_readonly("rc", &IO::EdgeCircle::rc)
        .def_readonly("zc", &IO::EdgeCircle::zc)
        .def_readonly("mc", &IO::EdgeCircle::mc)
        .def_readonly("thc", &IO::EdgeCircle::thc);
    py::class_<IO::Rotor67Section>(m, "Rotor67Section", "One published NASA Rotor 67 section (TP-2879)")
        .def_readonly("index", &IO::Rotor67Section::index)
        .def_readonly("le", &IO::Rotor67Section::le)
        .def_readonly("te", &IO::Rotor67Section::te)
        .def_readonly("streamline", &IO::Rotor67Section::streamline)
        .def_property_readonly("z", [](const IO::Rotor67Section& s) { return to_np(s.z); })
        .def_property_readonly("r", [](const IO::Rotor67Section& s) { return to_np(s.r); })
        .def_property_readonly("theta1", [](const IO::Rotor67Section& s) { return to_np(s.th1); })
        .def_property_readonly("theta2", [](const IO::Rotor67Section& s) { return to_np(s.th2); })
        .def_property_readonly("m", [](const IO::Rotor67Section& s) { return to_np(s.m); })
        .def("profile_m_theta", [](const IO::Rotor67Section& s, int a, int b) { return to_np<2>(s.profile_m_theta(a, b)); },
             py::arg("n_surface") = 161, py::arg("n_arc") = 33)
        .def("profile_conformal", [](const IO::Rotor67Section& s, int a, int b) { return to_np<2>(s.profile_conformal(a, b)); },
             py::arg("n_surface") = 161, py::arg("n_arc") = 33)
        .def("profile_xyz", [](const IO::Rotor67Section& s, int a, int b) { return to_np<3>(s.profile_xyz(a, b)); },
             py::arg("n_surface") = 161, py::arg("n_arc") = 33);
    py::class_<IO::Rotor67Data>(m, "Rotor67Data")
        .def_readonly("sections", &IO::Rotor67Data::sections)
        .def_readonly("n_blades", &IO::Rotor67Data::n_blades)
        .def_readonly("tip_clearance", &IO::Rotor67Data::tip_clearance)
        .def_property_readonly("hub", [](const IO::Rotor67Data& d) { return to_np<2>(d.hub); })
        .def_property_readonly("shroud", [](const IO::Rotor67Data& d) { return to_np<2>(d.shroud); })
        .def("hub_curve", &IO::Rotor67Data::hub_curve)
        .def("shroud_curve", &IO::Rotor67Data::shroud_curve);
    m.def("read_rotor67", &IO::read_rotor67, py::arg("path"), "NASA Rotor 67 sections and flowpath (metres)");
}

}  // namespace bladegen_py
