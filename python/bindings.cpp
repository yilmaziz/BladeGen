// Python bindings of the BladeGen C++ engine (module bladegen._core).
//
// A thin layer: no geometry is computed here, every call goes to the C++
// library and every precondition is the library's own.  std::invalid_argument
// arrives in Python as ValueError, std::runtime_error as RuntimeError.
// Point lists are (N, 2) / (N, 3) float64 numpy arrays.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "bindings_common.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
    m.doc() = "BladeGen C++ engine: C2 blade profiles, 3-D blades, inverse fitting, OpenCASCADE solids";
    m.attr("__version__") = BLADEGEN_VERSION;
    bladegen_py::bind_profile(m);
    bladegen_py::bind_naca(m);
    bladegen_py::bind_blade(m);
    bladegen_py::bind_fit(m);
    bladegen_py::bind_occt(m);
    bladegen_py::bind_design(m);
}
