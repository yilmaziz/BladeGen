// Shared helpers of the BladeGen Python bindings.
#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "math/Vec.hpp"
#include "profile/Quality.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace bladegen_py {

namespace py = pybind11;
using BladeGen::Math::Vec2;
using BladeGen::Math::Vec3;
using InArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

inline py::array_t<double> to_np(const std::vector<double>& v)
{
    return py::array_t<double>(static_cast<py::ssize_t>(v.size()), v.data());
}

inline py::array_t<int> to_np(const std::vector<int>& v)
{
    return py::array_t<int>(static_cast<py::ssize_t>(v.size()), v.data());
}

template <int D>
py::array_t<double> to_np(const std::vector<BladeGen::Math::Vec<D>>& pts)
{
    py::array_t<double> a({static_cast<py::ssize_t>(pts.size()), static_cast<py::ssize_t>(D)});
    auto r = a.mutable_unchecked<2>();
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (int k = 0; k < D; ++k) r(static_cast<py::ssize_t>(i), k) = pts[i][k];
    return a;
}

template <int D>
std::vector<BladeGen::Math::Vec<D>> from_np(const InArray& a, const char* what)
{
    if (a.ndim() != 2 || a.shape(1) != D)
        throw std::invalid_argument(std::string(what) + ": expected an (N, " + std::to_string(D) + ") array");
    std::vector<BladeGen::Math::Vec<D>> out(static_cast<std::size_t>(a.shape(0)));
    auto r = a.unchecked<2>();
    for (py::ssize_t i = 0; i < a.shape(0); ++i)
        for (int k = 0; k < D; ++k) out[static_cast<std::size_t>(i)][k] = r(i, k);
    return out;
}

inline py::array_t<double> vec_np(const Vec2& v) { return to_np<2>(std::vector<Vec2>{v}).attr("reshape")(2); }
inline py::array_t<double> vec_np(const Vec3& v) { return to_np<3>(std::vector<Vec3>{v}).attr("reshape")(3); }

/// A curve sample and a quality report as dicts (bind_profile.cpp).
py::dict sample_dict(const BladeGen::Profile::CurveSample& c);
py::dict quality_dict(const BladeGen::Profile::QualityReport& r);

void bind_profile(py::module_& m);
void bind_naca(py::module_& m);
void bind_blade(py::module_& m);
void bind_fit(py::module_& m);
void bind_occt(py::module_& m);
void bind_design(py::module_& m);

}  // namespace bladegen_py
