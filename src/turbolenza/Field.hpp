/// @file turbolenza/Field.hpp
/// @brief Node fields of the (m, eta) grid of TURBOlenza's blade-to-blade solver.
///
/// TURBOlenza (the blade-to-blade compressible potential-flow solver after Aungier 2003,
/// Ch. 5) is a C++ port of the verified Python package turbolenza.b2b (v1, commit 67599a6).
/// Every field is a node array of shape (n_m, n_eta) stored row major, i (along m) slowest
/// and j (along eta) fastest: index i * n_eta + j.  That is the C order of v1's NumPy arrays
/// and the node number k = i n_eta + j of its sparse system (matrix_assembly.py), so the port
/// can be compared with v1 entry by entry.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// A scalar field on the (n_m, n_eta) nodes, row major (j fastest).
struct Field
{
    int n_m = 0, n_eta = 0;
    std::vector<double> v;

    Field() = default;
    Field(int nm, int ne, double value = 0.0)
        : n_m(nm), n_eta(ne), v(static_cast<std::size_t>(nm) * static_cast<std::size_t>(ne), value)
    {
        if (nm < 0 || ne < 0) throw std::invalid_argument("Field: negative shape");
    }

    [[nodiscard]] std::size_t index(int i, int j) const noexcept
    {
        return static_cast<std::size_t>(i) * static_cast<std::size_t>(n_eta) + static_cast<std::size_t>(j);
    }
    [[nodiscard]] double& operator()(int i, int j) noexcept { return v[index(i, j)]; }
    [[nodiscard]] double operator()(int i, int j) const noexcept { return v[index(i, j)]; }
    [[nodiscard]] std::size_t size() const noexcept { return v.size(); }
    [[nodiscard]] bool same_shape(const Field& o) const noexcept { return n_m == o.n_m && n_eta == o.n_eta; }
};

/// A point of the cascade plane: m (axial, the grid's first coordinate) and y (tangential,
/// r theta of the planar cascade), in the active length unit.
struct Point
{
    double m = 0.0, y = 0.0;
};

}  // namespace BladeGen::Turbolenza
