/// @file turbolenza/Discretization.hpp
/// @brief Finite differences of the stream-function equation: the normalized nine-point
///        equation at every interior node and the rows of the periodic PDE line.
///
/// Port of v1 discretization.py, production form "expanded" (TURBOlenza v1 commit 67599a6;
/// Aungier 2003 Eqs. 5-35 to 5-50).  Second-order central differences (5-35 .. 5-39) in Eq. 5-21,
/// divided by -H, H = 2 A / deta^2 + 2 C / dm^2 (Eq. 5-44):
///     west  = -(C/dm^2 - E/(2 dm)) / H    [5-45]     east  = -(C/dm^2 + E/(2 dm)) / H    [5-46]
///     south = -(A/deta^2 - D/(2 deta)) / H [5-47]    north = -(A/deta^2 + D/(2 deta)) / H [5-48]
///     mixed = B / (2 dm deta H)            [5-49]    rhs = -Q / H = 0                     [5-50]
/// center 1, NE = SW = +mixed, SE = NW = -mixed; NaN off the interior (i = 1 .. n_m - 2,
/// j = 1 .. n_eta - 2).  dm = length_scale grid.dm, deta = grid.deta.
/// The periodic PDE line (v1 build_periodic_pde_rows): Eq. 5-44 at (i, J), J = n_eta - 1, on the
/// periodic stations 1 .. f - 1 and l + 1 .. n_m - 2, the nodes above the line replaced by the
/// periodic continuation psi[k, J + 1] = psi[k, 1] + 1; the mixed term's m difference central
/// (standard), backward (3, -4, 1) at f - 1 (LE transition), forward (-3, 4, -1) at l + 1 (TE
/// transition); the TE row's psi_m and psi_m,m forward one-sided along j = J ((-3, 4, -1) / (2 dm),
/// (2, -5, 4, -1) / dm^2), its center then (2 A / deta^2 - 2 C / dm^2 + 3 E / (2 dm)) / H.
/// The v1 study option transition_form "flux" is not ported.
#pragma once

#include "turbolenza/BoundaryConditions.hpp"
#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Field.hpp"
#include "turbolenza/Metrics.hpp"

#include <functional>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 DiscretizationError.
class DiscretizationError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// v1 InteriorStencil (stencil weights, not the PDE coefficients).
struct InteriorStencil
{
    Field center, west, east, south, north, mixed, rhs, H;
    double dm = 0.0, deta = 0.0;
    [[nodiscard]] bool interior(int i, int j) const noexcept
    {
        return i >= 1 && i <= center.n_m - 2 && j >= 1 && j <= center.n_eta - 2;
    }
};

enum class PeriodicRowKind { Standard, LeTransition, TeTransition };
[[nodiscard]] const char* periodic_row_kind_name(PeriodicRowKind k) noexcept;  ///< v1 strings

struct PDETerm
{
    Node node;
    double coefficient = 0.0;
};

/// v1 PeriodicPDERow: sum coefficient psi[node] = rhs; H the divisor of the row.
struct PeriodicPDERow
{
    PeriodicRowKind kind = PeriodicRowKind::Standard;
    Node target;
    std::vector<PDETerm> terms;
    double rhs = 0.0;
    double H = 0.0;
    [[nodiscard]] double coefficient(Node n) const noexcept;
    [[nodiscard]] double residual(const Field& psi) const;
};

/// v1 build_interior_stencil (coefficients checked to belong to these metrics: C = A (L S*)^2
/// cos^2(beta), B = C tan(beta) / (L S*) to ROUND_OFF).
[[nodiscard]] InteriorStencil build_interior_stencil(const CoefficientFields& coefficients, const GridMetrics& metrics);

/// v1 build_periodic_pde_rows, in station order.
[[nodiscard]] std::vector<PeriodicPDERow> build_periodic_pde_rows(const CoefficientFields& coefficients,
                                                                  const GridMetrics& metrics,
                                                                  const BoundaryConditions& boundary_conditions);

/// v1 _interior_stencil on raw fields (manufactured-field tests).
[[nodiscard]] InteriorStencil interior_stencil(const Field& A, const Field& B, const Field& C, const Field& D,
                                               const Field& E, double dm, double deta);

/// v1 _periodic_pde_rows on raw fields; ghost(k) is the Elimination of psi[k, n_eta].
[[nodiscard]] std::vector<PeriodicPDERow> periodic_pde_rows(const Field& A, const Field& B, const Field& C,
                                                            const Field& D, const Field& E, double dm, double deta,
                                                            int first_blade_i, int last_blade_i,
                                                            const std::function<Elimination(int)>& ghost);

}  // namespace BladeGen::Turbolenza
