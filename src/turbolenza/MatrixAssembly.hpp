/// @file turbolenza/MatrixAssembly.hpp
/// @brief The local equations gathered into one square sparse system F Psi = Q.
///
/// Port of v1 matrix_assembly.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Sec. 5.4; Aungier
/// 2000 Ch. 13: the node numbering, [F_kn][Psi_n] = [Q_k], the far-boundary values eliminated into
/// the adjacent equations and rebuilt after the solve).  No mathematics of its own:
///  * eliminated: every node of the two far stations by its meridional relation solved for it
///    (psi[0, j] = 4/3 psi[1, j] - 1/3 psi[2, j], likewise at the exit); the far stations' periodic
///    jumps are then implied (kept in `implied`, carried by no row);
///  * every other node is one unknown with one equation: interior Eq. 5-44, wall Dirichlet,
///    periodic PDE row at (i, n_eta - 1), periodic jump at (i, 0);
///  * numbering: full index k = i n_eta + j; the unknowns in increasing k (active index r);
///  * substitution psi_b = c_b + sum alpha_br psi_r, F = L T, Q = q - L c; entries meeting on one
///    column summed, exact zeros dropped; no row scaled.
/// The v1 study option weight_override (transition_form "flux") is not ported.
#pragma once

#include "turbolenza/BoundaryConditions.hpp"
#include "turbolenza/Coefficients.hpp"
#include "turbolenza/Discretization.hpp"

#include <Eigen/SparseCore>

#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// v1 MatrixAssemblyError.
class MatrixAssemblyError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// The equation of a row: v1 row_kind strings ("interior", "lower wall", "upper wall",
/// "periodic", "periodic standard", "periodic LE transition", "periodic TE transition").
enum class RowKind { Interior, LowerWall, UpperWall, Periodic, PeriodicStandard, PeriodicLeTransition, PeriodicTeTransition };
[[nodiscard]] const char* row_kind_name(RowKind k) noexcept;

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

/// v1 AssemblyResult (without the expansion matrix T, which only v1's own checks read).
struct AssemblyResult
{
    SparseMatrix matrix;                ///< F (N, N), compressed, no explicit zeros
    std::vector<double> rhs;            ///< Q (N)
    int n_m = 0, n_eta = 0;
    std::vector<int> active_nodes;      ///< full index of unknown r
    std::vector<int> full_to_active;    ///< unknown of a node, -1 if eliminated
    std::vector<RowKind> row_kind;      ///< the equation of row r
    std::vector<char> substituted;      ///< row r read an eliminated node
    std::vector<Elimination> eliminations;  ///< every eliminated node, dependencies first
    std::vector<Relation> implied;      ///< relations carried by no row (the far stations' jumps)
    std::vector<double> constant;       ///< c (n_m n_eta)
    int substitutions = 0;

    [[nodiscard]] int n_unknowns() const noexcept { return static_cast<int>(matrix.rows()); }
};

/// v1 full_index / node_of (MatrixAssemblyError off the grid, no wrap).
[[nodiscard]] int full_index(int i, int j, int n_m, int n_eta);
[[nodiscard]] Node node_of(int k, int n_m, int n_eta);

/// v1 assemble_system: the interior rows of stencil, the periodic PDE rows and the relations.
[[nodiscard]] AssemblyResult assemble_system(const InteriorStencil& stencil, const std::vector<PeriodicPDERow>& rows,
                                             const BoundaryConditions& boundary_conditions);
/// v1 assemble_from_coefficients ("expanded").
[[nodiscard]] AssemblyResult assemble_from_coefficients(const CoefficientFields& coefficients, const GridMetrics& metrics,
                                                        const BoundaryConditions& boundary_conditions);

/// v1 reconstruct: psi on the whole grid (unknowns copied, eliminated nodes rebuilt).
[[nodiscard]] Field reconstruct(const AssemblyResult& result, const std::vector<double>& solution);
/// v1 implied_residuals.
[[nodiscard]] std::vector<double> implied_residuals(const AssemblyResult& result, const Field& psi);

}  // namespace BladeGen::Turbolenza
