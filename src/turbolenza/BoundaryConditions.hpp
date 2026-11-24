/// @file turbolenza/BoundaryConditions.hpp
/// @brief The boundary conditions of the stream-function problem as local linear relations.
///
/// Port of v1 boundary_conditions.py (TURBOlenza v1 commit 67599a6; Aungier 2003 Sec. 5.3,
/// Fig. 5-3, Eqs. 5-29 to 5-33, the one-sided differences of Eqs. 5-40 to 5-43; the matrix
/// strategy of Aungier 2000 Ch. 13).  Every condition is sum a_k psi[node_k] = rhs:
///     lower wall   blade stations, j = 0:          psi = 0
///     upper wall   blade stations, j = n_eta - 1:  psi = 1
///     periodic     periodic stations:              psi[i, J] - psi[i, 0] = 1
///     upstream     i = 0, every j:                 -3 psi[0, j] + 4 psi[1, j] - psi[2, j] = 0   (psi_m = 0)
///     downstream   i = N, every j:                 3 psi[N, j] - 4 psi[N-1, j] + psi[N-2, j] = 0
/// Eliminations (Relation::solve_for): psi[0, j] = 4/3 psi[1, j] - 1/3 psi[2, j], likewise at the
/// exit.  The periodic continuation (periodic_ghost): psi[i, -1] = psi[i, J - 1] - 1 and
/// psi[i, n_eta] = psi[i, 1] + 1, on periodic stations only.  Requirements: n_eta >= 3, at least
/// three periodic stations each side, the topology sources agreeing (BoundaryConditionError).
#pragma once

#include "turbolenza/Field.hpp"
#include "turbolenza/GridData.hpp"

#include <compare>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

inline constexpr double kPsiLowerWall = 0.0;
inline constexpr double kPsiUpperWall = 1.0;
inline constexpr double kPsiPeriodicJump = kPsiUpperWall - kPsiLowerWall;

/// v1 BoundaryConditionError.
class BoundaryConditionError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// A node (i, j); j may be virtual (-1 or n_eta) in an Elimination of a periodic ghost.
struct Node
{
    int i = 0, j = 0;
    auto operator<=>(const Node&) const = default;
};

/// coefficient * psi[node].
struct Term
{
    Node node;
    double coefficient = 0.0;
};

/// psi[node] = sum of the terms + constant.
struct Elimination
{
    Node node;
    std::vector<Term> terms;
    double constant = 0.0;
    [[nodiscard]] double value(const Field& psi) const;
};

enum class RelationKind { LowerWall, UpperWall, Periodic, Upstream, Downstream };
[[nodiscard]] const char* relation_kind_name(RelationKind k) noexcept;  ///< v1 strings ("lower wall", ...)

/// sum coefficient * psi[node] = rhs; target the node it determines (none for the jump).
struct Relation
{
    RelationKind kind = RelationKind::LowerWall;
    std::vector<Term> terms;
    double rhs = 0.0;
    std::optional<Node> target;

    /// Throws BoundaryConditionError unless the nodes are distinct, coefficients finite non-zero.
    Relation(RelationKind kind_, std::vector<Term> terms_, double rhs_, std::optional<Node> target_);
    [[nodiscard]] bool has(Node n) const noexcept;
    [[nodiscard]] double residual(const Field& psi) const;
    [[nodiscard]] Elimination solve_for(Node node) const;
};

/// v1 BoundaryConditions: masks (row major, n_m x n_eta) and the relations of one grid.
struct BoundaryConditions
{
    std::shared_ptr<const Grid> grid;
    int n_m = 0, n_eta = 0;
    std::vector<char> lower_wall, upper_wall, periodic_lower, periodic_upper, upstream, downstream;
    std::vector<Relation> relations;  ///< lower walls, upper walls, periodic jumps, upstream, downstream

    [[nodiscard]] std::size_t at(int i, int j) const noexcept
    {
        return static_cast<std::size_t>(i) * static_cast<std::size_t>(n_eta) + static_cast<std::size_t>(j);
    }
    [[nodiscard]] bool blade_wall(int i, int j) const noexcept { return lower_wall[at(i, j)] || upper_wall[at(i, j)]; }
    [[nodiscard]] bool periodic_side(int i, int j) const noexcept
    {
        return periodic_lower[at(i, j)] || periodic_upper[at(i, j)];
    }
    [[nodiscard]] bool meridional_boundary(int i, int j) const noexcept { return upstream[at(i, j)] || downstream[at(i, j)]; }
    [[nodiscard]] bool boundary(int i, int j) const noexcept
    {
        return blade_wall(i, j) || periodic_side(i, j) || meridional_boundary(i, j);
    }
    [[nodiscard]] bool interior(int i, int j) const noexcept { return !boundary(i, j); }
    [[nodiscard]] bool periodic_line(int i, int j) const noexcept { return periodic_side(i, j) && !meridional_boundary(i, j); }
    [[nodiscard]] bool blade_station(int i) const noexcept { return lower_wall[at(i, 0)] != 0; }

    [[nodiscard]] Elimination upstream_elimination(int j) const;
    [[nodiscard]] Elimination downstream_elimination(int j) const;
    [[nodiscard]] Relation periodic_relation(int i, double jump = kPsiPeriodicJump) const;
    /// j = -1: f[i, -1] = f[i, n_eta - 2] - jump; j = n_eta: f[i, n_eta] = f[i, 1] + jump.
    [[nodiscard]] Elimination periodic_ghost(int i, int j, double jump = kPsiPeriodicJump) const;
};

/// v1 build_boundary_conditions.
[[nodiscard]] BoundaryConditions build_boundary_conditions(std::shared_ptr<const Grid> grid);

}  // namespace BladeGen::Turbolenza
