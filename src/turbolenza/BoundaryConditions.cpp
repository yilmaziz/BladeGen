// Port of v1 boundary_conditions.py (see BoundaryConditions.hpp): Relation, Elimination,
// BoundaryConditions and build_boundary_conditions / _boundary_conditions / _periodic.
#include "turbolenza/BoundaryConditions.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <string>

namespace BladeGen::Turbolenza {

using detail::sz;

namespace {

constexpr double kForward[3] = {-3.0, 4.0, -1.0};   // x 1/(2 dm) at i = 0, 1, 2 [Eqs. 5-40 to 5-43]
constexpr double kBackward[3] = {3.0, -4.0, 1.0};   // x 1/(2 dm) at i = N, N-1, N-2

std::string nodes_text(const std::vector<Term>& terms)
{
    std::string s = "(";
    for (std::size_t k = 0; k < terms.size(); ++k)
        s += std::format("{}({}, {})", k ? ", " : "", terms[k].node.i, terms[k].node.j);
    return s + ")";
}

Relation periodic(int i, int n_eta, double jump)
{
    return Relation(RelationKind::Periodic, {Term{{i, n_eta - 1}, 1.0}, Term{{i, 0}, -1.0}}, jump, std::nullopt);
}

double real(double value, const char* name)
{
    if (!std::isfinite(value))
        throw BoundaryConditionError(std::format("{} must be a finite real number, got {}", name, detail::py_repr(value)));
    return value;
}

}  // namespace

const char* relation_kind_name(RelationKind k) noexcept
{
    switch (k) {
    case RelationKind::LowerWall: return "lower wall";
    case RelationKind::UpperWall: return "upper wall";
    case RelationKind::Periodic: return "periodic";
    case RelationKind::Upstream: return "upstream";
    case RelationKind::Downstream: return "downstream";
    }
    return "?";
}

double Elimination::value(const Field& psi) const
{
    double s = 0.0;
    for (const Term& t : terms) s += t.coefficient * psi(t.node.i, t.node.j);
    return s + constant;
}

Relation::Relation(RelationKind kind_, std::vector<Term> terms_, double rhs_, std::optional<Node> target_)
    : kind(kind_), terms(std::move(terms_)), rhs(rhs_), target(target_)
{
    bool distinct = !terms.empty();
    for (std::size_t a = 0; a < terms.size(); ++a)
        for (std::size_t b = a + 1; b < terms.size(); ++b)
            if (terms[a].node == terms[b].node) distinct = false;
    if (!distinct)
        throw BoundaryConditionError(
            std::format("a {} relation needs distinct nodes, got {}", relation_kind_name(kind), nodes_text(terms)));
    for (const Term& t : terms)
        if (!(std::isfinite(t.coefficient) && t.coefficient != 0.0) || !std::isfinite(rhs))
            throw BoundaryConditionError(std::format("a {} relation needs finite non-zero coefficients and a finite rhs",
                                                     relation_kind_name(kind)));
    if (target && !has(*target))
        throw BoundaryConditionError(
            std::format("the target ({}, {}) is not a node of the relation", target->i, target->j));
}

bool Relation::has(Node n) const noexcept
{
    for (const Term& t : terms)
        if (t.node == n) return true;
    return false;
}

double Relation::residual(const Field& psi) const
{
    double s = 0.0;
    for (const Term& t : terms) s += t.coefficient * psi(t.node.i, t.node.j);
    return s - rhs;
}

// v1 Relation.solve_for: psi[node] = sum (-a_k / a) psi[k] + rhs / a.
Elimination Relation::solve_for(Node node) const
{
    const Term* own = nullptr;
    for (const Term& t : terms)
        if (t.node == node) {
            own = &t;
            break;
        }
    if (!own)
        throw BoundaryConditionError(std::format("node ({}, {}) is not in this {} relation ({})", node.i, node.j,
                                                 relation_kind_name(kind), nodes_text(terms)));
    const double a = own->coefficient;
    Elimination e{node, {}, rhs / a + 0.0};
    for (const Term& t : terms)
        if (!(t.node == node)) e.terms.push_back(Term{t.node, -t.coefficient / a});
    return e;
}

Elimination BoundaryConditions::upstream_elimination(int j) const
{
    if (!(0 <= j && j < n_eta))
        throw BoundaryConditionError(std::format("j = {} is not one of the {} eta nodes", j, n_eta));
    for (const Relation& r : relations)
        if (r.kind == RelationKind::Upstream && r.target && *r.target == Node{0, j}) return r.solve_for({0, j});
    throw BoundaryConditionError(std::format("no upstream relation at (0, {})", j));
}

Elimination BoundaryConditions::downstream_elimination(int j) const
{
    if (!(0 <= j && j < n_eta))
        throw BoundaryConditionError(std::format("j = {} is not one of the {} eta nodes", j, n_eta));
    const int N = n_m - 1;
    for (const Relation& r : relations)
        if (r.kind == RelationKind::Downstream && r.target && *r.target == Node{N, j}) return r.solve_for({N, j});
    throw BoundaryConditionError(std::format("no downstream relation at ({}, {})", N, j));
}

namespace {

void periodic_station(const BoundaryConditions& b, int i)
{
    if (!(0 <= i && i < b.n_m))
        throw BoundaryConditionError(std::format("station {} is not one of the {} stations", i, b.n_m));
    if (!b.periodic_lower[b.at(i, 0)])
        throw BoundaryConditionError(std::format(
            "station {} is a blade station: its eta = 0 and 1 nodes are walls, and below the lower wall and above "
            "the upper wall lies the blade, so there is no periodic relation or continuation there",
            i));
}

}  // namespace

Relation BoundaryConditions::periodic_relation(int i, double jump) const
{
    periodic_station(*this, i);
    return periodic(i, n_eta, real(jump, "jump"));
}

Elimination BoundaryConditions::periodic_ghost(int i, int j, double jump) const
{
    periodic_station(*this, i);
    jump = real(jump, "jump");
    if (j == -1) return Elimination{{i, -1}, {Term{{i, n_eta - 2}, 1.0}}, -jump};
    if (j == n_eta) return Elimination{{i, n_eta}, {Term{{i, 1}, 1.0}}, jump};
    throw BoundaryConditionError(
        std::format("a periodic ghost node is j = -1 or j = n_eta = {}, got j = {}", n_eta, j));
}

BoundaryConditions build_boundary_conditions(std::shared_ptr<const Grid> grid)
{
    if (!grid) throw BoundaryConditionError("build_boundary_conditions takes a Grid, got None");
    const Grid& g = *grid;
    for (std::size_t k = 0; k < g.m.size(); ++k)
        if (!std::isfinite(g.m[k])) throw BoundaryConditionError(std::format("m is not finite at index {}", k));
    for (std::size_t k = 0; k < g.eta.size(); ++k)
        if (!std::isfinite(g.eta[k])) throw BoundaryConditionError(std::format("eta is not finite at index {}", k));
    const int n_m = g.n_m(), n_eta = g.n_eta();
    if (n_eta < 3)
        throw BoundaryConditionError(std::format("{} eta nodes: the walls and the periodic continuation need an interior "
                                                 "eta node, n_eta >= 3",
                                                 n_eta));
    const int f = g.first_blade_i, l = g.last_blade_i;
    if (!(1 <= f && f <= l && l <= n_m - 2))
        throw BoundaryConditionError(std::format("blade stations {} to {} do not fit {} stations with periodic stations "
                                                 "on both sides",
                                                 f, l, n_m));
    bool same = static_cast<int>(g.wall_identity.stations.size()) == l - f + 1;
    for (std::size_t k = 0; same && k < g.wall_identity.stations.size(); ++k)
        same = g.wall_identity.stations[k] == f + static_cast<int>(k);
    if (!same)
        throw BoundaryConditionError(std::format("the wall identity's stations are not first_blade_i .. last_blade_i = {} "
                                                 ".. {}: the topology sources disagree",
                                                 f, l));
    // grid.region agrees with f, l by construction (Grid::region derives from them)
    if (f < 3 || n_m - 1 - l < 3) {
        const bool inlet = f < 3;
        throw BoundaryConditionError(std::format(
            "the {} extension holds {} periodic station(s): the second-order one-sided psi_m = 0 relation reaches "
            "three stations, which must all be periodic (lengthen the extension or raise n_passage); no first-order "
            "fallback",
            inlet ? "inlet" : "exit", inlet ? f : n_m - 1 - l));
    }

    BoundaryConditions b;
    b.grid = grid;
    b.n_m = n_m;
    b.n_eta = n_eta;
    const std::size_t n = sz(n_m) * sz(n_eta);
    for (auto* mask : {&b.lower_wall, &b.upper_wall, &b.periodic_lower, &b.periodic_upper, &b.upstream, &b.downstream})
        mask->assign(n, 0);
    for (int i = 0; i < n_m; ++i) {
        const bool blade = f <= i && i <= l;
        (blade ? b.lower_wall : b.periodic_lower)[b.at(i, 0)] = 1;
        (blade ? b.upper_wall : b.periodic_upper)[b.at(i, n_eta - 1)] = 1;
    }
    for (int j = 0; j < n_eta; ++j) {
        b.upstream[b.at(0, j)] = 1;
        b.downstream[b.at(n_m - 1, j)] = 1;
    }
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (b.blade_wall(i, j) && (b.periodic_side(i, j) || b.meridional_boundary(i, j)))
                throw BoundaryConditionError(
                    "a wall node is also periodic or on a far station: contradictory classification");

    const int N = n_m - 1;
    for (int i = f; i <= l; ++i)
        b.relations.emplace_back(RelationKind::LowerWall, std::vector<Term>{Term{{i, 0}, 1.0}}, kPsiLowerWall, Node{i, 0});
    for (int i = f; i <= l; ++i)
        b.relations.emplace_back(RelationKind::UpperWall, std::vector<Term>{Term{{i, n_eta - 1}, 1.0}}, kPsiUpperWall,
                                 Node{i, n_eta - 1});
    for (int i = 0; i < n_m; ++i)
        if (i < f || i > l) b.relations.push_back(periodic(i, n_eta, kPsiPeriodicJump));
    for (int j = 0; j < n_eta; ++j)
        b.relations.emplace_back(RelationKind::Upstream,
                                 std::vector<Term>{Term{{0, j}, kForward[0]}, Term{{1, j}, kForward[1]},
                                                   Term{{2, j}, kForward[2]}},
                                 0.0, Node{0, j});
    for (int j = 0; j < n_eta; ++j)
        b.relations.emplace_back(RelationKind::Downstream,
                                 std::vector<Term>{Term{{N, j}, kBackward[0]}, Term{{N - 1, j}, kBackward[1]},
                                                   Term{{N - 2, j}, kBackward[2]}},
                                 0.0, Node{N, j});
    return b;
}

}  // namespace BladeGen::Turbolenza
