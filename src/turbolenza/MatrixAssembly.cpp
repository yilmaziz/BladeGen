// Port of v1 matrix_assembly.py, production path (see MatrixAssembly.hpp): full_index, node_of,
// assemble_system, assemble_from_coefficients, reconstruct, implied_residuals, _assemble, _resolve.
#include "turbolenza/MatrixAssembly.hpp"

#include "turbolenza/Detail.hpp"

#include <cmath>
#include <format>
#include <functional>
#include <map>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::sz;

const char* row_kind_name(RowKind k) noexcept
{
    switch (k) {
    case RowKind::Interior: return "interior";
    case RowKind::LowerWall: return "lower wall";
    case RowKind::UpperWall: return "upper wall";
    case RowKind::Periodic: return "periodic";
    case RowKind::PeriodicStandard: return "periodic standard";
    case RowKind::PeriodicLeTransition: return "periodic LE transition";
    case RowKind::PeriodicTeTransition: return "periodic TE transition";
    }
    return "?";
}

int full_index(int i, int j, int n_m, int n_eta)
{
    if (!(i >= 0 && i < n_m && j >= 0 && j < n_eta))
        throw MatrixAssemblyError(std::format("node ({}, {}) is off the {} x {} grid", i, j, n_m, n_eta));
    return i * n_eta + j;
}

Node node_of(int k, int n_m, int n_eta)
{
    if (k < 0 || k >= n_m * n_eta)
        throw MatrixAssemblyError(
            std::format("full index {} is off the {} x {} grid (0 .. {})", k, n_m, n_eta, n_m * n_eta - 1));
    return {k / n_eta, k % n_eta};
}

namespace {

struct LocalEquation
{
    RowKind kind;
    Node owner;
    std::vector<Term> terms;
    double rhs;
};

RowKind kind_of(PeriodicRowKind k)
{
    switch (k) {
    case PeriodicRowKind::Standard: return RowKind::PeriodicStandard;
    case PeriodicRowKind::LeTransition: return RowKind::PeriodicLeTransition;
    case PeriodicRowKind::TeTransition: return RowKind::PeriodicTeTransition;
    }
    return RowKind::PeriodicStandard;
}

RowKind kind_of(RelationKind k)
{
    switch (k) {
    case RelationKind::LowerWall: return RowKind::LowerWall;
    case RelationKind::UpperWall: return RowKind::UpperWall;
    default: return RowKind::Periodic;
    }
}

struct Resolved
{
    std::vector<std::pair<Node, double>> terms;  // unknowns only, insertion order
    double constant = 0.0;
};

// v1 _resolve: every elimination in terms of unknowns, dependencies first (depth first; a cycle
// is a MatrixAssemblyError).
void resolve(const std::map<Node, Elimination>& by_node, std::vector<Node>& order, std::map<Node, Resolved>& resolved)
{
    std::map<Node, int> state;  // 1 open, 2 done
    std::function<void(const Node&)> visit = [&](const Node& node) {
        const int s = state[node];
        if (s == 2) return;
        if (s == 1)
            throw MatrixAssemblyError(std::format("the eliminations form a cycle through ({}, {})", node.i, node.j));
        state[node] = 1;
        Resolved r;
        r.constant = by_node.at(node).constant;
        auto add = [&](const Node& u, double v) {
            for (auto& [n, a] : r.terms)
                if (n == u) {
                    a = a + v;
                    return;
                }
            r.terms.emplace_back(u, 0.0 + v);
        };
        for (const Term& t : by_node.at(node).terms) {
            if (by_node.count(t.node)) {
                visit(t.node);
                const Resolved& sub = resolved.at(t.node);
                for (const auto& [u, b] : sub.terms) add(u, t.coefficient * b);
                r.constant += t.coefficient * sub.constant;
            } else {
                add(t.node, t.coefficient);
            }
        }
        resolved[node] = std::move(r);
        state[node] = 2;
        order.push_back(node);
    };
    for (const auto& [node, e] : by_node) visit(node);
}

}  // namespace

AssemblyResult assemble_system(const InteriorStencil& stencil, const std::vector<PeriodicPDERow>& rows,
                               const BoundaryConditions& bcs)
{
    const int n_m = bcs.n_m, n_eta = bcs.n_eta;
    if (stencil.center.n_m != n_m || stencil.center.n_eta != n_eta)
        throw MatrixAssemblyError(std::format("the stencil is for a ({}, {}) grid, the boundary conditions for ({}, {})",
                                              stencil.center.n_m, stencil.center.n_eta, n_m, n_eta));
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (stencil.interior(i, j) != bcs.interior(i, j))
                throw MatrixAssemblyError("the stencil's interior nodes are not the interior nodes of the boundary "
                                          "conditions: the two were built for different grids");
    for (const PeriodicPDERow& row : rows) {
        (void)full_index(row.target.i, row.target.j, n_m, n_eta);
        if (!bcs.periodic_line(row.target.i, row.target.j))
            throw MatrixAssemblyError(std::format("the periodic PDE row at ({}, {}) is not on a periodic line of these "
                                                  "boundary conditions",
                                                  row.target.i, row.target.j));
    }

    std::map<Node, Elimination> by_node;
    for (const Relation& r : bcs.relations)
        if (r.kind == RelationKind::Upstream || r.kind == RelationKind::Downstream) {
            Elimination e = r.solve_for(*r.target);
            if (by_node.count(e.node))
                throw MatrixAssemblyError(std::format("node ({}, {}) is eliminated twice", e.node.i, e.node.j));
            by_node.emplace(e.node, std::move(e));
        }

    // the equations: the interior block first (v1 blocks [interior, others]), then the periodic
    // rows, the walls and the jumps
    std::vector<LocalEquation> equations;
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j) {
            if (!stencil.interior(i, j)) continue;
            LocalEquation eq{RowKind::Interior, {i, j}, {}, stencil.rhs(i, j)};
            for (int di = -1; di <= 1; ++di)
                for (int dj = -1; dj <= 1; ++dj) {
                    double w = 0.0;
                    if (di == 0 && dj == 0) w = stencil.center(i, j);
                    else if (di == -1 && dj == 0) w = stencil.west(i, j);
                    else if (di == 1 && dj == 0) w = stencil.east(i, j);
                    else if (di == 0 && dj == -1) w = stencil.south(i, j);
                    else if (di == 0 && dj == 1) w = stencil.north(i, j);
                    else if (di == dj) w = stencil.mixed(i, j);   // NE, SW
                    else w = -stencil.mixed(i, j);                // SE, NW
                    eq.terms.push_back(Term{{i + di, j + dj}, w});
                }
            equations.push_back(std::move(eq));
        }
    std::vector<char> owned(sz(n_m * n_eta), 0);
    for (int i = 0; i < n_m; ++i)
        for (int j = 0; j < n_eta; ++j)
            if (stencil.interior(i, j)) owned[sz(i * n_eta + j)] = 1;
    for (const PeriodicPDERow& row : rows) {
        LocalEquation eq{kind_of(row.kind), row.target, {}, row.rhs};
        for (const PDETerm& t : row.terms) eq.terms.push_back(Term{t.node, t.coefficient});
        owned[sz(row.target.i * n_eta + row.target.j)] = 1;
        equations.push_back(std::move(eq));
    }
    for (const Relation& r : bcs.relations)
        if (r.kind == RelationKind::LowerWall || r.kind == RelationKind::UpperWall) {
            equations.push_back(LocalEquation{kind_of(r.kind), *r.target, r.terms, r.rhs});
            owned[sz(r.target->i * n_eta + r.target->j)] = 1;
        }
    std::vector<Relation> implied;
    for (const Relation& r : bcs.relations) {
        if (r.kind != RelationKind::Periodic) continue;
        std::vector<Node> free;
        bool all_eliminated = true;
        for (const Term& t : r.terms) {
            const bool elim = by_node.count(t.node) > 0;
            if (!elim) all_eliminated = false;
            if (!elim && !owned[sz(t.node.i * n_eta + t.node.j)]) free.push_back(t.node);
        }
        if (all_eliminated) {
            implied.push_back(r);
        } else if (free.size() == 1) {
            equations.push_back(LocalEquation{RowKind::Periodic, free[0], r.terms, r.rhs});
            owned[sz(free[0].i * n_eta + free[0].j)] = 1;
        } else {
            throw MatrixAssemblyError(std::format("a periodic relation has {} image(s) without an equation: its owner "
                                                  "is not determined by the rows supplied",
                                                  free.size()));
        }
    }
    for (const LocalEquation& eq : equations) {
        if (!std::isfinite(eq.rhs)) throw MatrixAssemblyError("a local equation has a non-finite rhs");
        (void)full_index(eq.owner.i, eq.owner.j, n_m, n_eta);
        for (const Term& t : eq.terms) {
            if (!std::isfinite(t.coefficient)) throw MatrixAssemblyError("a local equation has a non-finite coefficient");
            (void)full_index(t.node.i, t.node.j, n_m, n_eta);
        }
    }
    for (const auto& [node, e] : by_node) {
        if (!std::isfinite(e.constant))
            throw MatrixAssemblyError(std::format("the elimination of ({}, {}) is not finite", node.i, node.j));
        for (const Term& t : e.terms) {
            if (!std::isfinite(t.coefficient))
                throw MatrixAssemblyError(std::format("the elimination of ({}, {}) is not finite", node.i, node.j));
            (void)full_index(t.node.i, t.node.j, n_m, n_eta);
        }
    }
    std::vector<Node> order;
    std::map<Node, Resolved> resolved;
    resolve(by_node, order, resolved);

    AssemblyResult out;
    out.n_m = n_m;
    out.n_eta = n_eta;
    const int n_full = n_m * n_eta;
    std::vector<char> active(sz(n_full), 1);
    for (const auto& [node, e] : by_node) active[sz(node.i * n_eta + node.j)] = 0;
    out.full_to_active.assign(sz(n_full), -1);
    for (int k = 0; k < n_full; ++k)
        if (active[sz(k)]) {
            out.full_to_active[sz(k)] = static_cast<int>(out.active_nodes.size());
            out.active_nodes.push_back(k);
        }
    const int n = static_cast<int>(out.active_nodes.size());

    // ownership: every unknown exactly one equation, none on an eliminated node
    std::vector<int> row_of_eq(equations.size());
    std::vector<int> count(sz(n), 0);
    for (std::size_t e = 0; e < equations.size(); ++e) {
        const Node o = equations[e].owner;
        const int r = out.full_to_active[sz(o.i * n_eta + o.j)];
        if (r < 0)
            throw MatrixAssemblyError(std::format("the {} equation of ({}, {}) stands on an eliminated node",
                                                  row_kind_name(equations[e].kind), o.i, o.j));
        row_of_eq[e] = r;
        ++count[sz(r)];
    }
    for (int r = 0; r < n; ++r) {
        const Node node = node_of(out.active_nodes[sz(r)], n_m, n_eta);
        if (count[sz(r)] > 1)
            throw MatrixAssemblyError(std::format("node ({}, {}) owns {} equations", node.i, node.j, count[sz(r)]));
        if (count[sz(r)] == 0)
            throw MatrixAssemblyError(std::format("the unknown ({}, {}) has no equation: the system would have {} "
                                                  "equations for {} unknowns",
                                                  node.i, node.j, equations.size(), n));
    }

    // F = L T, Q = q - L c
    out.constant.assign(sz(n_full), 0.0);
    for (const Node& node : order) out.constant[sz(node.i * n_eta + node.j)] = resolved.at(node).constant;
    std::vector<Eigen::Triplet<double, int>> triplets;
    triplets.reserve(equations.size() * 10);
    out.rhs.assign(sz(n), 0.0);
    out.row_kind.assign(sz(n), RowKind::Interior);
    out.substituted.assign(sz(n), 0);
    for (std::size_t e = 0; e < equations.size(); ++e) {
        const LocalEquation& eq = equations[e];
        const int r = row_of_eq[e];
        out.row_kind[sz(r)] = eq.kind;
        double Lc = 0.0;
        for (const Term& t : eq.terms) {
            const int k = t.node.i * n_eta + t.node.j;
            if (active[sz(k)]) {
                triplets.emplace_back(r, out.full_to_active[sz(k)], t.coefficient);
            } else {
                const Resolved& res = resolved.at(t.node);
                for (const auto& [u, alpha] : res.terms)
                    triplets.emplace_back(r, out.full_to_active[sz(u.i * n_eta + u.j)], t.coefficient * alpha);
                Lc += t.coefficient * res.constant;
                if (t.coefficient != 0.0) {
                    out.substituted[sz(r)] = 1;
                    ++out.substitutions;
                }
            }
        }
        out.rhs[sz(r)] = eq.rhs - Lc;
    }
    out.matrix.resize(n, n);
    out.matrix.setFromTriplets(triplets.begin(), triplets.end());
    out.matrix.prune([](auto, auto, const double& v) { return v != 0.0; });
    out.matrix.makeCompressed();

    // checks on the result
    for (int c = 0; c < out.matrix.outerSize(); ++c)
        for (SparseMatrix::InnerIterator it(out.matrix, c); it; ++it)
            if (!std::isfinite(it.value())) throw MatrixAssemblyError("the assembled matrix or right-hand side is not finite");
    for (double q : out.rhs)
        if (!std::isfinite(q)) throw MatrixAssemblyError("the assembled matrix or right-hand side is not finite");
    std::vector<int> row_nnz(sz(n), 0);
    for (int c = 0; c < out.matrix.outerSize(); ++c) {
        int col_nnz = 0;
        for (SparseMatrix::InnerIterator it(out.matrix, c); it; ++it) {
            ++row_nnz[sz(static_cast<int>(it.row()))];
            ++col_nnz;
        }
        if (col_nnz == 0) {
            const Node node = node_of(out.active_nodes[sz(c)], n_m, n_eta);
            throw MatrixAssemblyError(
                std::format("the unknown ({}, {}) appears in no equation (an empty column)", node.i, node.j));
        }
    }
    for (int r = 0; r < n; ++r)
        if (row_nnz[sz(r)] == 0) {
            const Node node = node_of(out.active_nodes[sz(r)], n_m, n_eta);
            throw MatrixAssemblyError(std::format("the row of ({}, {}) ({}) is empty", node.i, node.j,
                                                  row_kind_name(out.row_kind[sz(r)])));
        }
    for (const Node& node : order) out.eliminations.push_back(by_node.at(node));
    out.implied = std::move(implied);
    return out;
}

AssemblyResult assemble_from_coefficients(const CoefficientFields& coefficients, const GridMetrics& metrics,
                                          const BoundaryConditions& boundary_conditions)
{
    const InteriorStencil stencil = build_interior_stencil(coefficients, metrics);
    const std::vector<PeriodicPDERow> rows = build_periodic_pde_rows(coefficients, metrics, boundary_conditions);
    return assemble_system(stencil, rows, boundary_conditions);
}

Field reconstruct(const AssemblyResult& result, const std::vector<double>& x)
{
    if (static_cast<int>(x.size()) != result.n_unknowns())
        throw MatrixAssemblyError(
            std::format("the solution has shape ({},), the system {} unknowns", x.size(), result.n_unknowns()));
    for (double v : x)
        if (!std::isfinite(v)) throw MatrixAssemblyError("the solution is not finite");
    Field psi(result.n_m, result.n_eta, detail::nan());
    for (std::size_t r = 0; r < x.size(); ++r) psi.v[sz(result.active_nodes[r])] = x[r];
    for (const Elimination& e : result.eliminations) psi(e.node.i, e.node.j) = e.value(psi);
    for (double v : psi.v)
        if (!std::isfinite(v)) throw MatrixAssemblyError("a node is neither an unknown nor rebuilt by an elimination");
    return psi;
}

std::vector<double> implied_residuals(const AssemblyResult& result, const Field& psi)
{
    if (psi.n_m != result.n_m || psi.n_eta != result.n_eta)
        throw MatrixAssemblyError(std::format("psi has shape ({}, {}), the grid ({}, {})", psi.n_m, psi.n_eta,
                                              result.n_m, result.n_eta));
    std::vector<double> out;
    for (const Relation& r : result.implied) out.push_back(r.residual(psi));
    return out;
}

}  // namespace BladeGen::Turbolenza
