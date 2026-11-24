// Port of v1 discretization.py, production form (see Discretization.hpp): build_interior_stencil,
// _check_belongs, build_periodic_pde_rows, _periodic_pde_rows, _interior_stencil, _spacing.
#include "turbolenza/Discretization.hpp"

#include "turbolenza/Detail.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza {

using detail::nan;
using detail::py_repr;
using detail::sz;

const char* periodic_row_kind_name(PeriodicRowKind k) noexcept
{
    switch (k) {
    case PeriodicRowKind::Standard: return "periodic standard";
    case PeriodicRowKind::LeTransition: return "periodic LE transition";
    case PeriodicRowKind::TeTransition: return "periodic TE transition";
    }
    return "?";
}

double PeriodicPDERow::coefficient(Node n) const noexcept
{
    for (const PDETerm& t : terms)
        if (t.node == n) return t.coefficient;
    return 0.0;
}

double PeriodicPDERow::residual(const Field& psi) const
{
    double s = 0.0;
    for (const PDETerm& t : terms) s += t.coefficient * psi(t.node.i, t.node.j);
    return s - rhs;
}

namespace {

double spacing(double value, const char* name)
{
    if (!(std::isfinite(value) && value > 0.0))
        throw DiscretizationError(
            std::format("{} = {}: a grid spacing must be finite and positive", name, py_repr(value)));
    return value;
}

// The pointwise weights of Eqs. 5-45 to 5-49 at one node.
struct Weights
{
    double H, west, east, south, north, mixed;
};

Weights weights_at(double A, double B, double C, double D, double E, double dm, double deta)
{
    const double dm2 = dm * dm, deta2 = deta * deta;
    Weights w{};
    w.H = 2.0 * A / deta2 + 2.0 * C / dm2;
    w.west = -(C / dm2 - E / (2.0 * dm)) / w.H;        // Eq. 5-45
    w.east = -(C / dm2 + E / (2.0 * dm)) / w.H;        // Eq. 5-46
    w.south = -(A / deta2 - D / (2.0 * deta)) / w.H;   // Eq. 5-47
    w.north = -(A / deta2 + D / (2.0 * deta)) / w.H;   // Eq. 5-48
    w.mixed = B / (2.0 * dm * deta * w.H);             // Eq. 5-49
    return w;
}

void check_fields(const Field& A, const Field& B, const Field& C, const Field& D, const Field& E)
{
    const std::pair<const char*, const Field*> all[] = {{"A", &A}, {"B", &B}, {"C", &C}, {"D", &D}, {"E", &E}};
    for (const auto& [name, f] : all) {
        if (!f->same_shape(A))
            throw DiscretizationError(std::format("{} has shape ({}, {}): A to E must be (n_m, n_eta) arrays of one shape",
                                                  name, f->n_m, f->n_eta));
        for (int i = 0; i < f->n_m; ++i)
            for (int j = 0; j < f->n_eta; ++j)
                if (!std::isfinite((*f)(i, j)))
                    throw DiscretizationError(std::format("{} is not finite at node ({}, {})", name, i, j));
    }
    if (A.n_m < 3 || A.n_eta < 3)
        throw DiscretizationError(std::format(
            "{} x {} nodes hold no interior node: the nine-point stencil needs three a side", A.n_m, A.n_eta));
}

Weights checked_weights(const Field& A, const Field& B, const Field& C, const Field& D, const Field& E, int i, int j,
                        double dm, double deta)
{
    const Weights w = weights_at(A(i, j), B(i, j), C(i, j), D(i, j), E(i, j), dm, deta);
    if (!(std::isfinite(w.H) && w.H > 0.0))
        throw DiscretizationError(std::format("H = 2 A / deta^2 + 2 C / dm^2 = {} at node ({}, {}) (A = {:.6g}, C = "
                                              "{:.6g}): it must be finite and positive",
                                              py_repr(w.H), i, j, A(i, j), C(i, j)));
    const std::pair<const char*, double> all[] = {
        {"west", w.west}, {"east", w.east}, {"south", w.south}, {"north", w.north}, {"mixed", w.mixed}};
    for (const auto& [name, v] : all)
        if (!std::isfinite(v))
            throw DiscretizationError(
                std::format("the {} weight is not finite at node ({}, {}) (H = {:.6g})", name, i, j, w.H));
    return w;
}

// v1 _check_belongs.
void check_belongs(const CoefficientFields& co, const GridMetrics& metrics)
{
    const Field &A = co.A, &B = co.B, &C = co.C;
    const double L = co.length_scale();
    for (int i = 0; i < A.n_m; ++i) {
        const double S = L * metrics.S[sz(i)];
        for (int j = 0; j < A.n_eta; ++j) {
            const double c = metrics.cos_beta(i, j), t = metrics.tan_beta(i, j);
            const double Sc = S * c;
            const double off_C = std::abs(A(i, j) * (Sc * Sc) - C(i, j)) / std::abs(C(i, j));
            const double scale = std::abs(C(i, j) / S) * std::max(1.0, std::abs(t));
            const double off_B = std::abs(C(i, j) * t / S - B(i, j)) / scale;
            const char* name = nullptr;
            double off = 0.0;
            if (!(off_C <= kRoundOff)) {
                name = "C / A = (L S*)^2 cos^2(beta)";
                off = off_C;
            }
            if (!name && !(off_B <= kRoundOff)) {
                name = "B / C = tan(beta) / (L S*)";
                off = off_B;
            }
            if (name)
                throw DiscretizationError(std::format(
                    "the coefficient fields do not belong to these metrics: {} fails at node ({}, {}) by {:.3g} "
                    "(relative). Rebuild the coefficients from these metrics (after grid.rebuild_exit too), at the "
                    "length scale meant",
                    name, i, j, off));
        }
    }
}

void check_shape(const CoefficientFields& co, const Grid& grid)
{
    if (co.A.n_m != grid.n_m() || co.A.n_eta != grid.n_eta())
        throw DiscretizationError(std::format("the coefficient fields are ({}, {}), the grid has ({}, {}) nodes",
                                              co.A.n_m, co.A.n_eta, grid.n_m(), grid.n_eta()));
}

}  // namespace

InteriorStencil interior_stencil(const Field& A, const Field& B, const Field& C, const Field& D, const Field& E,
                                 double dm, double deta)
{
    dm = spacing(dm, "dm");
    deta = spacing(deta, "deta");
    check_fields(A, B, C, D, E);
    const int n_m = A.n_m, n_eta = A.n_eta;
    InteriorStencil s;
    for (Field* f : {&s.center, &s.west, &s.east, &s.south, &s.north, &s.mixed, &s.rhs, &s.H})
        *f = Field(n_m, n_eta, nan());
    s.dm = dm;
    s.deta = deta;
    for (int i = 1; i < n_m - 1; ++i)
        for (int j = 1; j < n_eta - 1; ++j) {
            const Weights w = checked_weights(A, B, C, D, E, i, j, dm, deta);
            s.H(i, j) = w.H;
            s.west(i, j) = w.west;
            s.east(i, j) = w.east;
            s.south(i, j) = w.south;
            s.north(i, j) = w.north;
            s.mixed(i, j) = w.mixed;
            s.center(i, j) = 1.0;
            s.rhs(i, j) = 0.0;  // Eq. 5-50: -Q / H with Q = 0
        }
    return s;
}

InteriorStencil build_interior_stencil(const CoefficientFields& coefficients, const GridMetrics& metrics)
{
    if (!metrics.grid) throw DiscretizationError("build_interior_stencil takes GridMetrics with their grid");
    const Grid& grid = *metrics.grid;
    check_shape(coefficients, grid);
    check_belongs(coefficients, metrics);
    return interior_stencil(coefficients.A, coefficients.B, coefficients.C, coefficients.D, coefficients.E,
                            coefficients.length_scale() * grid.dm(), grid.deta());
}

std::vector<PeriodicPDERow> periodic_pde_rows(const Field& A, const Field& B, const Field& C, const Field& D,
                                              const Field& E, double dm, double deta, int first_blade_i,
                                              int last_blade_i, const std::function<Elimination(int)>& ghost)
{
    dm = spacing(dm, "dm");
    deta = spacing(deta, "deta");
    check_fields(A, B, C, D, E);
    const int n_m = A.n_m, n_eta = A.n_eta;
    const int J = n_eta - 1;
    const int f = first_blade_i, l = last_blade_i;
    if (!(1 <= f && f <= l && l <= n_m - 2))
        throw DiscretizationError(std::format("blade stations {} to {} do not fit {} stations", f, l, n_m));
    if (f < 3 || n_m - 1 - l < 4)
        throw DiscretizationError(std::format(
            "{} inlet and {} exit periodic stations: the second-order one-sided differences of the transition rows "
            "need three on the inlet side and four on the exit side (the TE row's psi_m,m reaches last_blade_i + 4; no "
            "first-order fallback)",
            f, n_m - 1 - l));
    // the pointwise weights at (i, J), every station 1 .. n_m - 2 (v1 evaluates the whole column)
    std::vector<Weights> line(sz(n_m));
    for (int i = 1; i < n_m - 1; ++i) line[sz(i)] = checked_weights(A, B, C, D, E, i, J, dm, deta);

    std::vector<int> stations;
    for (int i = 1; i < f; ++i) stations.push_back(i);
    for (int i = l + 1; i < n_m - 1; ++i) stations.push_back(i);
    std::vector<PeriodicPDERow> rows;
    for (int i : stations) {
        const PeriodicRowKind kind = i == f - 1   ? PeriodicRowKind::LeTransition
                                     : i == l + 1 ? PeriodicRowKind::TeTransition
                                                  : PeriodicRowKind::Standard;
        const Weights& w = line[sz(i)];
        const double mixed = w.mixed;
        std::vector<PDETerm> terms;  // insertion order of v1's dict
        double rhs = 0.0;            // -Q / H, Q = 0
        auto add = [&](Node node, double weight) {
            for (PDETerm& t : terms)
                if (t.node == node) {
                    t.coefficient = t.coefficient + weight;
                    return;
                }
            terms.push_back(PDETerm{node, 0.0 + weight});
        };
        auto above = [&](int k, double weight) {
            const Elimination rule = ghost(k);
            if (!(rule.node == Node{k, n_eta}))
                throw DiscretizationError(std::format("the ghost rule for station {} is for ({}, {}), not ({}, {})", k,
                                                      rule.node.i, rule.node.j, k, n_eta));
            for (const Term& t : rule.terms) add(t.node, weight * t.coefficient);
            return -weight * rule.constant;
        };
        if (kind == PeriodicRowKind::TeTransition) {
            // psi_m, psi_m,m forward one-sided along j = J: raw weights C c / dm^2 + E c / (2 dm) over -H, and the
            // eta part of the center, 2 A / deta^2 over H
            const double a_ij = A(i, J), c_ij = C(i, J), e_ij = E(i, J), H = w.H;
            double along[4] = {0.0, 0.0, 0.0, 0.0};
            const double mm[4] = {2.0, -5.0, 4.0, -1.0};  // TE_MM_DIFFERENCE
            const double m1[3] = {-3.0, 4.0, -1.0};       // TE_M_DIFFERENCE
            for (int o = 0; o < 4; ++o) along[o] = along[o] + mm[o] * c_ij / (dm * dm);
            for (int o = 0; o < 3; ++o) along[o] = along[o] + m1[o] * e_ij / (2.0 * dm);
            add({i, J}, 2.0 * a_ij / (deta * deta) / H);
            for (int o = 0; o < 4; ++o) add({i + o, J}, -along[o] / H);
        } else {
            add({i, J}, 1.0);
            add({i - 1, J}, w.west);
            add({i + 1, J}, w.east);
        }
        add({i, J - 1}, w.south);
        rhs += above(i, w.north);
        static const std::vector<std::pair<int, double>> standard = {{-1, -1.0}, {1, 1.0}};
        static const std::vector<std::pair<int, double>> le = {{0, 3.0}, {-1, -4.0}, {-2, 1.0}};
        static const std::vector<std::pair<int, double>> te = {{0, -3.0}, {1, 4.0}, {2, -1.0}};
        const auto& diff = kind == PeriodicRowKind::Standard ? standard : kind == PeriodicRowKind::LeTransition ? le : te;
        for (const auto& [offset, c] : diff) {
            rhs += above(i + offset, c * mixed);
            add({i + offset, J - 1}, -c * mixed);
        }
        PeriodicPDERow row{kind, {i, J}, {}, rhs, w.H};
        for (const PDETerm& t : terms)
            if (t.coefficient != 0.0 || t.node == row.target) row.terms.push_back(t);
        // v1 PeriodicPDERow.__post_init__
        for (const PDETerm& t : row.terms)
            if (!std::isfinite(t.coefficient))
                throw DiscretizationError(
                    std::format("the periodic PDE row at ({}, {}) has a non-finite coefficient, rhs or H", i, J));
        if (!(std::isfinite(row.rhs) && std::isfinite(row.H) && row.H > 0.0))
            throw DiscretizationError(
                std::format("the periodic PDE row at ({}, {}) has a non-finite coefficient, rhs or H", i, J));
        if (kind != PeriodicRowKind::TeTransition && row.coefficient(row.target) != 1.0)
            throw DiscretizationError(
                std::format("the periodic PDE row at ({}, {}) must have the center coefficient 1", i, J));
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<PeriodicPDERow> build_periodic_pde_rows(const CoefficientFields& coefficients, const GridMetrics& metrics,
                                                    const BoundaryConditions& boundary_conditions)
{
    if (!metrics.grid) throw DiscretizationError("build_periodic_pde_rows takes GridMetrics with their grid");
    if (boundary_conditions.grid != metrics.grid)
        throw DiscretizationError("the boundary conditions were built from another grid than the metrics");
    const Grid& grid = *metrics.grid;
    check_shape(coefficients, grid);
    check_belongs(coefficients, metrics);
    const int n_eta = grid.n_eta();
    return periodic_pde_rows(coefficients.A, coefficients.B, coefficients.C, coefficients.D, coefficients.E,
                             coefficients.length_scale() * grid.dm(), grid.deta(), grid.first_blade_i,
                             grid.last_blade_i,
                             [&](int k) {  // the station's jump: 1 (v1), or the wake's with transpiration
                                 return boundary_conditions.periodic_ghost(k, n_eta, boundary_conditions.station_jump(k));
                             });
}

}  // namespace BladeGen::Turbolenza
