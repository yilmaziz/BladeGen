// Discretization (port of v1 discretization.py, expanded rows): the nine-point Eq. 5-44 and the
// periodic PDE rows reproduce the exact derivatives of quadratics (the TE row's psi_m,m of a
// cubic), the row sums, the stencil of v1 pass 0 through its periodic rows, the refusals.  The
// assembled v1 matrix (test_matrix_assembly) checks every weight against v1 entry by entry.
#include "SolverRef.hpp"

#include "turbolenza/Discretization.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace BladeGen::Turbolenza;

namespace {

struct Quadratic
{
    double a, b, c, d, e, f, g;  // psi = a m^2 + b m eta + c eta^2 + d m + e eta + f + g m^3
    double operator()(double m, double eta) const
    {
        return a * m * m + b * m * eta + c * eta * eta + d * m + e * eta + f + g * m * m * m;
    }
};

}  // namespace

TEST_CASE("discretization: Eq. 5-44 is exact for quadratics", "[turbolenza][discretization]")
{
    const int n_m = 9, n_eta = 6;
    const double dm = 0.13, deta = 0.2;
    const double A = 3.1, B = -0.7, C = 0.9, D = 4.0, E = -2.5;
    const InteriorStencil s = interior_stencil(Field(n_m, n_eta, A), Field(n_m, n_eta, B), Field(n_m, n_eta, C),
                                               Field(n_m, n_eta, D), Field(n_m, n_eta, E), dm, deta);
    const Quadratic q{0.8, -1.3, 2.2, 0.4, -0.6, 0.1, 0.0};
    for (int i = 1; i < n_m - 1; ++i)
        for (int j = 1; j < n_eta - 1; ++j) {
            const double m = i * dm, eta = j * deta;
            auto psi = [&](int di, int dj) { return q((i + di) * dm, (j + dj) * deta); };
            const double sum = s.center(i, j) * psi(0, 0) + s.west(i, j) * psi(-1, 0) + s.east(i, j) * psi(1, 0) +
                               s.south(i, j) * psi(0, -1) + s.north(i, j) * psi(0, 1) +
                               s.mixed(i, j) * (psi(1, 1) - psi(1, -1) - psi(-1, 1) + psi(-1, -1));
            const double pde = A * 2.0 * q.c - 2.0 * B * q.b + C * 2.0 * q.a + D * (q.b * m + 2.0 * q.c * eta + q.e) +
                               E * (2.0 * q.a * m + q.b * eta + q.d);
            CHECK(std::abs(sum - (-pde / s.H(i, j))) <= 1e-13);
            // the four neighbours sum to -1: a constant psi solves the row
            CHECK(std::abs(s.west(i, j) + s.east(i, j) + s.south(i, j) + s.north(i, j) + 1.0) <= 1e-15);
            CHECK(s.H(i, j) == 2.0 * A / (deta * deta) + 2.0 * C / (dm * dm));
        }
    CHECK(std::isnan(s.center(0, 2)));
    CHECK(std::isnan(s.west(3, 0)));
    CHECK(std::isnan(s.mixed(n_m - 1, 2)));
}

TEST_CASE("discretization: the periodic PDE rows are exact on quasi-periodic quadratics", "[turbolenza][discretization]")
{
    const int n_m = 16, n_eta = 5, f = 4, l = 9;
    const double dm = 0.1, deta = 0.25;
    const double A = 2.0, B = 0.6, C = 1.4, D = -3.0, E = 1.7;
    auto ghost = [&](int k) { return Elimination{{k, n_eta}, {Term{{k, 1}, 1.0}}, 1.0}; };
    const std::vector<PeriodicPDERow> rows =
        periodic_pde_rows(Field(n_m, n_eta, A), Field(n_m, n_eta, B), Field(n_m, n_eta, C), Field(n_m, n_eta, D),
                          Field(n_m, n_eta, E), dm, deta, f, l, ghost);
    REQUIRE(rows.size() == static_cast<std::size_t>((f - 1) + (n_m - 2 - l)));
    // psi = eta + a m^2 + c m (+ g m^3 for the TE row's psi_m,m): psi[k, J + 1] = psi[k, 1] + 1 exactly
    for (const double g3 : {0.0, 0.9}) {
        const Quadratic q{0.7, 0.0, 0.0, -0.4, 1.0, 0.0, g3};
        Field psi(n_m, n_eta);
        for (int i = 0; i < n_m; ++i)
            for (int j = 0; j < n_eta; ++j) psi(i, j) = q(i * dm, j * deta);
        for (const PeriodicPDERow& r : rows) {
            if (g3 != 0.0 && r.kind != PeriodicRowKind::TeTransition) continue;  // only the TE row is exact for cubics
            const double m = r.target.i * dm;
            const double psi_mm = 2.0 * q.a + 6.0 * g3 * m;
            // the TE row's forward (-3, 4, -1) / (2 dm) is exact for quadratics; on a cubic its
            // error is exactly -dm^2 f(3) / 3 = -2 g dm^2 (the Taylor remainder, no higher terms)
            const double psi_m = 2.0 * q.a * m + q.d + 3.0 * g3 * m * m - 2.0 * g3 * dm * dm;
            const double pde = C * psi_mm + D * 1.0 + E * psi_m;
            CHECK(std::abs(r.residual(psi) - (-pde / r.H)) <= 1e-13);
        }
    }
    CHECK(rows.front().kind == PeriodicRowKind::Standard);
    CHECK(rows[static_cast<std::size_t>(f - 2)].kind == PeriodicRowKind::LeTransition);
    CHECK(rows[static_cast<std::size_t>(f - 1)].kind == PeriodicRowKind::TeTransition);
    for (const PeriodicPDERow& r : rows) {
        if (r.kind != PeriodicRowKind::TeTransition) CHECK(r.coefficient(r.target) == 1.0);
        // the north weight times the jump 1 is the row's right-hand side (the mixed jumps cancel)
        CHECK(std::abs(r.rhs - (A / (deta * deta) + D / (2.0 * deta)) / r.H) <= 1e-15);
    }
    // the TE row's center (2 A / deta^2 - 2 C / dm^2 + 3 E / (2 dm)) / H
    const PeriodicPDERow& te = rows[static_cast<std::size_t>(f - 1)];
    CHECK(std::abs(te.coefficient(te.target) -
                   (2.0 * A / (deta * deta) - 2.0 * C / (dm * dm) + 3.0 * E / (2.0 * dm)) / te.H) <= 1e-15);
}

TEST_CASE("discretization: refusals", "[turbolenza][discretization]")
{
    const Field one(9, 6, 1.0);
    CHECK_THROWS_AS(interior_stencil(Field(9, 6, -1.0), one, Field(9, 6, -1.0), one, one, 0.1, 0.2), DiscretizationError);
    CHECK_THROWS_AS(interior_stencil(one, one, one, one, one, 0.0, 0.2), DiscretizationError);
    CHECK_THROWS_AS(interior_stencil(one, one, one, one, Field(9, 5, 1.0), 0.1, 0.2), DiscretizationError);
    Field bad(9, 6, 1.0);
    bad(3, 3) = NAN;
    CHECK_THROWS_AS(interior_stencil(one, bad, one, one, one, 0.1, 0.2), DiscretizationError);
    auto ghost = [](int k) { return Elimination{{k, 6}, {Term{{k, 1}, 1.0}}, 1.0}; };
    // three inlet and four exit periodic stations are needed
    CHECK_THROWS_AS(periodic_pde_rows(one, one, one, one, one, 0.1, 0.2, 2, 4, ghost), DiscretizationError);
    CHECK_THROWS_AS(periodic_pde_rows(one, one, one, one, one, 0.1, 0.2, 3, 5, ghost), DiscretizationError);
    // coefficients of another grid's metrics are refused
    const auto g = std::make_shared<const Grid>(sref::channel(31, 7, 8, 20, 0.2));
    const GridMetrics m = build_metrics(g);
    CoefficientFields c = build_coefficients(m, Field(31, 7, 1.0), 1.0, 0.1);
    c.c_ref = 0.2;
    CHECK_THROWS_AS(build_interior_stencil(c, m), DiscretizationError);
}
