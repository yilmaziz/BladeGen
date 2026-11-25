// Convergence (port of v1 convergence.py): the change norms in closed form and the stop tests.
#include "turbolenza/Convergence.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <limits>

using namespace BladeGen::Turbolenza;

TEST_CASE("convergence: change norms and stop tests", "[turbolenza][convergence]")
{
    Field a(1, 3), b(1, 3);
    a.v = {1.0, 2.0, 4.0};
    b.v = {1.5, 2.0, 3.0};
    CHECK(change_norm(b, a) == 1.0 / 4.0);  // max |b - a| / max |a|
    Field u(1, 2), v(1, 2), U(1, 2), V(1, 2);
    u.v = {3.0, 0.0};
    v.v = {4.0, 1.0};
    U.v = {0.0, 0.0};
    V.v = {0.0, 1.0};
    CHECK(velocity_change(U, V, u, v) == 5.0 / 5.0);  // |(3, 4)| / max(5, 1)
    CHECK_THROWS_AS(change_norm(a, Field(1, 3, 0.0)), std::invalid_argument);
    CHECK_THROWS_AS(change_norm(a, Field(1, 2, 1.0)), std::invalid_argument);
    Field n = a;
    n.v[1] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(change_norm(n, a), std::invalid_argument);

    IterationSettings s;
    s.density_tolerance = 1e-5;
    s.kutta_tolerance = 1e-4;
    History h;
    CHECK(!density_converged(h, s));
    h.density_change = {1e-6};
    h.velocity_change = {std::numeric_limits<double>::quiet_NaN()};  // pass 0: no previous velocity, never converged
    CHECK(!density_converged(h, s));
    h.density_change.push_back(1e-6);
    h.velocity_change.push_back(2e-5);
    CHECK(!density_converged(h, s));
    s.velocity_tolerance = 3e-5;
    CHECK(density_converged(h, s));
    CHECK(!coupled_converged(h, s));
    h.kutta_error = {1e-3, std::numeric_limits<double>::quiet_NaN()};
    CHECK(!kutta_converged(h, s));
    h.kutta_error.back() = 1e-4;
    CHECK(coupled_converged(h, s));
}
