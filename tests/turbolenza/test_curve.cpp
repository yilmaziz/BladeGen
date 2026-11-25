// Curve.hpp: the clamped B-spline surface curves (port of v1 ParametricBladeCurve).
#include "RefData.hpp"

#include "turbolenza/Curve.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

namespace T = BladeGen::Turbolenza;

TEST_CASE("a quintic Bezier curve and its derivatives match the closed form", "[turbolenza][curve]")
{
    // Degree 5, no interior knot: a Bezier curve, checked against the Bernstein sum of its control
    // points and of its hodographs (forward differences) for every derivative order.
    std::vector<double> t(12, 0.0);
    for (int k = 6; k < 12; ++k) t[static_cast<std::size_t>(k)] = 1.0;
    std::vector<T::Point> P;
    const double ys[6] = {0.0, 0.3, -0.2, 0.7, 0.1, 0.5};
    for (int i = 0; i < 6; ++i) P.push_back({i / 5.0, ys[i]});
    const T::ParametricBladeCurve c(t, P, 1.0, 5);
    auto bern = [&](double u, int order) {
        // derivative of order `order` of the Bezier curve by the hodograph (forward differences)
        std::vector<T::Point> Q = P;
        double f = 1.0;
        for (int r = 0; r < order; ++r) {
            const int n = static_cast<int>(Q.size()) - 1;
            std::vector<T::Point> D;
            for (int i = 0; i < n; ++i) D.push_back({Q[static_cast<std::size_t>(i + 1)].m - Q[static_cast<std::size_t>(i)].m,
                                                     Q[static_cast<std::size_t>(i + 1)].y - Q[static_cast<std::size_t>(i)].y});
            f *= n;
            Q = D;
        }
        const int n = static_cast<int>(Q.size()) - 1;
        T::Point s{0, 0};
        for (int i = 0; i <= n; ++i) {
            double b = std::pow(u, i) * std::pow(1 - u, n - i);
            for (int a = 1; a <= i; ++a) b *= static_cast<double>(n - i + a) / a;
            s.m += b * Q[static_cast<std::size_t>(i)].m;
            s.y += b * Q[static_cast<std::size_t>(i)].y;
        }
        return T::Point{f * s.m, f * s.y};
    };
    for (double u : {0.0, 0.1, 0.37, 0.5, 0.93, 1.0})
        for (int k = 0; k <= 5; ++k) {
            const auto a = c.derivative(u, k), b = bern(u, k);
            CHECK(a.m == Catch::Approx(b.m).margin(1e-12));
            CHECK(a.y == Catch::Approx(b.y).margin(1e-12));
        }
    CHECK(c.evaluate(0.0).y == 0.0);
    CHECK(c.evaluate(1.0).y == 0.5);
    CHECK_THROWS_AS(c.derivative(1.1, 0), std::invalid_argument);
    CHECK_THROWS_AS(c.derivative(0.5, 6), std::invalid_argument);
    CHECK_THROWS_AS(T::ParametricBladeCurve(t, std::vector<T::Point>(P.begin(), P.end() - 1), 1.0, 5), std::invalid_argument);
}

TEST_CASE("the v1 C2 surface curves evaluate to v1's wall-node derivatives", "[turbolenza][curve][v1ref]")
{
    TURBOLENZA_REF_OR_SKIP(fam, "grid_family");
    const auto curves = ref::surface_curves(fam.at("cascade").at("active_curves"));
    double worst[3] = {0, 0, 0};
    for (const char* level : {"very_coarse", "coarse", "medium", "fine", "very_fine"}) {
        const auto grid = ref::grid(fam.at("levels").at(level));
        const double pitch = grid.pitch();
        for (int k = 0; k < 3; ++k) {
            // the scale of each order: the largest |X^(k)| over the curve's wall nodes
            double scale = 0.0;
            for (const auto* w : {&grid.wall_d[k].lower, &grid.wall_d[k].upper})
                for (const auto& p : *w) scale = std::max(scale, std::hypot(p.m, p.y));
            for (std::size_t b = 0; b < grid.wall_d[k].lower.size(); ++b) {
                const auto& wl = grid.wall_identity.lower;
                const auto& wu = grid.wall_identity.upper;
                const auto& cl = wl.surface[b] == T::Surface::Suction ? curves.suction : curves.pressure;
                const auto& cu = wu.surface[b] == T::Surface::Suction ? curves.suction : curves.pressure;
                const auto a = cl.derivative(wl.u[b], k);
                auto c = cu.derivative(wu.u[b], k);
                if (k == 0) c.y += pitch;
                const auto& lo = grid.wall_d[k].lower[b];
                const auto& up = grid.wall_d[k].upper[b];
                for (double d : {a.m - lo.m, a.y - lo.y, c.m - up.m, c.y - up.y})
                    worst[k] = std::max(worst[k], std::abs(d) / scale);
            }
        }
    }
    // measured (MSVC 19.51, 5 levels): order 0 exact, order 1 2.2e-13, order 2 2.5e-11 relative to the
    // order's largest value (round-off of the basis recurrence on the short nose spans, against
    // scipy's de Boor); the bounds cover the measurement with a margin of about 10
    INFO("max relative differences: " << worst[0] << ", " << worst[1] << ", " << worst[2]);
    CHECK(worst[0] <= 1e-15);
    CHECK(worst[1] <= 1e-12);
    CHECK(worst[2] <= 1e-10);
}
