// Tests of naca/Naca65: Table 4-1, the A10 mean line, the C4 thickness law, the monotone
// LE fair, both TE variants, quality, and the agreement with TURBOlenza's frozen C2 geometry.
#include "math/Constrained.hpp"
#include "naca/Naca65.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>

using namespace BladeGen;
using namespace BladeGen::Naca;
using Catch::Approx;
using Math::Vec2;

namespace {
double seg_dist(const Vec2& p, const Vec2& a, const Vec2& b)
{
    const Vec2 ab = b - a;
    const double t = std::clamp((p - a).dot(ab) / ab.squaredNorm(), 0.0, 1.0);
    return (p - (a + t * ab)).norm();
}
double poly_dist(const Vec2& p, const std::vector<Vec2>& poly)
{
    double d = 1e300;
    for (std::size_t i = 0; i + 1 < poly.size(); ++i) d = std::min(d, seg_dist(p, poly[i], poly[i + 1]));
    return d;
}
}  // namespace

TEST_CASE("Gauss-Legendre of any order and the NNLS / LSEI solvers", "[naca][math]")
{
    std::vector<double> x, w;
    Math::gauss_legendre(16, x, w);
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += w[i] * std::pow(x[i], 30);
    REQUIRE(s == Approx(2.0 / 31.0).epsilon(1e-13));
    // NNLS: min |A x - b|, x >= 0, with one active bound
    Eigen::MatrixXd A(3, 2);
    A << 1, 0, 0, 1, 1, 1;
    Eigen::VectorXd b(3);
    b << 2, -1, 1;
    const auto xn = Math::nnls(A, b);
    REQUIRE(xn[1] == Approx(0.0).margin(1e-14));
    REQUIRE(xn[0] == Approx(1.5).epsilon(1e-12));
    // LSEI: min |x - (1, 2, 3)| subject to x1 + x2 + x3 = 3, x >= 0
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3), E(1, 3);
    E << 1, 1, 1;
    Eigen::VectorXd t(3), f(1);
    t << 1, 2, 3;
    f << 3;
    const auto xl = Math::lsei(I, t, E, f);
    REQUIRE(xl);
    REQUIRE((*xl - Eigen::Vector3d(0, 1, 2)).norm() < 1e-12);
    f << -1;  // x >= 0 with a negative sum: infeasible
    REQUIRE_FALSE(Math::lsei(I, t, E, f));
}

TEST_CASE("Table 4-1 and the closed-form A10 mean line", "[naca]")
{
    const auto& T = table41();
    REQUIRE(T.thickness_printed[4] == Approx(0.03480));
    REQUIRE(T.thickness[4] == Approx(0.03148));  // the misprint, corrected
    // the tabulated Cl0 = 1 mean line agrees with the closed form to 2.3e-5 c
    double worst = 0.0;
    for (std::size_t i = 1; i + 1 < 26; ++i) worst = std::max(worst, std::abs(a10(T.x[i], 1.0)[0] - T.camber[i]));
    REQUIRE(worst < 2.5e-5);
    // derivatives against finite differences
    const double x = 0.3, h = 1e-6;
    const auto d = a10(x, 1.2);
    REQUIRE(d[1] == Approx((a10(x + h, 1.2)[0] - a10(x - h, 1.2)[0]) / (2 * h)).epsilon(1e-8));
    REQUIRE(d[2] == Approx((a10(x + h, 1.2)[1] - a10(x - h, 1.2)[1]) / (2 * h)).epsilon(1e-7));
    REQUIRE_THROWS_AS(a10(0.0, 1.0), std::invalid_argument);
    const auto ang = equivalent_arc_angles(1.2);
    REQUIRE(std::tan(ang[0] / 4.0) == Approx(0.1103 * 1.2));
}

TEST_CASE("thickness law: exact at every station (scaled), C4 quintic, maximum as TURBOlenza", "[naca]")
{
    const Naca65Section s({1.2, 0.08, Variant::EmeryBase});
    const auto& T = table41();
    for (std::size_t i = 0; i < 26; ++i)
        REQUIRE(s.table_half_thickness(T.x[i]) == Approx(0.5 * T.thickness[i] * 0.8).margin(1e-15));
    // quintic_natural_sqrt maximum: 0.99961 t at station 0.4098 (TURBOlenza docstring)
    REQUIRE(s.metadata().max_thickness / 0.08 == Approx(0.99961).epsilon(2e-5));
    REQUIRE(s.metadata().max_thickness_station == Approx(0.4098).margin(5e-4));
    // station: 0 at the LE, 1 at the TE, monotone
    REQUIRE(s.station(0.0) == Approx(0.0).margin(1e-15));
    REQUIRE(s.station(1.0) == Approx(1.0).epsilon(1e-14));
    REQUIRE(s.x_of_station(s.station(0.37)) == Approx(0.37).epsilon(1e-12));
}

TEST_CASE("the C2 leading edge: nose circle, monotone curvature, G2 onto the surface", "[naca]")
{
    for (auto var : {Variant::KovachSandercock, Variant::EmeryBase}) {
        const Naca65Section s({1.2, 0.10, var});
        const auto& m = s.metadata();
        REQUIRE(m.le_radius == Approx(0.00687));
        REQUIRE(m.le_center.norm() == Approx(m.le_radius).epsilon(1e-14));
        for (auto side : {Profile::Side::Suction, Profile::Side::Pressure}) {
            const auto c = s.sample_side(side, 400, 6000);  // the first body sample ~1e-7 c past the anchor
            REQUIRE(c.xy.front().norm() < 1e-15);
            // curvature at the apex is 1/R_LE and falls monotonically along the fair
            REQUIRE(c.kappa.front() == Approx(1.0 / m.le_radius).epsilon(1e-12));
            for (std::size_t i = 1; i < c.kappa.size() && c.piece[i] == 0; ++i)
                REQUIRE(c.kappa[i] <= c.kappa[i - 1] * (1.0 + 1e-12));
            // G2 at the anchor: the fair's end curvature equals the surface's there
            std::size_t j = 0;
            while (c.piece[j] == 0) ++j;
            REQUIRE((c.xy[j] - c.xy[j - 1]).norm() < 2e-3);
            REQUIRE(c.kappa[j - 1] == Approx(c.kappa[j]).epsilon(2e-4));
        }
        REQUIRE(m.le_fairness[0] < 4.0);
        REQUIRE(m.le_fairness[1] < 4.0);
    }
}

TEST_CASE("Kovach TE: circle through (1, 0) on the camber line, G2 fair with its curvature", "[naca]")
{
    const Naca65Section s({1.2, 0.10, Variant::KovachSandercock});
    const auto& m = s.metadata();
    REQUIRE(m.te_radius == Approx(0.008));
    REQUIRE((m.te_center - Vec2(1, 0)).norm() == Approx(0.008).epsilon(1e-12));
    REQUIRE(std::abs(m.te_center.y() - s.camber_point(m.te_center.x()).y()) < 1e-15);
    const auto ss = s.sample_side(Profile::Side::Suction, 400, 200), ps = s.sample_side(Profile::Side::Pressure, 400, 200);
    REQUIRE((ss.xy.back() - Vec2(1, 0)).norm() < 1e-15);
    REQUIRE((ps.xy.back() - Vec2(1, 0)).norm() < 1e-15);
    REQUIRE(ss.kappa.back() == Approx(1.0 / 0.008).epsilon(1e-6));
    REQUIRE(ps.kappa.back() == Approx(1.0 / 0.008).epsilon(1e-6));
    INFO("TE fair deviation from the Kovach circle " << m.te_circle_deviation);
    REQUIRE(m.te_circle_deviation < 2e-4);
}

TEST_CASE("Emery TE: sharp, the closing piece C2 with the surface", "[naca]")
{
    const Naca65Section s({1.2, 0.10, Variant::EmeryBase});
    REQUIRE(s.sharp_te());
    const auto ss = s.sample_side(Profile::Side::Suction, 20000, 400);  // TE step 2.5e-7 c
    REQUIRE((ss.xy.back() - Vec2(1, 0)).norm() < 1e-15);
    std::size_t j = 0;
    while (ss.piece[j] != 2) ++j;
    REQUIRE(ss.kappa[j] == Approx(ss.kappa[j - 1]).epsilon(2e-3));
    const auto q = check_quality(s);
    INFO(q.summary());
    REQUIRE(q.ok());
}

TEST_CASE("NACA 65 sections of the reference range have no geometry errors", "[naca]")
{
    for (double cl0 : {0.0, 0.4, 1.2, 1.8})
        for (double t : {0.06, 0.10, 0.12})
            for (auto var : {Variant::KovachSandercock, Variant::EmeryBase}) {
                const Naca65Section s({cl0, t, var});
                const auto q = check_quality(s);
                INFO(s.params().designation() << " " << to_string(var) << ": " << q.summary());
                REQUIRE(q.ok());
            }
    REQUIRE_THROWS_AS(Naca65Section({3.0, 0.1, Variant::EmeryBase}), std::invalid_argument);
    REQUIRE_THROWS_AS(Naca65Section({1.0, 0.5, Variant::EmeryBase}), std::invalid_argument);
    REQUIRE_THROWS_AS(variant_from_string("dca"), std::invalid_argument);
    REQUIRE(Naca65Params{1.2, 0.10, Variant::KovachSandercock}.designation() == "NACA 65-(12A10)10");
}

TEST_CASE("BladeGen NACA 65 = TURBOlenza frozen C2 geometry (reference cases)", "[naca][reference]")
{
    std::ifstream f(std::string(BLADEGEN_SOURCE_DIR) + "/tests/data/turbolenza_naca65_reference.csv");
    REQUIRE(f.good());
    struct Case { double cl0 = 0, t = 0; std::string var; std::vector<Vec2> ss, ps; };
    std::map<int, Case> cases;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("case", 0) == 0) continue;
        std::stringstream st(line);
        std::string c[7];
        for (auto& v : c) std::getline(st, v, ',');
        auto& k = cases[std::stoi(c[0])];
        k.cl0 = std::stod(c[1]);
        k.t = std::stod(c[2]);
        k.var = c[3];
        (c[4] == "SS" ? k.ss : k.ps).emplace_back(std::stod(c[5]), std::stod(c[6]));
    }
    REQUIRE(cases.size() == 5);
    for (const auto& [id, k] : cases) {
        const Naca65Section s({k.cl0, k.t, variant_from_string(k.var)});
        for (int sd = 0; sd < 2; ++sd) {
            const auto mine = s.sample_side(sd == 0 ? Profile::Side::Suction : Profile::Side::Pressure, 1500, 3000).xy;
            double body = 0.0, le = 0.0, te = 0.0;
            for (const auto& p : (sd == 0 ? k.ss : k.ps)) {
                const double d = poly_dist(p, mine);
                if (p.x() < 0.04) le = std::max(le, d);
                else if (p.x() > 0.97) te = std::max(te, d);
                else body = std::max(body, d);
            }
            INFO("case " << id << " " << k.var << " cl0 " << k.cl0 << " t " << k.t << (sd == 0 ? " SS" : " PS")
                         << ": body " << body << ", LE zone " << le << ", TE zone " << te);
            const bool emery = k.var == "emery_base";
            // the same table, mean line, station, thickness law and LE fair rule: Emery agrees to 5e-8 c.
            // Kovach differs by design at the taper start (C2 blend vs TURBOlenza's kernel smoothing,
            // <= 3e-6 c) and at the TE (BladeGen's G2 fair connector vs the smoothed radial fair)
            REQUIRE(body < (emery ? 2e-7 : 5e-6));
            REQUIRE(le < 1e-7);
            REQUIRE(te < (emery ? 2e-5 : 2.5e-4));
        }
    }
}

#include "fit/NacaFit.hpp"

TEST_CASE("NACA 65 + Custom: custom parameters fitted to the NACA section", "[naca][fit][slow]")
{
    for (const Naca65Params p : {Naca65Params{1.2, 0.10, Variant::KovachSandercock}, Naca65Params{0.4, 0.06, Variant::KovachSandercock},
                                 Naca65Params{1.8, 0.12, Variant::EmeryBase}}) {
        const Naca65Section s(p);
        const auto r = Fit::fit_naca65(s);
        const Profile::Profile2D prof(r.params);
        const auto q = Profile::check_quality(prof);
        INFO(p.designation() << " " << to_string(p.variant) << ": rms " << r.rms << " max " << r.max << " in " << r.seconds
                             << " s, " << q.summary());
        REQUIRE(r.rms < (p.variant == Variant::EmeryBase ? 1e-3 : 3e-4));  // a round TE cannot be sharp
        REQUIRE(q.ok());
    }
}
