// The general curvilinear multi-block operator (GeneralOperator.hpp), verification stages 1 and 2 of the
// plan: the nine-point table re-derived, the metric identities (analytic and discrete), the discrete
// geometric conservation law and the uniform flow, the conservation of the face fluxes across interfaces,
// manufactured solutions on six mesh families with the observed order of accuracy, and r, b, omega in the
// operator and the velocity recovery.  Every measured number is printed as a MEASURE line.
#include "GeneralFixture.hpp"

#include "turbolenza/GeneralOperator.hpp"
#include "turbolenza/GeneralSolver.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <format>
#include <map>

using namespace BladeGen::Turbolenza;
using gfix::measure;

namespace {

std::size_t z(int k) { return static_cast<std::size_t>(k); }

/// The flux of face f of node P in its +axis direction at psi.
double face_flux(const GeneralOperator& op, int P, const FaceGeometry& f, const std::vector<double>& psi,
                 const std::vector<double>& c)
{
    const double kf = 0.5 * (c[z(P)] + c[z(f.q.node)]);
    const double dpsi = static_cast<double>(f.sign) * (psi[z(f.q.node)] - f.q.images - psi[z(P)]);
    return kf * (f.G_nn * dpsi + f.G_nt * op.apply(f.psi_t, psi, 1.0));
}

bool all_central(const GeneralOperator& op, int g)
{
    auto central = [&](int n) { return op.derivative(n, 0).size() == 2 && op.derivative(n, 1).size() == 2; };
    if (!central(g)) return false;
    for (const FaceGeometry& f : op.faces(g))
        if (!central(f.q.node)) return false;
    return !op.faces(g).empty();
}

}  // namespace

TEST_CASE("general: the nine-point table of the flux form, derived again", "[turbolenza][general]")
{
    // A curvilinear block, variable coefficient: the row of an interior node read from the residual of
    // unit fields against the table P: -(pE + pW + pN + pS), E: pE + cN - cS, W: pW - cN + cS,
    // N: pN + cE - cW, S: pS - cE + cW, NE: cE + cN, NW: -cW - cN, SE: -cE - cS, SW: cW + cS,
    // p = K_nn of the face, c = K_nt / 4 (h = k = 1).
    const GeneralOperator op(gfix::dirichlet_box(gfix::block(gfix::family(5), 17, 15)));
    const gfix::Manufactured M = gfix::smooth(false);
    OperatorFields fields = gfix::mms_fields(op, M);
    fields.source.clear();
    const int i = 7, j = 6;
    const int P = op.node(0, i, j);
    REQUIRE(op.equation(P) == NodeEquation::Pde);
    std::map<std::pair<int, int>, double> p, c;  // by (axis, sign)
    for (const FaceGeometry& f : op.faces(P)) {
        const double kf = 0.5 * (fields.coefficient[z(P)] + fields.coefficient[z(f.q.node)]);
        p[{f.axis, f.sign}] = kf * f.G_nn;
        c[{f.axis, f.sign}] = kf * f.G_nt / 4.0;
    }
    const double pW = p[{0, -1}], pE = p[{0, 1}], pS = p[{1, -1}], pN = p[{1, 1}];
    const double cW = c[{0, -1}], cE = c[{0, 1}], cS = c[{1, -1}], cN = c[{1, 1}];
    const std::map<std::pair<int, int>, double> table = {
        {{0, 0}, -(pE + pW + pN + pS)}, {{1, 0}, pE + cN - cS},   {{-1, 0}, pW - cN + cS},
        {{0, 1}, pN + cE - cW},         {{0, -1}, pS - cE + cW},  {{1, 1}, cE + cN},
        {{-1, 1}, -cW - cN},            {{1, -1}, -cE - cS},      {{-1, -1}, cW + cS}};
    double worst = 0.0, scale = 0.0;
    for (const auto& [d, want] : table) {
        std::vector<double> unit(z(op.n_nodes()), 0.0);
        unit[z(op.node(0, i + d.first, j + d.second))] = 1.0;
        const double got = op.residual(unit, fields)[z(P)];
        worst = std::max(worst, std::abs(got - want));
        scale = std::max(scale, std::abs(want));
    }
    // and nothing else in the row
    double others = 0.0;
    for (int g = 0; g < op.n_nodes(); ++g) {
        const auto l = op.location(g);
        if (std::abs(l[1] - i) <= 1 && std::abs(l[2] - j) <= 1) continue;
        std::vector<double> unit(z(op.n_nodes()), 0.0);
        unit[z(g)] = 1.0;
        others = std::max(others, std::abs(op.residual(unit, fields)[z(P)]));
    }
    CHECK(measure("general nine-point table vs assembled row (relative)", worst / scale) <= 1e-14);  // measured ~1e-16
    CHECK(others == 0.0);
    // the cross coefficients are not zero on this mesh (the table is exercised)
    CHECK(std::abs(cE) > 1e-3 * pE);
}

TEST_CASE("general: metric identities of Aungier's mapping (analytic and discrete)", "[turbolenza][general]")
{
    // x = m'(m) = ln(1 + 0.3 m) / 0.3 (r = 1 + 0.3 m), y = theta0(m) + eta Theta(m), (m, eta) computational:
    // J = S / r^2 and J g = [[S, -T], [-T, (1 + T^2) / S]], T = r dtheta/dm|eta, S = r Theta (Eqs. 5-3, 5-9).
    auto r = [](double m) { return 1.0 + 0.3 * m; };
    auto th0 = [](double m) { return 0.2 * std::sin(m); };
    auto dth0 = [](double m) { return 0.2 * std::cos(m); };
    auto Th = [](double m) { return 0.5 + 0.1 * m; };
    double worst = 0.0;
    for (double m : {0.0, 0.37, 1.1, 2.0})
        for (double eta : {0.0, 0.3, 1.0}) {
            const double x_m = 1.0 / r(m), y_m = dth0(m) + eta * 0.1, y_e = Th(m);
            const double T = r(m) * y_m, S = r(m) * Th(m);
            const auto K = contravariant_tensor(x_m, 0.0, y_m, y_e, 1.0);
            const double J = x_m * y_e;
            worst = std::max({worst, std::abs(J - S / (r(m) * r(m))) / J, std::abs(K[0] - S) / S, std::abs(K[1] + T) / std::max(1.0, std::abs(T)),
                              std::abs(K[2] - (1.0 + T * T) / S) / K[2]});
        }
    CHECK(measure("general Aungier map: J = S/r^2, Jg = [[S,-T],[-T,(1+T^2)/S]] (rel)", worst) <= 1e-15);  // measured 1.1e-16

    // discrete: the nodal metrics of the mesh of that map against the analytic J, second order
    double previous = 0.0;
    for (int N : {16, 32, 64, 128}) {
        const double m1 = 2.0;
        GeneralMesh mesh = gfix::dirichlet_box(gfix::block(
            [&](double u, double v) {
                const double m = m1 * u;
                return std::array<double, 2>{std::log(r(m)) / 0.3, th0(m) + v * Th(m)};
            },
            N + 1, N + 1));
        const GeneralOperator op(mesh);
        double e = 0.0;
        for (int g = 0; g < op.n_nodes(); ++g) {
            const auto l = op.location(g);
            const double m = m1 * l[1] / N;
            const double J = op.metric(g).J / ((m1 / N) * (1.0 / N));  // per unit (m, eta)
            e = std::max(e, std::abs(J - Th(m) / r(m)) / (Th(m) / r(m)));
        }
        measure(std::format("general Aungier map discrete J error N={}", N), e);
        if (previous > 0.0) {
            const double p = measure(std::format("general Aungier map discrete J order N={}", N), gfix::order(previous, e));
            CHECK(p > 1.9);  // measured 2.0
        }
        previous = e;
    }
}

TEST_CASE("general: discrete metric identities and the uniform flow", "[turbolenza][general]")
{
    for (int fam : {4, 5}) {
        const GeneralOperator op(gfix::dirichlet_box(gfix::block(gfix::family(fam), 25, 21)));
        const auto gcl = op.gcl_residual();
        double inner = 0.0, edge = 0.0;
        for (int g = 0; g < op.n_nodes(); ++g) {
            const double v = std::max(std::abs(gcl[z(g)][0]), std::abs(gcl[z(g)][1]));
            (all_central(op, g) ? inner : edge) = std::max(all_central(op, g) ? inner : edge, v);
        }
        CHECK(measure(std::format("general GCL residual, all-central nodes, {}", gfix::family_name(fam)), inner) <= 1e-15);  // measured ~1e-17
        measure(std::format("general GCL residual, nodes next to one-sided differences, {}", gfix::family_name(fam)), edge);

        // psi = a x + b y with a constant coefficient: zero residual wherever the GCL holds
        OperatorFields f;
        f.coefficient.assign(z(op.n_nodes()), 1.3);
        f.dirichlet = [&](int g) { return 0.7 * op.x(g) - 1.9 * op.y(g); };
        std::vector<double> psi(z(op.n_nodes()));
        for (int g = 0; g < op.n_nodes(); ++g) psi[z(g)] = f.dirichlet(g);
        const auto res = op.residual(psi, f);
        double rin = 0.0, rout = 0.0;
        for (int g = 0; g < op.n_nodes(); ++g) {
            if (op.equation(g) != NodeEquation::Pde) continue;
            (all_central(op, g) ? rin : rout) = std::max(all_central(op, g) ? rin : rout, std::abs(res[z(g)]));
            // the residual of a uniform flow is the GCL: a R_y - b R_x ... (F = c (a y_t - b x_t) sigma)
            const double expect = 1.3 * (0.7 * gcl[z(g)][0] + (-1.9) * gcl[z(g)][1]);
            CHECK(std::abs(res[z(g)] - expect) <= 1e-13);
        }
        CHECK(measure(std::format("general uniform-flow residual, all-central nodes, {}", gfix::family_name(fam)), rin) <= 1e-13);
        measure(std::format("general uniform-flow residual, other PDE nodes, {}", gfix::family_name(fam)), rout);
        const auto psi_h = gfix::mms_solve(op, f);
        double err = 0.0;
        for (int g = 0; g < op.n_nodes(); ++g) err = std::max(err, std::abs(psi_h[z(g)] - psi[z(g)]));
        measure(std::format("general uniform flow solved, max error, {}", gfix::family_name(fam)), err);
    }
}

TEST_CASE("general: manufactured solutions on single blocks, observed order", "[turbolenza][general][mms]")
{
    for (bool constant : {true, false})
        for (int fam = 1; fam <= 5; ++fam) {
            const gfix::Manufactured M = gfix::smooth(constant);
            std::vector<gfix::MmsError> e;
            for (int N : {16, 32, 64, 128}) {
                const GeneralOperator op(gfix::dirichlet_box(gfix::block(gfix::family(fam), N + 1, N + 1)));
                const OperatorFields f = gfix::mms_fields(op, M);
                e.push_back(gfix::mms_error(op, M, gfix::mms_solve(op, f)));
            }
            const std::string name = std::format("{}, {} coefficient", gfix::family_name(fam), constant ? "constant" : "variable");
            for (std::size_t k = 0; k < e.size(); ++k)
                measure(std::format("general mms psi error N={} {}", 16 << k, name), e[k].psi);
            const double p_psi = measure("general mms psi order 64->128 " + name, gfix::order(e[2].psi, e[3].psi));
            const double p_grad = measure("general mms grad order 64->128 " + name, gfix::order(e[2].grad, e[3].grad));
            CHECK(p_psi > 1.85);   // measured 2.0 on every family
            CHECK(p_grad > 1.5);   // measured ~2 (one-sided at the boundary)
        }
}

TEST_CASE("general: manufactured solutions on three turned blocks with periodic and interface faces",
          "[turbolenza][general][mms]")
{
    const double P = 0.8;
    const gfix::Manufactured M = gfix::quasi_periodic(P);
    for (int variant = 0; variant < 2; ++variant) {
        // 0: Dirichlet inlet and outlet; 1: Dirichlet inlet, the flow-angle condition at the outlet
        const SegmentKind outlet = variant == 0 ? SegmentKind::Dirichlet : SegmentKind::Outlet;
        std::vector<gfix::MmsError> e;
        for (int Mc : {8, 16, 32, 64}) {
            const GeneralOperator op(gfix::three_block_periodic(Mc, P, SegmentKind::Dirichlet, outlet));
            const OperatorFields f = gfix::mms_fields(op, M);
            const auto psi = gfix::mms_solve(op, f);
            e.push_back(gfix::mms_error(op, M, psi));
            // conservation: every face between two PDE nodes enters one control volume and leaves the
            // other (across interfaces and the periodic cut, where a control volume is split between
            // blocks), so the sum of all class residuals is the flux through the faces to the boundary
            if (Mc == 16) {
                OperatorFields nf = f;
                nf.source.clear();
                const auto res = op.residual(psi, nf);
                double total = 0.0, boundary = 0.0, scale = 0.0;
                for (int g = 0; g < op.n_nodes(); ++g) {
                    if (op.equation(g) == NodeEquation::Pde && op.representative(g) == g) total += res[z(g)];
                    for (const FaceGeometry& fa : op.faces(g)) {
                        const double F = fa.sign * fa.weight * face_flux(op, g, fa, psi, f.coefficient);
                        scale += std::abs(F);
                        if (op.equation(fa.q.node) != NodeEquation::Pde) boundary += F;
                    }
                }
                CHECK(measure("general telescoping of the face fluxes (relative)", std::abs(total - boundary) / scale) <= 1e-13);
                CHECK(scale > 0.0);
            }
        }
        const std::string name = variant == 0 ? "3 blocks, Dirichlet ends" : "3 blocks, flow-angle outlet";
        for (std::size_t k = 0; k < e.size(); ++k) measure(std::format("general mms psi error M={} {}", 8 << k, name), e[k].psi);
        const double p = measure("general mms psi order 32->64 " + name, gfix::order(e[2].psi, e[3].psi));
        const double pg = measure("general mms grad order 32->64 " + name, gfix::order(e[2].grad, e[3].grad));
        CHECK(p > 1.85);
        CHECK(pg > 1.5);
    }
}

TEST_CASE("general: manufactured solution on an O-block round a body", "[turbolenza][general][mms]")
{
    const gfix::Manufactured M = gfix::smooth(false);
    std::vector<gfix::MmsError> e;
    for (int N : {16, 32, 64, 128}) {
        const GeneralOperator op(gfix::oblock(2 * N + 1, N + 1));
        e.push_back(gfix::mms_error(op, M, gfix::mms_solve(op, gfix::mms_fields(op, M))));
    }
    for (std::size_t k = 0; k < e.size(); ++k) measure(std::format("general mms psi error O-block N={}", 16 << k), e[k].psi);
    const double p = measure("general mms psi order 64->128 O-block", gfix::order(e[2].psi, e[3].psi));
    const double pg = measure("general mms grad order 64->128 O-block", gfix::order(e[2].grad, e[3].grad));
    CHECK(p > 1.85);
    CHECK(pg > 1.5);
}

TEST_CASE("general: r, b and omega in the operator and the velocity recovery", "[turbolenza][general][mms]")
{
    // mdot' div(kappa grad psi) = omega (r^2)_x + G, kappa = 1 / (rho b), with r = 1 + 0.5 x, b = 1 + 0.2 x,
    // rho = 1.2 + 0.1 sin(x) cos(y), mdot' = 1.7, omega = 3: psi, W_m = mdot' kappa psi_y / r and
    // W_theta = -mdot' kappa psi_x / r against the exact ones.
    const double mdot = 1.7, omega = 3.0;
    auto r = [](double x, double) { return 1.0 + 0.5 * x; };
    auto b = [](double x, double) { return 1.0 + 0.2 * x; };
    auto rho = [](double x, double y) { return 1.2 + 0.1 * std::sin(x) * std::cos(y); };
    auto rho_x = [](double x, double y) { return 0.1 * std::cos(x) * std::cos(y); };
    auto rho_y = [](double x, double y) { return -0.1 * std::sin(x) * std::sin(y); };
    gfix::Manufactured M = gfix::smooth(true);
    M.c = [=](double x, double y) { return mdot / (rho(x, y) * b(x, y)); };
    M.c_x = [=](double x, double y) {
        const double rb = rho(x, y) * b(x, y);
        return -mdot * (rho_x(x, y) * b(x, y) + rho(x, y) * 0.2) / (rb * rb);
    };
    M.c_y = [=](double x, double y) { return -mdot * rho_y(x, y) / (rho(x, y) * rho(x, y) * b(x, y)); };
    std::vector<double> ep, ew;
    for (int N : {16, 32, 64, 128}) {
        const GeneralOperator op(gfix::dirichlet_box(gfix::block(gfix::family(5), N + 1, N + 1, r, b)));
        const std::vector<double> rot = op.rotation_source(omega);
        // G = mdot' div(kappa grad psi) - omega (r^2)_x; the source is (J f) = rot + J G, so the solve's
        // source is J div(C grad psi) - J omega (r^2)_x + rot
        std::vector<double> extra(z(op.n_nodes()));
        for (int g = 0; g < op.n_nodes(); ++g) {
            const double x = op.x(g);
            extra[z(g)] = rot[z(g)] - op.metric(g).J * omega * 2.0 * r(x, 0.0) * 0.5;  // (r^2)_x = 2 r r_x
        }
        const OperatorFields f = gfix::mms_fields(op, M, &extra);
        const auto psi = gfix::mms_solve(op, f);
        const StreamGradient G = op.gradient(psi);
        double e1 = 0.0, e2 = 0.0;
        for (int g = 0; g < op.n_nodes(); ++g) {
            const double x = op.x(g), y = op.y(g);
            e1 = std::max(e1, std::abs(psi[z(g)] - M.psi(x, y)));
            const double k = mdot / (rho(x, y) * b(x, y) * r(x, y));
            const double Wm = k * G.psi_y[z(g)], Wt = -k * G.psi_x[z(g)];
            e2 = std::max(e2, std::hypot(Wm - k * M.psi_y(x, y), Wt + k * M.psi_x(x, y)));
        }
        ep.push_back(e1);
        ew.push_back(e2);
        measure(std::format("general mms r,b,omega psi error N={}", N), e1);
        measure(std::format("general mms r,b,omega W error N={}", N), e2);
    }
    CHECK(measure("general mms r,b,omega psi order 64->128", gfix::order(ep[2], ep[3])) > 1.85);
    CHECK(measure("general mms r,b,omega W order 64->128", gfix::order(ew[2], ew[3])) > 1.5);
    // the rotation source is not zero (omega and dr/dx enter)
    const GeneralOperator op(gfix::dirichlet_box(gfix::block(gfix::family(5), 9, 9, r, b)));
    const auto rot = op.rotation_source(omega);
    double big = 0.0;
    for (double v : rot) big = std::max(big, std::abs(v));
    CHECK(big > 0.0);
}

TEST_CASE("general: the rotation sense", "[turbolenza][general]")
{
    CHECK(theta_angular_speed(100.0, RotationSense::AgainstTheta) == -100.0);
    CHECK(theta_angular_speed(100.0, RotationSense::AlongTheta) == 100.0);
}

TEST_CASE("general: refusals of the operator", "[turbolenza][general]")
{
    const GeneralMesh box = gfix::dirichlet_box(gfix::block(gfix::family(1), 5, 5));
    const GeneralOperator op(box);
    OperatorFields f;
    f.coefficient.assign(z(op.n_nodes()), 1.0);
    CHECK_THROWS_AS(op.assemble(f), GeneralOperatorError);  // Dirichlet segments without values
    f.dirichlet = [](int) { return 0.0; };
    f.coefficient.pop_back();
    CHECK_THROWS_AS(op.assemble(f), GeneralOperatorError);  // wrong size
    f.coefficient.assign(z(op.n_nodes()), 1.0);
    f.coefficient[3] = -1.0;
    CHECK_THROWS_AS(op.assemble(f), GeneralOperatorError);  // non-positive coefficient
    f.coefficient[3] = std::nan("");
    CHECK_THROWS_AS(op.assemble(f), GeneralOperatorError);
    CHECK_THROWS_AS(contravariant_tensor(1.0, 2.0, 0.5, 1.0, 1.0), GeneralOperatorError);  // J = 0
    // an outlet without the flow-angle tangents
    GeneralMesh m2 = box;
    m2.segments[1].kind = SegmentKind::Outlet;
    const GeneralOperator op2(m2);
    OperatorFields f2;
    f2.coefficient.assign(z(op2.n_nodes()), 1.0);
    f2.dirichlet = [](int) { return 0.0; };
    CHECK_THROWS_AS(op2.assemble(f2), GeneralOperatorError);
}

TEST_CASE("general: manufactured solution on Aungier's H-grid layout with blade walls and transitions",
          "[turbolenza][general][mms]")
{
    // The stations next to the leading and trailing edges are the transitions of Aungier's grid: the two
    // periodic images of the last inlet station see the suction and the pressure wall.  A smooth
    // manufactured field (Dirichlet on the walls): the order shows whether the transition treatment is
    // consistent.
    const double P = 0.8;
    const gfix::Manufactured M = gfix::quasi_periodic(P);
    for (bool rounded : {false, true}) {
        std::vector<gfix::MmsError> e;
        for (int N : {10, 20, 40, 80, 160}) {
            const GeneralOperator op(gfix::hgrid_with_blade(N, static_cast<int>(0.8 * N) + 1, P, 0.12, rounded));
            e.push_back(gfix::mms_error(op, M, gfix::mms_solve(op, gfix::mms_fields(op, M))));
        }
        const std::string name = rounded ? "H-grid with a round-edged blade" : "H-grid with a sharp-edged blade";
        for (std::size_t k = 0; k < e.size(); ++k) measure(std::format("general mms psi error {} N={}", name, 10 << k), e[k].psi);
        const double p = measure("general mms psi order 80->160 " + name, gfix::order(e[3].psi, e[4].psi));
        measure("general mms grad order 80->160 " + name, gfix::order(e[3].grad, e[4].grad));
        CHECK(p > 0.5);  // measured: see the MEASURE lines (first order at the transitions)
    }
}
