#include "naca/Naca65.hpp"

#include "math/Constrained.hpp"
#include "math/FairCurve.hpp"
#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>

namespace BladeGen::Naca {

using Math::cross2;
using Math::deg2rad;
using Math::kPi;
using Math::left_normal;

namespace {

// ── constants (TURBOlenza geometry.py names in brackets) ─────────────────────
constexpr double TABLE_LE_RADIUS = 0.00687;   // [TABLE_LE_RADIUS] R_LE / c of the 65-010 form
constexpr double TABLE_T_MAX = 0.10;          // [TABLE_T_MAX]
constexpr double ARC_CAMBER_COEFF = 0.1103;   // [ARC_CAMBER_COEFF] Aungier Eq. 4-6
constexpr double XI_MIN = 1e-3;               // [XI_MIN] straight camber end pieces
constexpr double FIRST_TABLE_STATION = 0.005; // [FIRST_TABLE_STATION]
constexpr double LE_TANGENT_X = 0.005;        // [LE_TANGENT_X] nose axis
constexpr double KOVACH_START = 0.60;         // [KOVACH_START_STATION]
constexpr double KOVACH_TE_RADIUS = 0.008;    // [KOVACH_TE_RADIUS]
constexpr double KOVACH_BLEND = 0.005;        // BladeGen: C2 blend half-width at station 0.60
constexpr double EMERY_CLOSE_X = 0.995;       // BladeGen: start of the C2 closing piece
constexpr double TE_NEAREST_RANGE = 0.005;    // [TE_NEAREST_RANGE]
const double LE_FAIR_MIN_ARC = deg2rad(20.0); // [LE_FAIR_MIN_ARC]
constexpr int LE_FAIR_SPANS = 24;             // [LE_FAIR_SPANS]
constexpr int LE_FAIR_GAUSS_POINTS = 16;      // [LE_FAIR_GAUSS_POINTS]
constexpr int LE_FAIR_ITERATIONS = 6;         // [LE_FAIR_ITERATIONS]
constexpr int LE_FAIR_FIXED_CONTROLS = 3;     // [LE_FAIR_FIXED_CONTROLS]
constexpr double LE_ANCHOR_OFFSET = 0.5;      // [LE_ANCHOR_OFFSET] R_LE
constexpr double LE_ANCHOR_STEP = 0.125;      // [LE_ANCHOR_STEP] R_LE
constexpr double LE_ANCHOR_RANGE = 0.1;       // [LE_ANCHOR_RANGE] fraction of the source length
constexpr int LE_ANCHOR_BISECTIONS = 48;      // [LE_ANCHOR_BISECTIONS]
constexpr double LE_FAIRNESS_LIMIT = 4.0;     // [LE_FAIRNESS_LIMIT]
constexpr double LE_CLOSURE_TOL = 1e-10;      // [LE_CLOSURE_TOL] (1e-12 there)

// Aungier 2003 Table 4-1, p. 63, NACA 65 series, percent of chord, as printed
constexpr std::array<double, 26> X_PCT = {0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0,
                                          40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0};
constexpr std::array<double, 26> CAMBER_PCT = {0.0,   0.250, 0.350, 0.535, 0.930, 1.580, 2.120, 2.585, 3.365,
                                               3.980, 4.475, 4.860, 5.150, 5.355, 5.475, 5.515, 5.475, 5.355,
                                               5.150, 4.860, 4.475, 3.980, 3.365, 2.585, 1.580, 0.0};
constexpr std::array<double, 26> THICK_PCT = {0.0,   1.544, 1.864, 2.338, 3.480, 4.354, 5.294, 6.080, 7.332,
                                              8.286, 9.006, 9.520, 9.848, 9.992, 9.926, 9.624, 9.060, 8.292,
                                              7.364, 6.312, 5.168, 3.974, 2.770, 1.620, 0.612, 0.0};

double smoothstep5(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t));
}
double smoothstep5_d(double t)
{
    if (t <= 0.0 || t >= 1.0) return 0.0;
    return 30.0 * t * t * (1.0 - t) * (1.0 - t);
}
double smoothstep5_dd(double t)
{
    if (t <= 0.0 || t >= 1.0) return 0.0;
    return 60.0 * t * (1.0 - t) * (1.0 - 2.0 * t);
}

/// Golden-section minimum of f on [a, b].
template <class F>
double golden_min(F&& f, double a, double b, double tol = 1e-13)
{
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double c = b - g * (b - a), d = a + g * (b - a), fc = f(c), fd = f(d);
    while (std::abs(b - a) > tol) {
        if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); }
        else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
    }
    return 0.5 * (a + b);
}

struct GaussRule
{
    std::vector<double> x, w;
    GaussRule() { Math::gauss_legendre(LE_FAIR_GAUSS_POINTS, x, w); }
};
const GaussRule& gauss16()
{
    static const GaussRule g;
    return g;
}

std::vector<double> le_fair_knots()
{
    std::vector<double> U = {0.0, 0.0, 0.0};
    for (int i = 0; i <= LE_FAIR_SPANS; ++i) U.push_back(static_cast<double>(i) / LE_FAIR_SPANS);
    U.insert(U.end(), {1.0, 1.0, 1.0});
    return U;
}

/// All cubic basis functions (value and first derivative) at eta.
void cubic_basis(const std::vector<double>& U, double eta, std::vector<double>& N, std::vector<double>& dN)
{
    const int n_ctrl = static_cast<int>(U.size()) - 4;
    N.assign(static_cast<std::size_t>(n_ctrl), 0.0);
    dN.assign(static_cast<std::size_t>(n_ctrl), 0.0);
    eta = std::clamp(eta, 0.0, 1.0);
    const int span = Math::find_span(n_ctrl - 1, 3, eta, U);
    const auto d = Math::ders_basis_funs(span, eta, 3, 1, U);
    for (int j = 0; j <= 3; ++j) {
        N[static_cast<std::size_t>(span - 3 + j)] = d[0][static_cast<std::size_t>(j)];
        dN[static_cast<std::size_t>(span - 3 + j)] = d[1][static_cast<std::size_t>(j)];
    }
}

}  // namespace

std::string to_string(Variant v) { return v == Variant::EmeryBase ? "emery_base" : "kovach_sandercock"; }

Variant variant_from_string(const std::string& s)
{
    if (s == "kovach_sandercock" || s == "kovach" || s == "Kovach") return Variant::KovachSandercock;
    if (s == "emery_base" || s == "emery" || s == "Emery") return Variant::EmeryBase;
    throw std::invalid_argument("unknown NACA 65 variant '" + s + "' (kovach_sandercock or emery_base)");
}

void Naca65Params::validate() const
{
    if (!(std::isfinite(cl0) && cl0 >= 0.0 && cl0 <= 2.5)) throw std::invalid_argument("NACA 65: cl0 must be in [0, 2.5]");
    if (!(std::isfinite(t_max_over_c) && t_max_over_c >= 0.02 && t_max_over_c <= 0.20))
        throw std::invalid_argument("NACA 65: t_max_over_c must be in [0.02, 0.20]");
}

std::string Naca65Params::designation() const
{
    char b[64];
    std::snprintf(b, sizeof(b), "NACA 65-(%.0fA10)%.0f", 10.0 * cl0, 100.0 * t_max_over_c);
    return b;
}

const Table41& table41()
{
    static const Table41 t = [] {
        Table41 r{};
        for (std::size_t i = 0; i < 26; ++i) {
            r.x[i] = X_PCT[i] / 100.0;
            r.camber[i] = CAMBER_PCT[i] / 100.0;
            r.thickness_printed[i] = THICK_PCT[i] / 100.0;
            r.thickness[i] = r.thickness_printed[i];
        }
        r.thickness[4] = 0.03148;  // x/c = 2.5 %: 3.480 is a misprint (twice the 65-010 half ordinate 1.574 %)
        return r;
    }();
    return t;
}

std::array<double, 4> a10(double x, double cl0)
{
    if (!(x > 0.0 && x < 1.0)) throw std::invalid_argument("a10: x must lie in (0, 1)");
    const double k = cl0 / (4.0 * kPi);
    const double y = -k * ((1.0 - x) * std::log(1.0 - x) + x * std::log(x));
    const double y1 = k * std::log((1.0 - x) / x);
    const double y2 = -k * (1.0 / (1.0 - x) + 1.0 / x);
    const double y3 = k * (1.0 / (x * x) - 1.0 / ((1.0 - x) * (1.0 - x)));
    return {y, y1, y2, y3};
}

std::array<double, 3> equivalent_arc_angles(double cl0)
{
    const double camber = 4.0 * std::atan(ARC_CAMBER_COEFF * cl0);
    return {camber, 0.5 * camber, 0.5 * camber};
}

// ── WhewellFair ──────────────────────────────────────────────────────────────

double WhewellFair::r(double phi) const
{
    if (phi <= min_arc || control.empty()) return radius;
    const double eta = std::clamp((phi - min_arc) / (anchor_angle - min_arc), 0.0, 1.0);
    const auto U = le_fair_knots();
    std::vector<double> N, dN;
    cubic_basis(U, eta, N, dN);
    double v = 0.0;
    for (std::size_t j = 0; j < N.size(); ++j) v += N[j] * control[j];
    return v;
}

double WhewellFair::r_slope(double phi) const
{
    if (phi <= min_arc || control.empty()) return 0.0;
    const double span = anchor_angle - min_arc;
    const double eta = std::clamp((phi - min_arc) / span, 0.0, 1.0);
    const auto U = le_fair_knots();
    std::vector<double> N, dN;
    cubic_basis(U, eta, N, dN);
    double v = 0.0;
    for (std::size_t j = 0; j < dN.size(); ++j) v += dN[j] * control[j];
    return v / span;
}

Vec2 WhewellFair::tangent(double phi) const
{
    const double a = apex_angle + sense * phi;
    return {std::cos(a), std::sin(a)};
}

Vec2 WhewellFair::point(double phi) const
{
    const double start = apex_angle;
    const double ang = start + sense * phi;
    // circle part over the whole [0, phi] with r = R, closed form
    Vec2 out = apex + radius * sense * Vec2(std::sin(ang) - std::sin(start), std::cos(start) - std::cos(ang));
    if (phi > min_arc && !control.empty()) {
        // + int_{min_arc}^{phi} (r - R) (cos, sin) dpsi, Gauss per span
        const auto& g = gauss16();
        const double span = anchor_angle - min_arc;
        const double eta_end = std::clamp((phi - min_arc) / span, 0.0, 1.0);
        for (int k = 0; k < LE_FAIR_SPANS; ++k) {
            const double lo = static_cast<double>(k) / LE_FAIR_SPANS;
            if (lo >= eta_end) break;
            const double hi = std::min(static_cast<double>(k + 1) / LE_FAIR_SPANS, eta_end);
            const double c = 0.5 * (lo + hi), h = 0.5 * (hi - lo);
            for (std::size_t q = 0; q < g.x.size(); ++q) {
                const double eta = c + h * g.x[q];
                const double ps = min_arc + span * eta;
                out += span * h * g.w[q] * (r(ps) - radius) * tangent(ps);
            }
        }
    }
    return out;
}

namespace {

/// The fairest monotone-curvature fair from the nose to the anchor (TURBOlenza _solve_le_fair).
std::optional<WhewellFair> solve_fair(const Vec2& apex, const Vec2& axis, double R, double side, const Vec2& point,
                                      const Vec2& tangent, double kappa, double anchor_s)
{
    const double sense = -side;
    const Vec2 normal = left_normal(axis);
    const Vec2 start_tangent = side * normal;
    const double anchor_angle = sense * std::atan2(cross2(start_tangent, tangent), start_tangent.dot(tangent));
    if (!(kappa * sense > 0.0 && anchor_angle > LE_FAIR_MIN_ARC && 1.0 / std::abs(kappa) > R)) return std::nullopt;
    const double end_radius = 1.0 / std::abs(kappa);
    const double span = anchor_angle - LE_FAIR_MIN_ARC;
    const auto U = le_fair_knots();
    const int n_control = LE_FAIR_SPANS + 3;
    const int first = LE_FAIR_FIXED_CONTROLS;
    const int nu = n_control - first;
    const auto& g = gauss16();
    const int nq = LE_FAIR_SPANS * LE_FAIR_GAUSS_POINTS;
    std::vector<double> eta(static_cast<std::size_t>(nq)), weight(static_cast<std::size_t>(nq));
    Eigen::MatrixXd cum(nq, nu), cum_slope(nq, nu);
    std::vector<double> N, dN;
    for (int k = 0, q = 0; k < LE_FAIR_SPANS; ++k) {
        const double lo = static_cast<double>(k) / LE_FAIR_SPANS, h = 0.5 / LE_FAIR_SPANS;
        for (int i = 0; i < LE_FAIR_GAUSS_POINTS; ++i, ++q) {
            const auto uq = static_cast<std::size_t>(q);
            eta[uq] = lo + h + h * g.x[static_cast<std::size_t>(i)];
            weight[uq] = h * g.w[static_cast<std::size_t>(i)] * span;  // d phi
            cubic_basis(U, eta[uq], N, dN);
            double sN = 0.0, sdN = 0.0;
            for (int j = n_control - 1; j >= first; --j) {  // S_j = sum_{k >= j} N_k
                sN += N[static_cast<std::size_t>(j)];
                sdN += dN[static_cast<std::size_t>(j)];
                cum(q, j - first) = sN;
                cum_slope(q, j - first) = sdN / span;
            }
        }
    }
    const double apex_angle = std::atan2(start_tangent.y(), start_tangent.x());
    Eigen::MatrixXd E(3, nu);
    E.setZero();
    for (int q = 0; q < nq; ++q) {
        const auto uq = static_cast<std::size_t>(q);
        const double ang = apex_angle + sense * (LE_FAIR_MIN_ARC + span * eta[uq]);
        E.row(0) += weight[uq] * std::cos(ang) * cum.row(q);
        E.row(1) += weight[uq] * std::sin(ang) * cum.row(q);
    }
    E.row(2).setOnes();
    const double end_angle = apex_angle + sense * anchor_angle;
    const Vec2 circle = R * sense * Vec2(std::sin(end_angle) - std::sin(apex_angle), std::cos(apex_angle) - std::cos(end_angle));
    Eigen::Vector3d f;
    f << point.x() - apex.x() - circle.x(), point.y() - apex.y() - circle.y(), end_radius - R;
    const double scale = end_radius - R;
    Eigen::MatrixXd Es = E * scale;
    Eigen::Vector3d fs = f;
    for (int i = 0; i < 3; ++i) {
        const double nrm = Es.row(i).norm();
        if (!(nrm > 0.0)) return std::nullopt;
        Es.row(i) /= nrm;
        fs[i] /= nrm;
    }
    Eigen::VectorXd r_hat(nq);
    for (int q = 0; q < nq; ++q) r_hat[q] = R + scale * eta[static_cast<std::size_t>(q)];
    Eigen::VectorXd c(nu);
    for (int it = 0; it < LE_FAIR_ITERATIONS; ++it) {
        Eigen::MatrixXd A(nq, nu);
        for (int q = 0; q < nq; ++q) {
            const double wq = std::sqrt(R * R * R * weight[static_cast<std::size_t>(q)] / std::pow(r_hat[q], 5));
            A.row(q) = wq * scale * cum_slope.row(q);
        }
        if (!A.allFinite()) return std::nullopt;
        const auto gamma = Math::lsei(A, Eigen::VectorXd::Zero(nq), Es, fs);
        if (!gamma) return std::nullopt;
        c = scale * gamma->cwiseMax(0.0);
        r_hat = (R + (cum * c).array()).matrix();
    }
    if (!c.allFinite()) return std::nullopt;
    WhewellFair fair;
    fair.apex = apex;
    fair.apex_angle = apex_angle;
    fair.sense = sense;
    fair.radius = R;
    fair.min_arc = LE_FAIR_MIN_ARC;
    fair.anchor_angle = anchor_angle;
    fair.anchor_s = anchor_s;
    fair.control.assign(static_cast<std::size_t>(first), R);
    double acc = R;
    for (int j = 0; j < nu; ++j) fair.control.push_back(acc += c[j]);
    double fairness = 0.0;
    for (int q = 0; q < nq; ++q) {
        const double rq = R + cum.row(q).dot(c), sl = cum_slope.row(q).dot(c);
        fairness += weight[static_cast<std::size_t>(q)] * sl * sl / std::pow(rq, 5);
    }
    fair.fairness = R * R * R * fairness;
    if (!((fair.point(anchor_angle) - point).norm() <= LE_CLOSURE_TOL)) return std::nullopt;
    return fair;
}

}  // namespace

// ── Naca65Section ────────────────────────────────────────────────────────────

Naca65Section::Naca65Section(Naca65Params p) : p_(p)
{
    p_.validate();
    build();
}

double Naca65Section::arc_station_raw(double x) const
{
    // arc length of the camber line from XI_MIN to x (x in [XI_MIN, 1 - XI_MIN])
    x = std::clamp(x, XI_MIN, 1.0 - XI_MIN);
    const auto it = std::upper_bound(arc_x_.begin(), arc_x_.end(), x);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - arc_x_.begin()), 1, arc_x_.size() - 1) - 1;
    auto speed = [&](double q) { const double y1 = a10(q, p_.cl0)[1]; return std::sqrt(1.0 + y1 * y1); };
    return arc_s_[k] + Math::gauss_legendre10(speed, arc_x_[k], x);
}

double Naca65Section::station(double x) const
{
    const double s = arc_station_raw(x) + end_piece_ * std::min(x / XI_MIN, 1.0) + end_piece_ * std::max((x - 1.0) / XI_MIN + 1.0, 0.0);
    return std::clamp(s / S_, 0.0, 1.0);
}

double Naca65Section::x_of_station(double sig) const
{
    if (sig <= 0.0) return 0.0;
    if (sig >= 1.0) return 1.0;
    return Math::brent_root([&](double q) { return station(q) - sig; }, 0.0, 1.0, 1e-15);
}

Vec2 Naca65Section::camber_point(double x) const
{
    if (x <= 0.0 || x >= 1.0) return {std::clamp(x, 0.0, 1.0), 0.0};
    return {x, a10(x, p_.cl0)[0]};
}

double Naca65Section::table_half_thickness(double sig) const
{
    sig = std::clamp(sig, 0.0, 1.0);
    return p_.t_max_over_c / TABLE_T_MAX * H_.eval(std::sqrt(sig))[0];
}

namespace {
/// h, dh/dsigma, d2h/dsigma2 of the table law (z = sqrt(sigma)).
std::array<double, 3> table_half3(const Math::BSpline<1>& H, double scale, double sig)
{
    sig = std::clamp(sig, 1e-12, 1.0);
    const double z = std::sqrt(sig);
    const auto d = H.derivs(z, 2);
    const double h1 = d[1][0] / (2.0 * z);
    const double h2 = d[2][0] / (4.0 * z * z) - d[1][0] / (4.0 * z * z * z);
    return {scale * d[0][0], scale * h1, scale * h2};
}
}  // namespace

double Naca65Section::half_thickness(double sig) const
{
    const double scale = p_.t_max_over_c / TABLE_T_MAX;
    const double ht = table_half3(H_, scale, sig)[0];
    if (p_.variant == Variant::EmeryBase || sig <= KOVACH_START - KOVACH_BLEND) return ht;
    const double line = h_start_ + (h_tan_ - h_start_) * (sig - KOVACH_START) / (sig_tan_ - KOVACH_START);
    if (sig >= KOVACH_START + KOVACH_BLEND) return line;
    const double b = smoothstep5((sig - (KOVACH_START - KOVACH_BLEND)) / (2.0 * KOVACH_BLEND));
    return (1.0 - b) * ht + b * line;
}

std::array<Vec2, 3> Naca65Section::surface(double x, double side) const
{
    if (!(x > 0.0 && x < 1.0)) throw std::invalid_argument("Naca65Section::surface: x must lie in (0, 1)");
    const auto c = a10(x, p_.cl0);
    const double y1 = c[1], y2 = c[2], y3 = c[3];
    const double q = 1.0 + y1 * y1, sq = std::sqrt(q);
    const double phi = std::atan(y1), ph1 = y2 / q, ph2 = (y3 * q - 2.0 * y1 * y2 * y2) / (q * q);
    const Vec2 n(-std::sin(phi), std::cos(phi)), tdir(std::cos(phi), std::sin(phi));
    const Vec2 n1 = -ph1 * tdir;
    const Vec2 n2 = -ph2 * tdir - ph1 * ph1 * n;
    const double sig = station(x);
    const double s1 = sq / S_, s2 = y1 * y2 / (sq * S_);
    // h(sigma) and its derivatives, Kovach taper and blend included
    const double scale = p_.t_max_over_c / TABLE_T_MAX;
    auto h3 = table_half3(H_, scale, sig);
    if (p_.variant == Variant::KovachSandercock && sig > KOVACH_START - KOVACH_BLEND) {
        const double sl = (h_tan_ - h_start_) / (sig_tan_ - KOVACH_START);
        const std::array<double, 3> L = {h_start_ + sl * (sig - KOVACH_START), sl, 0.0};
        if (sig >= KOVACH_START + KOVACH_BLEND) h3 = L;
        else {
            const double w = 2.0 * KOVACH_BLEND, t = (sig - (KOVACH_START - KOVACH_BLEND)) / w;
            const double b = smoothstep5(t), b1 = smoothstep5_d(t) / w, b2 = smoothstep5_dd(t) / (w * w);
            const std::array<double, 3> T = h3;
            h3 = {(1 - b) * T[0] + b * L[0], (1 - b) * T[1] + b * L[1] + b1 * (L[0] - T[0]),
                  (1 - b) * T[2] + b * L[2] + 2.0 * b1 * (L[1] - T[1]) + b2 * (L[0] - T[0])};
        }
    }
    const double h = h3[0], hx = h3[1] * s1, hxx = h3[2] * s1 * s1 + h3[1] * s2;
    const Vec2 P = Vec2(x, c[0]) + side * h * n;
    const Vec2 P1 = Vec2(1.0, y1) + side * (hx * n + h * n1);
    const Vec2 P2 = Vec2(0.0, y2) + side * (hxx * n + 2.0 * hx * n1 + h * n2);
    return {P, P1, P2};
}

void Naca65Section::build()
{
    const auto& T = table41();
    // thickness law: quintic in z = sqrt(station), h''' = h'''' = 0 at both ends
    {
        std::vector<double> z(26), U;
        std::vector<Math::Constraint<1>> cons;
        for (std::size_t i = 0; i < 26; ++i) z[i] = std::sqrt(T.x[i]);
        for (int i = 0; i < 6; ++i) U.push_back(z.front());
        for (std::size_t i = 1; i + 1 < 26; ++i) U.push_back(z[i]);
        for (int i = 0; i < 6; ++i) U.push_back(z.back());
        cons.push_back({z.front(), 0, Math::Vec<1>(0.5 * T.thickness.front())});
        cons.push_back({z.front(), 3, Math::Vec<1>(0.0)});
        cons.push_back({z.front(), 4, Math::Vec<1>(0.0)});
        for (std::size_t i = 1; i + 1 < 26; ++i) cons.push_back({z[i], 0, Math::Vec<1>(0.5 * T.thickness[i])});
        cons.push_back({z.back(), 4, Math::Vec<1>(0.0)});
        cons.push_back({z.back(), 3, Math::Vec<1>(0.0)});
        cons.push_back({z.back(), 0, Math::Vec<1>(0.5 * T.thickness.back())});
        H_ = Math::solve_bspline<1>(5, U, cons);
    }
    // camber arc length table (cosine stations on [XI_MIN, 1 - XI_MIN])
    {
        const int n = 2001;
        const double edge = std::acos(1.0 - 2.0 * XI_MIN);
        arc_x_.resize(n);
        arc_s_.assign(n, 0.0);
        auto speed = [&](double q) { const double y1 = a10(q, p_.cl0)[1]; return std::sqrt(1.0 + y1 * y1); };
        for (int i = 0; i < n; ++i)
            arc_x_[static_cast<std::size_t>(i)] = 0.5 * (1.0 - std::cos(edge + (kPi - 2.0 * edge) * i / (n - 1)));
        for (int i = 1; i < n; ++i) {
            const auto ui = static_cast<std::size_t>(i);
            arc_s_[ui] = arc_s_[ui - 1] + Math::gauss_legendre10(speed, arc_x_[ui - 1], arc_x_[ui]);
        }
        end_piece_ = std::hypot(XI_MIN, a10(XI_MIN, p_.cl0)[0]);
        S_ = arc_s_.back() + 2.0 * end_piece_;
    }
    meta_.designation = p_.designation();
    meta_.camber_length = S_;
    const auto ang = equivalent_arc_angles(p_.cl0);
    meta_.camber_angle = ang[0];
    meta_.chi_le = ang[1];
    meta_.chi_te = ang[2];
    // Kovach taper (TURBOlenza _kovach_closure)
    if (p_.variant == Variant::KovachSandercock) {
        const double r = KOVACH_TE_RADIUS;
        const double xc = Math::brent_root([&](double q) { return std::hypot(1.0 - q, a10(q, p_.cl0)[0]) - r; },
                                           1.0 - 50.0 * r, 1.0 - XI_MIN, 1e-15);
        const double sc = station(xc);
        const double run = (sc - KOVACH_START) * S_;
        h_start_ = table_half_thickness(KOVACH_START);
        const double taper = std::acos(r / std::hypot(h_start_, run)) - std::atan2(run, h_start_);
        sig_tan_ = sc + r * std::sin(taper) / S_;
        h_tan_ = r * std::cos(taper);
        meta_.te_radius = r;
        meta_.te_center = Vec2(xc, a10(xc, p_.cl0)[0]);
    }
    // maximum thickness
    {
        const double sm = golden_min([&](double s) { return -half_thickness(s); }, 0.35, 0.45, 1e-12);
        meta_.max_thickness_station = sm;
        meta_.max_thickness = 2.0 * half_thickness(sm);
    }
    // LE nose (TURBOlenza _naca65_nose)
    meta_.le_radius = TABLE_LE_RADIUS * std::pow(p_.t_max_over_c / TABLE_T_MAX, 2);
    const double phi_le = std::atan(a10(LE_TANGENT_X, p_.cl0)[1]);
    meta_.le_center = meta_.le_radius * Vec2(std::cos(phi_le), std::sin(phi_le));
    build_le(0);
    build_le(1);
    build_te();
}

void Naca65Section::build_le(int k)
{
    const double side = k == 0 ? 1.0 : -1.0;
    const double R = meta_.le_radius;
    const Vec2 apex(0.0, 0.0), axis = meta_.le_center / R;
    // the source: the table surface from station 0.005, in arc length
    const double x_first = x_of_station(FIRST_TABLE_STATION);
    std::vector<double> xs, ss;
    {
        const int n = 3000;
        const double x_hi = 0.4;
        for (int i = 0; i < n; ++i) xs.push_back(x_first * std::pow(x_hi / x_first, static_cast<double>(i) / (n - 1)));
        ss.assign(xs.size(), 0.0);
        auto speed = [&](double q) { return surface(q, side)[1].norm(); };
        for (std::size_t i = 1; i < xs.size(); ++i) ss[i] = ss[i - 1] + Math::gauss_legendre10(speed, xs[i - 1], xs[i]);
    }
    // the whole source length (to the TE), for the anchor range
    const double total = ss.back() + Math::integrate([&](double q) { return surface(q, side)[1].norm(); }, xs.back(),
                                                     1.0 - 10.0 * XI_MIN, 1e-10);
    auto x_at = [&](double s) {
        const auto it = std::upper_bound(ss.begin(), ss.end(), s);
        const std::size_t i = std::clamp<std::size_t>(static_cast<std::size_t>(it - ss.begin()), 1, ss.size() - 1);
        double x = xs[i - 1] + (xs[i] - xs[i - 1]) * (s - ss[i - 1]) / (ss[i] - ss[i - 1]);
        for (int it2 = 0; it2 < 3; ++it2) {
            const double si = ss[i - 1] + Math::gauss_legendre10([&](double q) { return surface(q, side)[1].norm(); }, xs[i - 1], x);
            x -= (si - s) / surface(x, side)[1].norm();
        }
        return x;
    };
    auto trial = [&](double s) {
        const auto P = surface(x_at(s), side);
        const double sp = P[1].norm();
        return solve_fair(apex, axis, R, side, P[0], P[1] / sp, cross2(P[1], P[2]) / (sp * sp * sp), s);
    };
    const double step = LE_ANCHOR_STEP * R, limit = LE_ANCHOR_RANGE * total;
    double prev = -1.0, s = 0.0;
    for (int i = 0;; ++i) {
        s = std::min(i * step, limit);
        if (trial(s)) break;
        if (s >= limit) throw std::invalid_argument("NACA 65: no monotone LE fair lands on the surface (" + p_.designation() + ")");
        prev = s;
    }
    double earliest = s;
    if (prev >= 0.0) {
        double lo = prev, hi = s;
        for (int i = 0; i < LE_ANCHOR_BISECTIONS; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (!(lo < mid && mid < hi)) break;
            if (trial(mid)) hi = mid;
            else lo = mid;
        }
        earliest = hi;
    }
    const double anchor = std::min(earliest + LE_ANCHOR_OFFSET * R, limit);
    auto fair = trial(anchor);
    if (!fair) throw std::invalid_argument("NACA 65: no LE fair at the offset anchor (" + p_.designation() + ")");
    fair->earliest_s = earliest;
    le_[k] = *fair;
    x_body_start_[k] = x_at(anchor);
    meta_.le_fairness[k] = fair->fairness;
    meta_.le_anchor_x[k] = x_body_start_[k];
    meta_.le_anchor_angle_deg[k] = Math::rad2deg(fair->anchor_angle);
}

void Naca65Section::build_te()
{
    const Vec2 te(1.0, 0.0);
    if (p_.variant == Variant::EmeryBase) {
        // C2 closing piece from x = 0.995 to the sharp TE, along TURBOlenza's closing direction
        for (int k = 0; k < 2; ++k) {
            const double side = k == 0 ? 1.0 : -1.0;
            const double xb = EMERY_CLOSE_X, D = 1.0 - xb;
            const auto P = surface(xb, side);
            const Vec2 dir = Math::unit<2>(Vec2(te - surface(1.0 - XI_MIN, side)[0]));
            std::vector<double> U = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
            std::vector<Math::Constraint<2>> c = {{0.0, 0, P[0]}, {0.0, 1, D * P[1]}, {0.0, 2, D * D * P[2]},
                                                  {1.0, 2, Vec2::Zero()}, {1.0, 1, D * P[1].norm() * dir}, {1.0, 0, te}};
            te_close_[k] = Math::solve_bspline<2>(5, U, c);
            x_body_end_[k] = xb;
            meta_.te_start_x[k] = xb;
        }
        return;
    }
    // Kovach: each surface ends at its nearest point to the TE circle centre
    const Vec2 cen = meta_.te_center;
    const double xt = x_of_station(sig_tan_);
    for (int k = 0; k < 2; ++k) {
        const double side = k == 0 ? 1.0 : -1.0;
        x_body_end_[k] = golden_min([&](double q) { return (surface(q, side)[0] - cen).norm(); }, xt - TE_NEAREST_RANGE,
                                    std::min(xt + TE_NEAREST_RANGE, 1.0 - XI_MIN), 1e-12);
        meta_.te_start_x[k] = x_body_end_[k];
    }
    // G2 fair connector SS end -> TE point -> PS end, with the circle's tangent and curvature at the TE point
    const auto A = surface(x_body_end_[0], 1.0), B = surface(x_body_end_[1], -1.0);
    const double sa = A[1].norm(), sb = B[1].norm();
    const Vec2 t0 = A[1] / sa, t1 = -B[1] / sb;
    const double k0 = cross2(A[1], A[2]) / (sa * sa * sa), k1 = -cross2(B[1], B[2]) / (sb * sb * sb);
    const Vec2 axis = Math::unit<2>(Vec2(cen - te));
    const Vec2 tm = left_normal(axis);
    const double km = -1.0 / meta_.te_radius;  // the loop runs clockwise
    const double L0 = (te - A[0]).norm(), L1 = (B[0] - te).norm(), L = L0 + L1;
    te_umid_ = L0 / L;
    std::vector<double> U = {0, 0, 0, 0, 0, 0};
    for (int i = 1; i <= 6; ++i) U.push_back(i / 7.0);
    U.insert(U.end(), {1, 1, 1, 1, 1, 1});
    const Vec2 ex(1, 0), ey(0, 1);
    std::vector<Math::LinearCondition> c;
    auto pos = [&](double u, const Vec2& q) { c.push_back({u, 0, ex, q.x()}); c.push_back({u, 0, ey, q.y()}); };
    auto g2 = [&](double u, const Vec2& t, double kap, double speed) {
        c.push_back({u, 1, ex, speed * t.x()});
        c.push_back({u, 1, ey, speed * t.y()});
        c.push_back({u, 2, left_normal(t), kap * speed * speed});
    };
    pos(0.0, A[0]);
    pos(te_umid_, te);
    pos(1.0, B[0]);
    g2(0.0, t0, k0, L);
    g2(te_umid_, tm, km, L);
    g2(1.0, t1, k1, L);
    te_conn_ = Math::fair_bspline(5, U, c, 3);
    double dev = 0.0;
    for (int i = 0; i <= 400; ++i) dev = std::max(dev, std::abs((te_conn_.eval(i / 400.0) - cen).norm() - meta_.te_radius));
    meta_.te_circle_deviation = dev;
}

Naca65Section::Point3 Naca65Section::te_piece(int k, double u) const
{
    if (p_.variant == Variant::EmeryBase) {
        const auto d = te_close_[k].derivs(u, 2);
        return {d[0], d[1], d[2]};
    }
    // SS: connector u -> u_mid; PS: 1 -> u_mid (reversed: the side runs LE -> TE)
    const double w = k == 0 ? u * te_umid_ : 1.0 - u * (1.0 - te_umid_);
    const auto d = te_conn_.derivs(w, 2);
    const double f = k == 0 ? te_umid_ : -(1.0 - te_umid_);
    return {d[0], f * d[1], f * f * d[2]};
}

Profile::CurveSample Naca65Section::sample_side(Profile::Side side_id, int n_edge, int n_body) const
{
    if (n_edge < 2 || n_body < 2) throw std::invalid_argument("Naca65Section sampling: need >= 2 points per piece");
    const int k = static_cast<int>(side_id);
    const double side = k == 0 ? 1.0 : -1.0;
    const double conv = k == 0 ? -1.0 : 1.0;  // convexity sign of the LE -> TE signed curvature
    Profile::CurveSample c;
    double t = 0.0;
    auto push = [&](const Vec2& p, double ksigned, int piece) {
        c.xy.push_back(p);
        c.kappa.push_back(conv * ksigned);
        c.piece.push_back(piece);
        c.t.push_back(t);
        t += 1.0;
    };
    const auto& F = le_[k];
    for (int i = 0; i < n_edge; ++i) {
        const double phi = F.anchor_angle * i / (n_edge - 1);
        push(F.point(phi), F.sense / F.r(phi), static_cast<int>(Profile::Piece::LeadingEdge));
    }
    const int body = k == 0 ? static_cast<int>(Profile::Piece::SuctionSide) : static_cast<int>(Profile::Piece::PressureSide);
    const double za = std::sqrt(x_body_start_[k]), zb = std::sqrt(x_body_end_[k]);
    for (int i = 1; i < n_body; ++i) {
        const double u = 0.5 * (1.0 - std::cos(kPi * i / (n_body - 1)));
        const double z = za + (zb - za) * u;
        const auto P = surface(z * z, side);
        const double sp = P[1].norm();
        push(P[0], cross2(P[1], P[2]) / (sp * sp * sp), body);
    }
    for (int i = 1; i < n_edge; ++i) {
        const auto P = te_piece(k, static_cast<double>(i) / (n_edge - 1));
        const double sp = P.d1.norm();
        push(i == n_edge - 1 ? Vec2(1.0, 0.0) : P.p, sp > 0 ? cross2(P.d1, P.d2) / (sp * sp * sp) : 0.0,
             static_cast<int>(Profile::Piece::TrailingEdge));
    }
    c.s.assign(c.xy.size(), 0.0);
    for (std::size_t i = 1; i < c.xy.size(); ++i) c.s[i] = c.s[i - 1] + (c.xy[i] - c.xy[i - 1]).norm();
    return c;
}

Profile::CurveSample Naca65Section::sample_loop(int n_edge, int n_body) const
{
    const auto a = sample_side(Profile::Side::Suction, n_edge, n_body);
    const auto b = sample_side(Profile::Side::Pressure, n_edge, n_body);
    Profile::CurveSample c = a;
    for (std::size_t i = b.xy.size() - 1; i-- > 0;) {
        c.xy.push_back(b.xy[i]);
        c.kappa.push_back(b.kappa[i]);
        c.piece.push_back(b.piece[i]);
        c.t.push_back(c.t.back() + 1.0);
    }
    c.s.assign(c.xy.size(), 0.0);
    for (std::size_t i = 1; i < c.xy.size(); ++i) c.s[i] = c.s[i - 1] + (c.xy[i] - c.xy[i - 1]).norm();
    return c;
}

Profile::QualityReport check_quality(const Naca65Section& s, const Profile::QualityOptions& opt)
{
    Profile::QualityReport rep;
    rep.ss = s.sample_side(Profile::Side::Suction, opt.n_edge, opt.n_body);
    rep.ps = s.sample_side(Profile::Side::Pressure, opt.n_edge, opt.n_body);
    const auto loop = s.sample_loop(opt.n_edge / 2 + 2, opt.n_body / 2 + 2);
    Profile::check_samples(rep, loop.xy, opt, !s.sharp_te());
    const auto& m = s.metadata();
    for (int k = 0; k < 2; ++k)
        if (m.le_fairness[k] > LE_FAIRNESS_LIMIT)
            rep.issues.push_back({Profile::Severity::Error, "le_fair_steep",
                                  std::string("LE fair too steep on the ") + (k == 0 ? "suction" : "pressure") +
                                      " side: R_LE^3 int kappa_s^2 ds = " + std::to_string(m.le_fairness[k]) + " > 4",
                                  Vec2(0.0, 0.0), k});
    if (s.sharp_te())
        rep.issues.push_back({Profile::Severity::Info, "sharp_te", "sharp trailing edge (emery_base): a design corner",
                              Vec2(1.0, 0.0)});
    rep.min_thickness = 2.0 * s.half_thickness(0.97);
    std::stable_sort(rep.issues.begin(), rep.issues.end(),
                     [](const Profile::Issue& a, const Profile::Issue& b) { return a.severity > b.severity; });
    return rep;
}

}  // namespace BladeGen::Naca
