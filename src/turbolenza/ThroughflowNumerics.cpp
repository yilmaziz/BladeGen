/// @file turbolenza/ThroughflowNumerics.cpp
/// @brief Aungier (2003) Eqs. 7-34 to 7-36 and 7-67 to 7-73 (see the header).
#include "turbolenza/ThroughflowNumerics.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace BladeGen::Turbolenza::Throughflow {

namespace {

std::size_t us(int k) { return static_cast<std::size_t>(k); }

void require_points(const std::vector<double>& x, const std::vector<double>& f, const char* who)
{
    if (x.size() < 3) throw std::invalid_argument(std::string(who) + ": at least three points are needed");
    if (f.size() != x.size()) throw std::invalid_argument(std::string(who) + ": x and f differ in length");
    for (std::size_t k = 0; k < x.size(); ++k) {
        if (!std::isfinite(x[k]) || !std::isfinite(f[k]))
            throw std::invalid_argument(std::string(who) + ": non-finite value at point " + std::to_string(k));
        if (k > 0 && !(x[k] > x[k - 1]))
            throw std::invalid_argument(std::string(who) + ": the abscissae must increase strictly (point " +
                                        std::to_string(k) + ")");
    }
}

/// Eq. 7-67 at interior point c.
double d1_interior(const std::vector<double>& x, const std::vector<double>& f, std::size_t c)
{
    const double m1 = x[c - 1], m2 = x[c], m3 = x[c + 1];
    const double f1 = f[c - 1], f2 = f[c], f3 = f[c + 1];
    return ((m3 - m2) * (f2 - f1) / (m2 - m1) + (m2 - m1) * (f3 - f2) / (m3 - m2)) / (m3 - m1);
}

/// Eq. 7-71 at interior point c.
double d2_interior(const std::vector<double>& x, const std::vector<double>& f, std::size_t c)
{
    const double y1 = x[c - 1], y2 = x[c], y3 = x[c + 1];
    const double f1 = f[c - 1], f2 = f[c], f3 = f[c + 1];
    return 2.0 / (y3 - y1) * ((f3 - f2) / (y3 - y2) - (f2 - f1) / (y2 - y1));
}

/// The Taylor centre of segment k (p. 172): its right end when interior (Eq. 7-72), else its left
/// end (Eq. 7-73, the last segment).
std::size_t centre_of(std::size_t n, int k)
{
    const std::size_t right = us(k) + 1;
    return right < n - 1 ? right : us(k);
}

}  // namespace

std::vector<double> derivative(const std::vector<double>& x, const std::vector<double>& f, EndRule ends)
{
    require_points(x, f, "derivative");
    const std::size_t n = x.size();
    std::vector<double> d(n);
    for (std::size_t c = 1; c + 1 < n; ++c) d[c] = d1_interior(x, f, c);
    if (ends == EndRule::TwoPoint) {
        d[0] = (f[1] - f[0]) / (x[1] - x[0]);                      // Eq. 7-70
        d[n - 1] = (f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2]);  // Eq. 7-70, mirrored
        return d;
    }
    {  // Eq. 7-68 at the first point
        const double m1 = x[0], m2 = x[1], m3 = x[2];
        const double f1 = f[0], f2 = f[1], f3 = f[2];
        d[0] = ((m3 - m1) * (f2 - f1) / (m2 - m1) - (m2 - m1) * (f3 - f1) / (m3 - m1)) / (m3 - m2);
    }
    {  // Eq. 7-69 at the last point
        const double m1 = x[n - 3], m2 = x[n - 2], m3 = x[n - 1];
        const double f1 = f[n - 3], f2 = f[n - 2], f3 = f[n - 1];
        d[n - 1] = ((m3 - m1) * (f3 - f2) / (m3 - m2) - (m3 - m2) * (f3 - f1) / (m3 - m1)) / (m2 - m1);
    }
    return d;
}

std::vector<double> second_derivative(const std::vector<double>& x, const std::vector<double>& f)
{
    require_points(x, f, "second_derivative");
    const std::size_t n = x.size();
    std::vector<double> d(n, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t c = 1; c + 1 < n; ++c) d[c] = d2_interior(x, f, c);
    return d;
}

double partial_segment_integral(const std::vector<double>& x, const std::vector<double>& f, int k, double xi)
{
    require_points(x, f, "partial_segment_integral");
    const std::size_t n = x.size();
    if (k < 0 || us(k) + 1 >= n) throw std::invalid_argument("partial_segment_integral: no segment " + std::to_string(k));
    const double a = x[us(k)], b = x[us(k) + 1];
    if (!(xi >= a && xi <= b)) throw std::invalid_argument("partial_segment_integral: xi outside the segment");
    const std::size_t c = centre_of(n, k);
    const double fc = f[c], d1 = d1_interior(x, f, c), d2 = d2_interior(x, f, c);
    const double ua = a - x[c], ub = xi - x[c];
    // int (fc + d1 u + d2 u^2 / 2) du from ua to ub: the series of Eqs. 7-72 / 7-73
    return fc * (ub - ua) + 0.5 * d1 * (ub * ub - ua * ua) + d2 * (ub * ub * ub - ua * ua * ua) / 6.0;
}

double segment_integral(const std::vector<double>& x, const std::vector<double>& f, int k)
{
    require_points(x, f, "segment_integral");
    const std::size_t n = x.size();
    if (k < 0 || us(k) + 1 >= n) throw std::invalid_argument("segment_integral: no segment " + std::to_string(k));
    const std::size_t c = centre_of(n, k);
    const double fc = f[c], d1 = d1_interior(x, f, c), d2 = d2_interior(x, f, c);
    if (c == us(k) + 1) {  // Eq. 7-72: point 2 = the right end (interior)
        const double h = x[c] - x[c - 1];
        return (fc - d1 * h / 2.0 + d2 * h * h / 6.0) * h;
    }
    const double h = x[c + 1] - x[c];  // Eq. 7-73: point 2 = the left end, point 3 the last point
    return (fc + d1 * h / 2.0 + d2 * h * h / 6.0) * h;
}

std::vector<double> cumulative_integral(const std::vector<double>& x, const std::vector<double>& f)
{
    require_points(x, f, "cumulative_integral");
    std::vector<double> I(x.size(), 0.0);
    for (std::size_t k = 0; k + 1 < x.size(); ++k) I[k + 1] = I[k] + segment_integral(x, f, static_cast<int>(k));
    return I;
}

double integral(const std::vector<double>& x, const std::vector<double>& f) { return cumulative_integral(x, f).back(); }

std::vector<double> LinearOdeSolution::evaluate(double w0) const
{
    std::vector<double> w(F.size());
    for (std::size_t k = 0; k < F.size(); ++k) w[k] = w0 * F[k] + G[k];
    return w;
}

LinearOdeSolution solve_linear_ode(const std::vector<double>& y, const std::vector<double>& f1,
                                   const std::vector<double>& f4)
{
    require_points(y, f1, "solve_linear_ode");
    require_points(y, f4, "solve_linear_ode");
    const std::size_t n = y.size();
    LinearOdeSolution s;
    const auto i1 = cumulative_integral(y, f1);
    s.F.resize(n);
    for (std::size_t k = 0; k < n; ++k) s.F[k] = std::exp(i1[k]);  // Eq. 7-36
    std::vector<double> q(n);
    for (std::size_t k = 0; k < n; ++k) q[k] = f4[k] / s.F[k];
    const auto iq = cumulative_integral(y, q);
    s.G.resize(n);
    for (std::size_t k = 0; k < n; ++k) s.G[k] = s.F[k] * iq[k];  // Eq. 7-35, second term
    return s;
}

}  // namespace BladeGen::Turbolenza::Throughflow
