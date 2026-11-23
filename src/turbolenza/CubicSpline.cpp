// Port of scipy.interpolate.CubicSpline, not-a-knot (see CubicSpline.hpp).
#include "turbolenza/CubicSpline.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Turbolenza::detail {

namespace {

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

/// LAPACK dgtsv for two right-hand sides (the m and y columns): dl, d, du are the sub-, main and
/// super-diagonal, overwritten; b holds the solution on return.
void gtsv(std::vector<double> dl, std::vector<double> d, std::vector<double> du, std::vector<std::array<double, 2>>& b)
{
    const int n = static_cast<int>(d.size());
    auto swap_step = [&](int i, bool fill) {
        const double fact = d[sz(i)] / dl[sz(i)];
        d[sz(i)] = dl[sz(i)];
        const double temp = d[sz(i + 1)];
        d[sz(i + 1)] = du[sz(i)] - fact * temp;
        if (fill) {
            dl[sz(i)] = du[sz(i + 1)];
            du[sz(i + 1)] = -fact * dl[sz(i)];
        }
        du[sz(i)] = temp;
        for (int j = 0; j < 2; ++j) {
            const double t = b[sz(i)][sz(j)];
            b[sz(i)][sz(j)] = b[sz(i + 1)][sz(j)];
            b[sz(i + 1)][sz(j)] = t - fact * b[sz(i + 1)][sz(j)];
        }
    };
    auto plain_step = [&](int i) {
        if (d[sz(i)] == 0.0) throw std::invalid_argument("CubicSpline: singular spline system");
        const double fact = dl[sz(i)] / d[sz(i)];
        d[sz(i + 1)] -= fact * du[sz(i)];
        for (int j = 0; j < 2; ++j) b[sz(i + 1)][sz(j)] -= fact * b[sz(i)][sz(j)];
    };
    for (int i = 0; i < n - 2; ++i) {
        if (std::abs(d[sz(i)]) >= std::abs(dl[sz(i)])) {
            plain_step(i);
            dl[sz(i)] = 0.0;
        } else {
            swap_step(i, true);
        }
    }
    if (n > 1) {
        const int i = n - 2;
        if (std::abs(d[sz(i)]) >= std::abs(dl[sz(i)]))
            plain_step(i);
        else
            swap_step(i, false);
    }
    if (d[sz(n - 1)] == 0.0) throw std::invalid_argument("CubicSpline: singular spline system");
    for (int j = 0; j < 2; ++j) {
        b[sz(n - 1)][sz(j)] /= d[sz(n - 1)];
        if (n > 1) b[sz(n - 2)][sz(j)] = (b[sz(n - 2)][sz(j)] - du[sz(n - 2)] * b[sz(n - 1)][sz(j)]) / d[sz(n - 2)];
        for (int i = n - 3; i >= 0; --i)
            b[sz(i)][sz(j)] =
                (b[sz(i)][sz(j)] - du[sz(i)] * b[sz(i + 1)][sz(j)] - dl[sz(i)] * b[sz(i + 2)][sz(j)]) / d[sz(i)];
    }
}

/// Dense 3 x 3 solve with partial pivoting (LAPACK gesv) for two right-hand sides.
void solve3(std::array<std::array<double, 3>, 3> A, std::vector<std::array<double, 2>>& b)
{
    for (int k = 0; k < 3; ++k) {
        int p = k;
        for (int r = k + 1; r < 3; ++r)
            if (std::abs(A[sz(r)][sz(k)]) > std::abs(A[sz(p)][sz(k)])) p = r;
        std::swap(A[sz(k)], A[sz(p)]);
        std::swap(b[sz(k)], b[sz(p)]);
        if (A[sz(k)][sz(k)] == 0.0) throw std::invalid_argument("CubicSpline: singular spline system");
        for (int r = k + 1; r < 3; ++r) {
            const double f = A[sz(r)][sz(k)] / A[sz(k)][sz(k)];
            for (int c = k; c < 3; ++c) A[sz(r)][sz(c)] -= f * A[sz(k)][sz(c)];
            for (int j = 0; j < 2; ++j) b[sz(r)][sz(j)] -= f * b[sz(k)][sz(j)];
        }
    }
    for (int j = 0; j < 2; ++j)
        for (int k = 2; k >= 0; --k) {
            double s = b[sz(k)][sz(j)];
            for (int c = k + 1; c < 3; ++c) s -= A[sz(k)][sz(c)] * b[sz(c)][sz(j)];
            b[sz(k)][sz(j)] = s / A[sz(k)][sz(k)];
        }
}

}  // namespace

CubicSpline::CubicSpline(std::vector<double> x, const std::vector<Point>& points) : x_(std::move(x))
{
    const int n = static_cast<int>(x_.size());
    if (n < 2 || points.size() != x_.size())
        throw std::invalid_argument("CubicSpline: needs at least two points and one x per point");
    for (int i = 0; i + 1 < n; ++i)
        if (!(x_[sz(i + 1)] > x_[sz(i)])) throw std::invalid_argument("CubicSpline: x must be strictly increasing");
    std::vector<std::array<double, 2>> y(sz(n));
    for (int i = 0; i < n; ++i) y[sz(i)] = {points[sz(i)].m, points[sz(i)].y};
    std::vector<double> dx(sz(n - 1));
    std::vector<std::array<double, 2>> slope(sz(n - 1));
    for (int i = 0; i + 1 < n; ++i) {
        dx[sz(i)] = x_[sz(i + 1)] - x_[sz(i)];
        for (int j = 0; j < 2; ++j) slope[sz(i)][sz(j)] = (y[sz(i + 1)][sz(j)] - y[sz(i)][sz(j)]) / dx[sz(i)];
    }
    std::vector<std::array<double, 2>> s(sz(n));
    if (n == 2) {
        s[0] = s[1] = slope[0];
    } else if (n == 3) {
        std::array<std::array<double, 3>, 3> A{};
        A[0] = {1.0, 1.0, 0.0};
        A[1] = {dx[1], 2.0 * (dx[0] + dx[1]), dx[0]};
        A[2] = {0.0, 1.0, 1.0};
        for (int j = 0; j < 2; ++j) {
            s[0][sz(j)] = 2.0 * slope[0][sz(j)];
            s[1][sz(j)] = 3.0 * (dx[0] * slope[1][sz(j)] + dx[1] * slope[0][sz(j)]);
            s[2][sz(j)] = 2.0 * slope[1][sz(j)];
        }
        solve3(A, s);
    } else {
        // banded rows of scipy: du[i] = A(i, i + 1), d[i] = A(i, i), dl[i] = A(i + 1, i)
        std::vector<double> dl(sz(n - 1), 0.0), d(sz(n), 0.0), du(sz(n - 1), 0.0);
        for (int i = 1; i + 1 < n; ++i) {
            d[sz(i)] = 2.0 * (dx[sz(i - 1)] + dx[sz(i)]);
            du[sz(i)] = dx[sz(i - 1)];
            dl[sz(i - 1)] = dx[sz(i)];
            for (int j = 0; j < 2; ++j)
                s[sz(i)][sz(j)] = 3.0 * (dx[sz(i)] * slope[sz(i - 1)][sz(j)] + dx[sz(i - 1)] * slope[sz(i)][sz(j)]);
        }
        {
            d[0] = dx[1];
            du[0] = x_[2] - x_[0];
            const double dd = x_[2] - x_[0];
            for (int j = 0; j < 2; ++j)
                s[0][sz(j)] =
                    ((dx[0] + 2.0 * dd) * dx[1] * slope[0][sz(j)] + dx[0] * dx[0] * slope[1][sz(j)]) / dd;
        }
        {
            d[sz(n - 1)] = dx[sz(n - 3)];
            dl[sz(n - 2)] = x_[sz(n - 1)] - x_[sz(n - 3)];
            const double dd = x_[sz(n - 1)] - x_[sz(n - 3)];
            const double a = dx[sz(n - 2)], c = dx[sz(n - 3)];
            for (int j = 0; j < 2; ++j)
                s[sz(n - 1)][sz(j)] = (a * a * slope[sz(n - 3)][sz(j)] + (2.0 * dd + a) * c * slope[sz(n - 2)][sz(j)]) / dd;
        }
        gtsv(dl, d, du, s);
    }
    c_.resize(sz(n - 1));
    for (int i = 0; i + 1 < n; ++i) {
        double cc[4][2];
        for (int j = 0; j < 2; ++j) {
            const double t = (s[sz(i)][sz(j)] + s[sz(i + 1)][sz(j)] - 2.0 * slope[sz(i)][sz(j)]) / dx[sz(i)];
            cc[0][j] = t / dx[sz(i)];
            cc[1][j] = (slope[sz(i)][sz(j)] - s[sz(i)][sz(j)]) / dx[sz(i)] - t;
            cc[2][j] = s[sz(i)][sz(j)];
            cc[3][j] = y[sz(i)][sz(j)];
        }
        for (int k = 0; k < 4; ++k) c_[sz(i)][sz(k)] = {cc[k][0], cc[k][1]};
    }
}

Point CubicSpline::operator()(double x, int nu) const
{
    if (nu < 0 || nu > 3) throw std::invalid_argument("CubicSpline: derivative order must lie in [0, 3]");
    if (std::isnan(x)) throw std::invalid_argument("CubicSpline: x is NaN");
    // scipy find_interval: x_i <= x < x_{i+1}, the last interval at x = x_{n-1}, extrapolation outside
    const int n = static_cast<int>(x_.size());
    int i = 0;
    if (x >= x_[sz(n - 1)])
        i = n - 2;
    else if (x > x_[0])
        i = static_cast<int>(std::upper_bound(x_.begin(), x_.end(), x) - x_.begin()) - 1;
    const double s = x - x_[sz(i)];
    // scipy evaluate_poly1: sum over ascending powers kp of c[3 - kp] z prefactor
    Point res{0.0, 0.0};
    double z = 1.0;
    for (int kp = 0; kp < 4; ++kp) {
        if (kp < nu) continue;
        double prefactor = 1.0;
        for (int k = kp; k > kp - nu; --k) prefactor *= static_cast<double>(k);
        const Point& c = c_[sz(i)][sz(3 - kp)];
        res.m = res.m + c.m * z * prefactor;
        res.y = res.y + c.y * z * prefactor;
        if (kp < 3) z *= s;
    }
    return res;
}

}  // namespace BladeGen::Turbolenza::detail
