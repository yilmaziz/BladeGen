// The axisymmetric stream sheet r(m), b(m) (see StreamSheet.hpp).  Not in v1.
#include "turbolenza/StreamSheet.hpp"

#include "turbolenza/Detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>

namespace BladeGen::Turbolenza {

using detail::py_repr;

namespace {

// 8-point Gauss-Legendre nodes and weights on [-1, 1] (Abramowitz & Stegun, Table 25.4).
constexpr std::array<double, 4> kX = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
constexpr std::array<double, 4> kW = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

double gauss(double a, double b, const std::function<double(double)>& f)
{
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    double s = 0.0;
    for (std::size_t k = 0; k < kX.size(); ++k) s += kW[k] * (f(c - h * kX[k]) + f(c + h * kX[k]));
    return h * s;
}

}  // namespace

StreamSheet::StreamSheet(std::vector<double> m, std::vector<double> r, std::vector<double> b, double r_ref)
    : m_(std::move(m)), r_(std::move(r)), b_(std::move(b)), r_ref_(r_ref)
{
    const std::size_t n = m_.size();
    if (n < 2 || r_.size() != n || b_.size() != n)
        throw std::invalid_argument(std::format("StreamSheet: {} stations with {} radii and {} thicknesses: at least two "
                                                "stations, one r and one b each",
                                                n, r_.size(), b_.size()));
    for (std::size_t k = 0; k < n; ++k) {
        if (!std::isfinite(m_[k]) || (k > 0 && !(m_[k] > m_[k - 1])))
            throw std::invalid_argument(std::format("StreamSheet: m must be finite and strictly increasing (station {})", k));
        if (!(std::isfinite(r_[k]) && r_[k] > 0.0))
            throw std::invalid_argument(std::format("StreamSheet: r = {} at station {}: a radius must be positive",
                                                    py_repr(r_[k]), k));
        if (!(std::isfinite(b_[k]) && b_[k] > 0.0))
            throw std::invalid_argument(std::format("StreamSheet: b = {} at station {}: the sheet thickness must be "
                                                    "positive",
                                                    py_repr(b_[k]), k));
    }
    if (!(std::isfinite(r_ref_) && r_ref_ > 0.0))
        throw std::invalid_argument(std::format("StreamSheet: r_ref = {} must be positive", py_repr(r_ref_)));
    cylinder_ = std::all_of(r_.begin(), r_.end(), [&](double v) { return v == r_.front(); }) &&
                std::all_of(b_.begin(), b_.end(), [&](double v) { return v == b_.front(); });
    if (!cylinder_) {
        std::vector<Point> rb(n);
        for (std::size_t k = 0; k < n; ++k) rb[k] = {r_[k], b_[k]};
        spline_.emplace(m_, rb);
    }
}

void StreamSheet::check(double m) const
{
    if (!(m >= m_.front() && m <= m_.back()))
        throw std::out_of_range(std::format("the stream sheet is given for m in [{:.10g}, {:.10g}] chords; m = {:.10g} lies "
                                            "outside it (no extrapolation: shorten the extensions)",
                                            m_.front(), m_.back(), m));
}

double StreamSheet::radius(double m) const
{
    check(m);
    return cylinder_ ? r_.front() : (*spline_)(m, 0).m;
}

double StreamSheet::sin_phi(double m) const
{
    check(m);
    return cylinder_ ? 0.0 : (*spline_)(m, 1).m;
}

double StreamSheet::thickness(double m) const
{
    check(m);
    return cylinder_ ? b_.front() : (*spline_)(m, 0).y;
}

double StreamSheet::integrate(double a, double b, const std::function<double(double)>& f) const
{
    if (a == b) return 0.0;
    if (a > b) return -integrate(b, a, f);
    check(a);
    check(b);
    // the breakpoints: a, the table stations strictly inside (a, b), b; each piece at most kSideLinePiece
    std::vector<double> cuts{a};
    for (double s : m_)
        if (s > a && s < b) cuts.push_back(s);
    cuts.push_back(b);
    double total = 0.0;
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
        const double lo = cuts[k], hi = cuts[k + 1];
        const int pieces = std::max(1, static_cast<int>(std::ceil((hi - lo) / kSideLinePiece)));
        for (int q = 0; q < pieces; ++q) {
            const double x0 = lo + (hi - lo) * static_cast<double>(q) / static_cast<double>(pieces);
            const double x1 = q + 1 == pieces ? hi : lo + (hi - lo) * static_cast<double>(q + 1) / static_cast<double>(pieces);
            total += gauss(x0, x1, f);
        }
    }
    return total;
}

double StreamSheet::side_line_rise(double m_j, double m, double beta_j, double beta_far, double m_far) const
{
    const double span = m_far - m_j;
    if (!(std::isfinite(span) && span != 0.0))
        throw std::invalid_argument("StreamSheet::side_line_rise: the far station must differ from the junction");
    const double turn = beta_far - beta_j;
    return integrate(m_j, m, [&](double s) { return std::tan(beta_j + (s - m_j) / span * turn) / scale(s); });
}

}  // namespace BladeGen::Turbolenza
