#include "blade/Meridional.hpp"

#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Blade {

namespace {
constexpr int kPerSpan = 48;  // table points per knot span
}

MeridionalCurve::MeridionalCurve(Math::BSpline2 curve) : c_(std::move(curve)) { tabulate(); }

MeridionalCurve MeridionalCurve::through_points(const std::vector<Vec2>& xr)
{
    if (xr.size() < 2) throw std::invalid_argument("MeridionalCurve: need at least 2 points");
    return MeridionalCurve(Math::interpolate_cubic<2>(Math::chord_length_params<2>(xr), xr));
}

MeridionalCurve MeridionalCurve::line(const Vec2& a, const Vec2& b)
{
    if (!((b - a).norm() > 0.0)) throw std::invalid_argument("MeridionalCurve::line: coincident points");
    return MeridionalCurve(Math::BSpline2::clamped(3, {a, a + (b - a) / 3.0, a + 2.0 * (b - a) / 3.0, b}));
}

double MeridionalCurve::integrate_m(double ta, double tb) const
{
    return Math::integrate([&](double t) { return c_.derivs(t, 1)[1].norm(); }, ta, tb, 1e-14);
}

double MeridionalCurve::integrate_mp(double ta, double tb) const
{
    return Math::integrate([&](double t) {
        const auto d = c_.derivs(t, 1);
        return d[1].norm() / d[0].y();
    }, ta, tb, 1e-14);
}

void MeridionalCurve::tabulate()
{
    const auto& U = c_.knots();
    std::vector<double> brk;
    for (double k : U)
        if (k >= c_.t0() && k <= c_.t1() && (brk.empty() || k > brk.back())) brk.push_back(k);
    t_.clear();
    t_.push_back(brk.front());
    for (std::size_t s = 0; s + 1 < brk.size(); ++s)
        for (int i = 1; i <= kPerSpan; ++i) t_.push_back(brk[s] + (brk[s + 1] - brk[s]) * i / kPerSpan);
    m_.assign(t_.size(), 0.0);
    mp_.assign(t_.size(), 0.0);
    for (std::size_t i = 0; i < t_.size(); ++i) {
        const auto d = c_.derivs(t_[i], 1);
        if (!(d[0].y() > 0.0)) throw std::invalid_argument("MeridionalCurve: radius must stay > 0");
        if (!(d[1].norm() > 0.0)) throw std::invalid_argument("MeridionalCurve: zero parametric speed");
        if (i > 0) {
            m_[i] = m_[i - 1] + integrate_m(t_[i - 1], t_[i]);
            mp_[i] = mp_[i - 1] + integrate_mp(t_[i - 1], t_[i]);
        }
    }
}

double MeridionalCurve::t_of_m(double m) const
{
    m = std::clamp(m, 0.0, length());
    const auto it = std::upper_bound(m_.begin(), m_.end(), m);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - m_.begin()), 1, m_.size() - 1);
    const double ta = t_[k - 1], tb = t_[k];
    const double base = m_[k - 1];
    auto f = [&](double t) { return base + integrate_m(ta, t) - m; };
    if (f(ta) >= 0.0) return ta;
    if (f(tb) <= 0.0) return tb;
    return Math::brent_root(f, ta, tb, 1e-15);
}

Vec2 MeridionalCurve::at_m(double m) const
{
    if (m < 0.0) return c_.eval(c_.t0()) + m * tangent_at_m(0.0);
    if (m > length()) return c_.eval(c_.t1()) + (m - length()) * tangent_at_m(length());
    return c_.eval(t_of_m(m));
}

Vec2 MeridionalCurve::tangent_at_m(double m) const
{
    return Math::unit<2>(c_.derivs(t_of_m(m), 1)[1]);
}

double MeridionalCurve::mprime_of_m(double m) const
{
    auto ext = [](double r_end, double dr, double dm) {  // \int_0^dm ds / (r_end + dr s)
        if (std::abs(dr * dm) < 1e-12 * r_end) return dm / r_end;
        const double r = r_end + dr * dm;
        if (!(r > 0.0)) throw std::invalid_argument("MeridionalCurve: extrapolation reaches r <= 0");
        return std::log(r / r_end) / dr;
    };
    if (m < 0.0) {
        const Vec2 p = c_.eval(c_.t0()), t = tangent_at_m(0.0);
        return -ext(p.y(), -t.y(), -m);
    }
    if (m > length()) {
        const Vec2 p = c_.eval(c_.t1()), t = tangent_at_m(length());
        return mprime_length() + ext(p.y(), t.y(), m - length());
    }
    const double t = t_of_m(m);
    const auto it = std::upper_bound(t_.begin(), t_.end(), t);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - t_.begin()), 1, t_.size() - 1);
    return mp_[k - 1] + integrate_mp(t_[k - 1], t);
}

double MeridionalCurve::m_of_mprime(double mprime) const
{
    if (mprime < 0.0 || mprime > mprime_length()) {
        const bool before = mprime < 0.0;
        const Vec2 p = before ? c_.eval(c_.t0()) : c_.eval(c_.t1());
        const Vec2 t = tangent_at_m(before ? 0.0 : length());
        const double dmp = before ? -mprime : mprime - mprime_length();
        const double dr = before ? -t.y() : t.y();
        double dm;
        if (std::abs(dr * dmp) < 1e-12) dm = dmp * p.y();
        else dm = p.y() * (std::exp(dr * dmp) - 1.0) / dr;
        return before ? -dm : length() + dm;
    }
    const auto it = std::upper_bound(mp_.begin(), mp_.end(), mprime);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - mp_.begin()), 1, mp_.size() - 1);
    const double ta = t_[k - 1], tb = t_[k];
    auto f = [&](double t) { return mp_[k - 1] + integrate_mp(ta, t) - mprime; };
    double t;
    if (f(ta) >= 0.0) t = ta;
    else if (f(tb) <= 0.0) t = tb;
    else t = Math::brent_root(f, ta, tb, 1e-15);
    return m_[k - 1] + integrate_m(ta, t);
}

std::vector<Vec2> MeridionalCurve::sample(int n) const
{
    if (n < 2) throw std::invalid_argument("MeridionalCurve::sample: n >= 2");
    std::vector<Vec2> out;
    for (int i = 0; i < n; ++i) out.push_back(at_m(length() * i / (n - 1)));
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────

MeridionalChannel::MeridionalChannel(MeridionalCurve hub, MeridionalCurve shroud, int n_samples)
    : hub_(std::move(hub)), shroud_(std::move(shroud)), n_(n_samples)
{
    if (n_ < 4) throw std::invalid_argument("MeridionalChannel: n_samples >= 4");
    const Vec2 h0 = hub_.at_m(0.0), s0 = shroud_.at_m(0.0);
    if (!(s0.y() > h0.y())) throw std::invalid_argument("MeridionalChannel: shroud must lie above the hub");
}

MeridionalCurve MeridionalChannel::streamline(double span) const
{
    std::vector<Vec2> pts;
    const auto h = hub_.sample(n_), s = shroud_.sample(n_);
    for (int i = 0; i < n_; ++i) pts.push_back((1.0 - span) * h[static_cast<std::size_t>(i)] + span * s[static_cast<std::size_t>(i)]);
    return MeridionalCurve::through_points(pts);
}

MeridionalChannel MeridionalChannel::straight(double x0, double x1, double r_hub0, double r_hub1, double r_tip0,
                                              double r_tip1)
{
    return MeridionalChannel(MeridionalCurve::line({x0, r_hub0}, {x1, r_hub1}),
                             MeridionalCurve::line({x0, r_tip0}, {x1, r_tip1}));
}

}  // namespace BladeGen::Blade
