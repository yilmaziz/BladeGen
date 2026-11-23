#include "blade/FlowPath.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Blade {

namespace {

constexpr int kDegree = 5;
constexpr int kTablePerPiece = 128;  // table intervals per u piece

void require(bool ok, const std::string& msg)
{
    if (!ok) throw std::invalid_argument("MeridionalConfig: " + msg);
}

/// Five-point Gauss-Legendre nodes and weights on [-1, 1].
constexpr std::array<double, 5> kGx = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                                       0.9061798459386640};
constexpr std::array<double, 5> kGw = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                       0.4786286704993665, 0.2369268850561891};

std::vector<double> wall_knots()
{
    std::vector<double> U;
    for (int i = 0; i <= kDegree; ++i) U.push_back(0.0);
    for (int i = 0; i < kDegree; ++i) U.push_back(1.0);
    for (int i = 0; i < kDegree; ++i) U.push_back(2.0);
    for (int i = 0; i <= kDegree; ++i) U.push_back(3.0);
    return U;
}

/// Cubic Hermite on [0, h] with values f0, f1 and slopes d0, d1 at local coordinate t.
template <class T>
T hermite(const T& f0, const T& f1, const T& d0, const T& d1, double h, double t)
{
    const double s = t / h, s2 = s * s, s3 = s2 * s;
    return (2 * s3 - 3 * s2 + 1) * f0 + (s3 - 2 * s2 + s) * h * d0 + (-2 * s3 + 3 * s2) * f1 + (s3 - s2) * h * d1;
}

/// Its derivative d/dt.
template <class T>
T hermite_d(const T& f0, const T& f1, const T& d0, const T& d1, double h, double t)
{
    const double s = t / h, s2 = s * s;
    return ((6 * s2 - 6 * s) * f0 + (-6 * s2 + 6 * s) * f1) / h + (3 * s2 - 4 * s + 1) * d0 + (3 * s2 - 2 * s) * d1;
}

}  // namespace

std::string to_string(WallLaw w)
{
    return w == WallLaw::Linear ? "linear" : "bell_mehta";
}

WallLaw wall_law_from_string(const std::string& s)
{
    if (s == "linear" || s == "Linear") return WallLaw::Linear;
    if (s == "bell_mehta" || s == "Bell & Mehta" || s == "bell-mehta") return WallLaw::BellMehta;
    throw std::invalid_argument("unknown wall law '" + s + "' (linear | bell_mehta)");
}

double bell_mehta(double xi) { return xi * xi * xi * (10.0 + xi * (-15.0 + 6.0 * xi)); }
double bell_mehta_d1(double xi) { return 30.0 * xi * xi * (1.0 - xi) * (1.0 - xi); }
double bell_mehta_d2(double xi) { return 60.0 * xi * (1.0 - xi) * (1.0 - 2.0 * xi); }

void MeridionalConfig::validate() const
{
    for (const auto* w : {&hub, &casing}) {
        const std::string n = w == &hub ? "hub" : "casing";
        require(std::isfinite(w->z_le) && std::isfinite(w->z_te) && std::isfinite(w->r_le) && std::isfinite(w->r_te),
                n + " coordinates must be finite");
        require(w->z_te > w->z_le, n + ": z_TE must be larger than z_LE");
        require(w->r_le > 0.0 && w->r_te > 0.0, n + ": radii must be > 0");
    }
    require(casing.r_le > hub.r_le && casing.r_te > hub.r_te, "the casing must lie above the hub at the LE and the TE");
    require(inlet_extension >= 0.05 && inlet_extension <= 5.0, "inlet extension must be in [0.05, 5]");
    require(outlet_extension >= 0.05 && outlet_extension <= 5.0, "outlet extension must be in [0.05, 5]");
    // the linear walls continue with their slope: radius must stay positive on the extensions
    for (const auto* w : {&hub, &casing}) {
        const WallLaw law = w == &hub ? hub_law : casing_law;
        if (law != WallLaw::Linear) continue;
        const double dz = w->z_te - w->z_le, k = (w->r_te - w->r_le) / dz;
        require(w->r_le - k * inlet_extension * dz > 0.0 && w->r_te + k * outlet_extension * dz > 0.0,
                std::string(w == &hub ? "hub" : "casing") + ": the linear wall reaches r <= 0 on an extension");
    }
}

Math::BSpline2 wall_spline(const WallPoints& w, WallLaw law, double inlet_ext, double outlet_ext)
{
    const double dz = w.z_te - w.z_le, dr = w.r_te - w.r_le;
    if (!(dz > 0.0)) throw std::invalid_argument("wall_spline: z_TE must be larger than z_LE");
    const double k = law == WallLaw::Linear ? dr / dz : 0.0;  // extension slope dr/dz
    std::vector<Vec2> P;
    P.reserve(16);
    // inlet piece: u in [0, 1]
    const double z0 = w.z_le - inlet_ext * dz;
    for (int i = 0; i <= kDegree; ++i) {
        const double z = z0 + (w.z_le - z0) * i / kDegree;
        P.emplace_back(z, w.r_le + k * (z - w.z_le));
    }
    // bladed piece: u in [1, 2]
    static constexpr std::array<double, 6> kBellMehtaBernstein = {0, 0, 0, 1, 1, 1};
    for (int i = 1; i <= kDegree; ++i) {
        const double z = w.z_le + dz * i / kDegree;
        const double f = law == WallLaw::Linear ? static_cast<double>(i) / kDegree : kBellMehtaBernstein[static_cast<std::size_t>(i)];
        P.emplace_back(z, w.r_le + dr * f);
    }
    // outlet piece: u in [2, 3]
    const double z3 = w.z_te + outlet_ext * dz;
    for (int i = 1; i <= kDegree; ++i) {
        const double z = w.z_te + (z3 - w.z_te) * i / kDegree;
        P.emplace_back(z, w.r_te + k * (z - w.z_te));
    }
    return Math::BSpline2(kDegree, wall_knots(), std::move(P));
}

// ─────────────────────────────────────────────────────────────────────────────

StreamSurface::StreamSurface(const Math::BSpline2& sl, double span) : c_(sl), span_(span)
{
    if (!(c_.t0() == 0.0 && c_.t1() == 3.0)) throw std::invalid_argument("StreamSurface: streamline must have u in [0, 3]");
    double rmin = c_.ctrl().front().y(), rmax = rmin;
    for (const auto& p : c_.ctrl()) {
        rmin = std::min(rmin, p.y());
        rmax = std::max(rmax, p.y());
    }
    if (!(rmin > 0.0)) throw std::invalid_argument("StreamSurface: radius must stay > 0");
    cylindrical_ = (rmax - rmin) <= 1e-12 * rmax;
    r0_ = cylindrical_ ? rmax : 0.0;

    const int n = 3 * kTablePerPiece;
    u_.resize(static_cast<std::size_t>(n + 1));
    for (int i = 0; i <= n; ++i) u_[static_cast<std::size_t>(i)] = 3.0 * i / n;
    m_.assign(u_.size(), 0.0);
    mp_.assign(u_.size(), 0.0);
    xr_.resize(u_.size());
    dxr_.resize(u_.size());
    for (std::size_t i = 0; i < u_.size(); ++i) {
        const auto d = c_.derivs(u_[i], 1);
        xr_[i] = d[0];
        const double sp = d[1].norm();
        if (!(sp > 0.0)) throw std::invalid_argument("StreamSurface: zero parametric speed");
        dxr_[i] = d[0].y() * d[1] / sp;  // d(x, r)/dm' = r * unit tangent
        if (i == 0) continue;
        const double a = u_[i - 1], b = u_[i], h = 0.5 * (b - a), c = 0.5 * (a + b);
        double dm = 0.0, dmp = 0.0;
        for (std::size_t g = 0; g < kGx.size(); ++g) {
            const auto q = c_.derivs(c + h * kGx[g], 1);
            const double v = q[1].norm();
            dm += kGw[g] * h * v;
            dmp += kGw[g] * h * v / q[0].y();
        }
        m_[i] = m_[i - 1] + dm;
        mp_[i] = mp_[i - 1] + dmp;
    }
    const auto ile = static_cast<std::size_t>(kTablePerPiece), ite = static_cast<std::size_t>(2 * kTablePerPiece);
    const double m0 = m_[ile], mp0 = mp_[ile];
    for (std::size_t i = 0; i < u_.size(); ++i) {
        m_[i] -= m0;
        mp_[i] -= mp0;
    }
    mp_te_ = mp_[ite];
    m_te_ = m_[ite];
    le_ = xr_[ile];
    te_ = xr_[ite];
    x_le_ = le_.x();
}

Vec2 StreamSurface::xr(double mp) const
{
    if (cylindrical_) return {x_le_ + r0_ * mp, r0_};
    if (mp <= mp_.front()) return xr_.front() + (mp - mp_.front()) * dxr_.front();
    if (mp >= mp_.back()) return xr_.back() + (mp - mp_.back()) * dxr_.back();
    const auto it = std::upper_bound(mp_.begin(), mp_.end(), mp);
    const auto k = static_cast<std::size_t>(it - mp_.begin());
    const double h = mp_[k] - mp_[k - 1];
    return hermite<Vec2>(xr_[k - 1], xr_[k], dxr_[k - 1], dxr_[k], h, mp - mp_[k - 1]);
}

Vec2 StreamSurface::dxr(double mp) const
{
    if (cylindrical_) return {r0_, 0.0};
    if (mp <= mp_.front()) return dxr_.front();
    if (mp >= mp_.back()) return dxr_.back();
    const auto it = std::upper_bound(mp_.begin(), mp_.end(), mp);
    const auto k = static_cast<std::size_t>(it - mp_.begin());
    const double h = mp_[k] - mp_[k - 1];
    return hermite_d<Vec2>(xr_[k - 1], xr_[k], dxr_[k - 1], dxr_[k], h, mp - mp_[k - 1]);
}

double StreamSurface::m_of_mprime(double mp) const
{
    if (cylindrical_) return r0_ * mp;
    if (mp <= mp_.front()) return m_.front() + (mp - mp_.front()) * xr_.front().y();
    if (mp >= mp_.back()) return m_.back() + (mp - mp_.back()) * xr_.back().y();
    const auto it = std::upper_bound(mp_.begin(), mp_.end(), mp);
    const auto k = static_cast<std::size_t>(it - mp_.begin());
    const double h = mp_[k] - mp_[k - 1];
    return hermite<double>(m_[k - 1], m_[k], xr_[k - 1].y(), xr_[k].y(), h, mp - mp_[k - 1]);  // dm/dm' = r
}

Vec3 StreamSurface::cartesian(double mp, double theta) const
{
    const Vec2 p = xr(mp);
    return {p.x(), p.y() * std::sin(theta), p.y() * std::cos(theta)};
}

std::vector<Vec2> StreamSurface::polyline(int n_per_piece) const
{
    if (n_per_piece < 1) throw std::invalid_argument("StreamSurface::polyline: n_per_piece >= 1");
    std::vector<Vec2> out;
    const int n = 3 * n_per_piece;
    for (int i = 0; i <= n; ++i) out.push_back(c_.eval(3.0 * i / n));
    return out;
}

double StreamSurface::curvature(double mp) const
{
    const double mpc = std::clamp(mp, mp_.front(), mp_.back());
    const auto it = std::upper_bound(mp_.begin(), mp_.end(), mpc);
    const auto k = std::clamp<std::size_t>(static_cast<std::size_t>(it - mp_.begin()), 1, mp_.size() - 1);
    const double w = (mpc - mp_[k - 1]) / (mp_[k] - mp_[k - 1]);
    const double u = u_[k - 1] + w * (u_[k] - u_[k - 1]);
    const auto d = c_.derivs(u, 2);
    return Math::signed_curvature(d[1], d[2]);
}

// ─────────────────────────────────────────────────────────────────────────────

FlowPath::FlowPath(MeridionalConfig cfg) : cfg_(cfg)
{
    cfg_.validate();
    hub_ = wall_spline(cfg_.hub, cfg_.hub_law, cfg_.inlet_extension, cfg_.outlet_extension);
    casing_ = wall_spline(cfg_.casing, cfg_.casing_law, cfg_.inlet_extension, cfg_.outlet_extension);
}

Math::BSpline2 FlowPath::streamline(double span) const
{
    if (!std::isfinite(span)) throw std::invalid_argument("FlowPath::streamline: non-finite span");
    std::vector<Vec2> P(hub_.ctrl().size());
    for (std::size_t i = 0; i < P.size(); ++i) P[i] = (1.0 - span) * hub_.ctrl()[i] + span * casing_.ctrl()[i];
    return Math::BSpline2(hub_.degree(), hub_.knots(), std::move(P));
}

Vec2 FlowPath::le_point(double span) const
{
    return (1.0 - span) * Vec2(cfg_.hub.z_le, cfg_.hub.r_le) + span * Vec2(cfg_.casing.z_le, cfg_.casing.r_le);
}

Vec2 FlowPath::te_point(double span) const
{
    return (1.0 - span) * Vec2(cfg_.hub.z_te, cfg_.hub.r_te) + span * Vec2(cfg_.casing.z_te, cfg_.casing.r_te);
}

}  // namespace BladeGen::Blade
