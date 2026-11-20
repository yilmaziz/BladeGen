#include "blade/Blade.hpp"

#include "math/Interpolation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace BladeGen::Blade {

using Math::deg2rad;
using Profile::ProfileParams;

namespace {

void require(bool ok, const std::string& msg)
{
    if (!ok) throw std::invalid_argument("BladeDefinition: " + msg);
}

/// Natural cubic spline through (x_i, y_i) evaluated at xq (xq clamped to the range).
double natural_spline(const std::vector<double>& x, const std::vector<double>& y, double xq)
{
    xq = std::clamp(xq, x.front(), x.back());
    if (x.size() < 3) {
        if (x.size() == 1) return y.front();
        const double w = (xq - x[0]) / (x[1] - x[0]);
        return y[0] + w * (y[1] - y[0]);
    }
    std::vector<Math::Vec<1>> pts;
    for (double v : y) pts.push_back(Math::Vec<1>(v));
    const auto c = Math::interpolate_cubic<1>(x, pts, Math::EndCondition<1>::natural(), Math::EndCondition<1>::natural());
    return c.eval(xq)[0];
}

double linear(const std::vector<double>& x, const std::vector<double>& y, double xq)
{
    if (x.size() == 1) return y.front();
    xq = std::clamp(xq, x.front(), x.back());
    const auto it = std::upper_bound(x.begin(), x.end(), xq);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - x.begin()), 1, x.size() - 1);
    const double w = (xq - x[k - 1]) / (x[k] - x[k - 1]);
    return y[k - 1] + w * (y[k] - y[k - 1]);
}

/// Pointwise blend of two streamlines at equal fractions of their lengths.
MeridionalCurve blend_streamlines(const MeridionalCurve& a, const MeridionalCurve& b, double w, int n = 161)
{
    const auto pa = a.sample(n), pb = b.sample(n);
    std::vector<Vec2> p(static_cast<std::size_t>(n));
    for (std::size_t i = 0; i < p.size(); ++i) p[i] = (1.0 - w) * pa[i] + w * pb[i];
    return MeridionalCurve::through_points(p);
}

}  // namespace

void BladeDefinition::validate() const
{
    require(!stations.empty(), "need at least one station");
    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto& s = stations[i];
        require(std::isfinite(s.span), "non-finite span");
        if (i > 0) require(s.span > stations[i - 1].span, "station spans must increase");
        require(s.placement.m_te > s.placement.m_le, "m_te must be larger than m_le");
        require(std::abs(s.placement.stagger_deg) < 85.0, "|stagger| must be below 85 deg");
        Profile::validate(s.profile);
    }
    const bool explicit_sl = std::all_of(stations.begin(), stations.end(),
                                         [](const BladeStation& s) { return s.streamline.has_value(); });
    const bool any_sl = std::any_of(stations.begin(), stations.end(),
                                    [](const BladeStation& s) { return s.streamline.has_value(); });
    require(explicit_sl || !any_sl, "either every station has its own streamline or none has");
    require(explicit_sl || !channel.hub().empty(), "a channel is needed when the stations have no streamline");
    for (const auto& s : stations)
        require(s.placement.mirror == stations.front().placement.mirror, "all stations must have the same mirror flag");
    require(n_blades >= 1, "n_blades >= 1");
    require(n_sections >= 2, "n_sections >= 2");
    require(span_extension >= 0.0 && span_extension < 0.5, "span_extension must be in [0, 0.5)");
    require(tip_clearance >= 0.0, "tip_clearance must be >= 0");
}

BladeDefinition BladeDefinition::demo_rotor()
{
    BladeDefinition d;
    d.name = "demo_rotor";
    d.channel = MeridionalChannel::straight(0.0, 0.12, 0.25, 0.27, 0.40, 0.395);
    d.n_blades = 24;
    d.n_sections = 11;
    auto base = ProfileParams::naca65_1210();
    base.te = {0.003, 1.2, 0.0, 10.0, 6.4, 4.0, 1.0};
    base.ps.handle_le_frac = 0.06;
    auto hub = base, mid = base, tip = base;
    hub.ss.ctrl2_theta = 0.20;  // more camber at the hub
    hub.le.beta_deg = 20.0;
    hub.thickness.t_max = 0.10;
    tip.ss.ctrl2_theta = 0.14;
    tip.le.beta_deg = 11.0;
    tip.thickness.t_max = 0.07;
    tip.thickness.t_le = 0.025;
    d.stations = {
        {0.0, hub, {0.30, 0.70, 30.0, 0.0, false}, std::nullopt},
        {0.5, mid, {0.30, 0.70, 42.0, 0.0, false}, std::nullopt},
        {1.0, tip, {0.30, 0.70, 52.0, 0.0, false}, std::nullopt},
    };
    return d;
}

Vec2 SectionMap::conformal(const Vec2& q) const
{
    const double c = std::cos(gamma), s = std::sin(gamma);
    return {mprime_le + chord * (q.x() * c - sign * q.y() * s), theta_le + chord * (q.x() * s + sign * q.y() * c)};
}

Vec3 SectionMap::cartesian(const Vec2& mt) const
{
    const Vec2 xr = streamline->at_m(streamline->m_of_mprime(mt.x()));
    return ::BladeGen::Blade::cartesian(xr.x(), xr.y(), mt.y());
}

SectionMap make_section_map(const MeridionalCurve& sl, const SectionPlacement& pl)
{
    SectionMap m;
    m.streamline = &sl;
    m.gamma = deg2rad(pl.stagger_deg);
    m.theta_le = deg2rad(pl.theta_le_deg);
    m.sign = pl.mirror ? -1.0 : 1.0;
    m.mprime_le = sl.mprime_of_m(pl.m_le * sl.length());
    const double mp_te = sl.mprime_of_m(pl.m_te * sl.length());
    if (!(mp_te > m.mprime_le)) throw std::invalid_argument("SectionPlacement: TE must lie downstream of the LE");
    m.chord = (mp_te - m.mprime_le) / std::cos(m.gamma);
    return m;
}

Section3D build_section(const ProfileParams& prof, const SectionPlacement& pl, const MeridionalCurve& sl,
                        double span, int n_edge, int n_body)
{
    Section3D s;
    s.span = span;
    s.params = prof;
    s.placement = pl;
    s.streamline = sl;
    const auto map = make_section_map(s.streamline, pl);
    s.mprime_le = map.mprime_le;
    s.chord_mprime = map.chord;
    const Profile::Profile2D p(prof);
    const auto loop = p.sample_loop(n_edge, n_body);
    s.unstaggered = loop.xy;
    s.piece = loop.piece;
    s.conformal.reserve(loop.xy.size());
    s.xyz.reserve(loop.xy.size());
    for (const auto& q : loop.xy) {
        s.conformal.push_back(map.conformal(q));
        s.xyz.push_back(map.cartesian(s.conformal.back()));
    }
    return s;
}

// ─────────────────────────────────────────────────────────────────────────────

Blade::Blade(BladeDefinition def) : def_(std::move(def))
{
    def_.validate();
    for (const auto& s : def_.stations) spans_.push_back(s.span);
}

double Blade::interp(const std::vector<double>& y, double span) const
{
    return def_.interpolation == SpanInterpolation::Cubic ? natural_spline(spans_, y, span) : linear(spans_, y, span);
}

ProfileParams Blade::profile_at(double span) const
{
    const auto& st = def_.stations;
    ProfileParams out = st.front().profile;
    if (st.size() == 1) return out;
    std::vector<std::vector<double>> cols(Profile::param_table().size());
    for (const auto& s : st) {
        const auto v = Profile::to_vector(s.profile);
        for (std::size_t k = 0; k < v.size(); ++k) cols[k].push_back(v[k]);
    }
    std::vector<double> v(cols.size());
    for (std::size_t k = 0; k < cols.size(); ++k) v[k] = interp(cols[k], span);
    Profile::from_vector(out, v);
    const bool same_ps = std::all_of(st.begin(), st.end(), [&](const BladeStation& s) {
        return s.profile.ps.mprime_ctrl.size() == st.front().profile.ps.mprime_ctrl.size();
    });
    if (same_ps)
        for (std::size_t j = 0; j < out.ps.mprime_ctrl.size(); ++j) {
            std::vector<double> y;
            for (const auto& s : st) y.push_back(s.profile.ps.mprime_ctrl[j]);
            out.ps.mprime_ctrl[j] = linear(spans_, y, span);
        }
    // a discrete option follows the nearest station
    const auto nearest = std::min_element(st.begin(), st.end(), [&](const BladeStation& a, const BladeStation& b) {
        return std::abs(a.span - span) < std::abs(b.span - span);
    });
    out.thickness.enforce_c1 = nearest->profile.thickness.enforce_c1;
    out.connector = nearest->profile.connector;
    return out;
}

SectionPlacement Blade::placement_at(double span) const
{
    const auto& st = def_.stations;
    SectionPlacement out = st.front().placement;
    if (st.size() == 1) return out;
    auto col = [&](auto get) {
        std::vector<double> y;
        for (const auto& s : st) y.push_back(get(s.placement));
        return interp(y, span);
    };
    out.m_le = col([](const SectionPlacement& p) { return p.m_le; });
    out.m_te = col([](const SectionPlacement& p) { return p.m_te; });
    out.stagger_deg = col([](const SectionPlacement& p) { return p.stagger_deg; });
    out.theta_le_deg = col([](const SectionPlacement& p) { return p.theta_le_deg; });
    return out;
}

MeridionalCurve Blade::streamline_at(double span) const
{
    const auto& st = def_.stations;
    const bool explicit_sl = std::all_of(st.begin(), st.end(), [](const BladeStation& s) { return s.streamline.has_value(); });
    if (!explicit_sl) return def_.channel.streamline(span);
    for (const auto& s : st)
        if (std::abs(s.span - span) < 1e-12) return *s.streamline;
    if (st.size() == 1) return *st.front().streamline;
    std::size_t k = 1;
    while (k + 1 < st.size() && st[k].span < span) ++k;
    const double w = (span - st[k - 1].span) / (st[k].span - st[k - 1].span);
    return blend_streamlines(*st[k - 1].streamline, *st[k].streamline, w);
}

Section3D Blade::section(double span, int n_edge, int n_body) const
{
    return build_section(profile_at(span), placement_at(span), streamline_at(span), span, n_edge, n_body);
}

std::vector<Section3D> Blade::sections(int n_edge, int n_body) const
{
    const auto& st = def_.stations;
    const bool explicit_sl = std::all_of(st.begin(), st.end(), [](const BladeStation& s) { return s.streamline.has_value(); });
    std::vector<double> spans;
    if (explicit_sl) {
        spans = spans_;  // the stations themselves (e.g. published sections)
    } else {
        const double e = def_.span_extension;
        if (e > 0.0) spans.push_back(-e);
        for (int i = 0; i < def_.n_sections; ++i) spans.push_back(static_cast<double>(i) / (def_.n_sections - 1));
        if (e > 0.0) spans.push_back(1.0 + e);
    }
    std::vector<Section3D> out;
    out.reserve(spans.size());
    for (double s : spans) out.push_back(section(s, n_edge, n_body));
    return out;
}

}  // namespace BladeGen::Blade
