#include "profile/Quality.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace BladeGen::Profile {

std::string to_string(Severity s)
{
    switch (s) {
        case Severity::Info: return "INFO";
        case Severity::Warning: return "WARNING";
        case Severity::Error: return "ERROR";
    }
    return "?";
}

int QualityReport::count(Severity s) const
{
    return static_cast<int>(std::count_if(issues.begin(), issues.end(),
                                          [s](const Issue& i) { return i.severity == s; }));
}

std::string QualityReport::summary() const
{
    std::ostringstream o;
    const int e = count(Severity::Error), w = count(Severity::Warning);
    if (e == 0 && w == 0) return "geometry OK";
    o << e << (e == 1 ? " error" : " errors") << ", " << w << (w == 1 ? " warning" : " warnings");
    for (const auto& i : issues)
        if (i.severity != Severity::Info) o << "\n  [" << to_string(i.severity) << "] " << i.message;
    return o.str();
}

namespace {

bool segments_cross(const Vec2& a1, const Vec2& a2, const Vec2& b1, const Vec2& b2)
{
    auto orient = [](const Vec2& p, const Vec2& q, const Vec2& r) { return Math::cross2(q - p, r - p); };
    const double d1 = orient(b1, b2, a1), d2 = orient(b1, b2, a2);
    const double d3 = orient(a1, a2, b1), d4 = orient(a1, a2, b2);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0;
}

/// Turning points of a sequence with hysteresis h (ripples smaller than h ignored).
std::vector<std::size_t> extrema(const std::vector<double>& v, double h)
{
    std::vector<std::size_t> out;
    if (v.size() < 3) return out;
    int dir = 0;  // +1 rising, -1 falling
    std::size_t ext = 0;
    for (std::size_t i = 1; i < v.size(); ++i) {
        if (dir >= 0 && v[i] > v[ext]) { ext = i; if (dir == 0 && v[i] - v[0] > h) dir = 1; continue; }
        if (dir <= 0 && v[i] < v[ext]) { ext = i; if (dir == 0 && v[0] - v[i] > h) dir = -1; continue; }
        if (dir == 1 && v[ext] - v[i] > h) { out.push_back(ext); dir = -1; ext = i; }
        else if (dir == -1 && v[i] - v[ext] > h) { out.push_back(ext); dir = 1; ext = i; }
        else if (dir == 0) {
            if (v[i] - v[0] > h) { dir = 1; ext = i; }
            else if (v[0] - v[i] > h) { dir = -1; ext = i; }
        }
    }
    return out;
}

void side_checks(QualityReport& rep, const CurveSample& c, int side, const QualityOptions& opt, double le_kmax)
{
    const std::string name = side == 0 ? "suction side" : "pressure side";
    const int body = side == 0 ? static_cast<int>(Piece::SuctionSide) : static_cast<int>(Piece::PressureSide);
    const double L = c.s.back();
    auto at = [&](std::size_t i, Severity sv, const std::string& code, const std::string& msg) {
        rep.issues.push_back({sv, code, msg, c.xy[i], side, L > 0 ? c.s[i] / L : 0.0});
    };
    // body: inflections and curvature oscillation
    std::vector<double> kb;
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < c.kappa.size(); ++i)
        if (c.piece[i] == body) { kb.push_back(c.kappa[i]); idx.push_back(i); }
    if (kb.size() < 3) return;
    const auto [mn, mx] = std::minmax_element(kb.begin(), kb.end());
    const double range = *mx - *mn;
    const double thr = 1e-3 * std::max(std::abs(*mn), std::abs(*mx));
    int infl = 0;
    double last = 0.0;
    for (std::size_t j = 0; j < kb.size(); ++j) {
        if (std::abs(kb[j]) <= thr) continue;
        if (last != 0.0 && (kb[j] > 0) != (last > 0)) {
            ++infl;
            const Severity sv = side == 0 ? Severity::Warning : Severity::Info;
            at(idx[j], sv, side == 0 ? "ss_inflection" : "ps_inflection",
               "inflection point on the " + name + " (curvature changes sign)");
        }
        last = kb[j];
    }
    if (side == 1 && infl > 1)
        at(idx[kb.size() / 2], Severity::Warning, "ps_waviness",
           "pressure side has " + std::to_string(infl) + " inflection points (wavy surface)");
    const auto ex = extrema(kb, opt.wiggle_rel * range);
    if (static_cast<int>(ex.size()) > opt.max_extrema_body)
        at(idx[ex[ex.size() / 2]], Severity::Warning, "curvature_oscillation",
           name + " curvature oscillates (" + std::to_string(ex.size()) +
               " extrema): the surface is not fair; check the handle lengths and CP2");
    const double kbody = std::max(std::abs(*mn), std::abs(*mx));
    if (le_kmax > 0.0 && kbody > le_kmax) {
        const auto j = static_cast<std::size_t>(std::distance(kb.begin(), std::abs(*mx) >= std::abs(*mn) ? mx : mn));
        at(idx[j], Severity::Warning, "body_curvature_spike",
           name + " body curvature exceeds the leading-edge curvature (local kink)");
    }
    if (side == 0) { rep.ss_inflections = infl; rep.ss_extrema = static_cast<int>(ex.size()); }
    else { rep.ps_inflections = infl; rep.ps_extrema = static_cast<int>(ex.size()); }
}

/// The LE / TE connector must rise from one junction curvature to a single peak
/// near the edge point and fall to the other: a dip below both junction values is
/// an S-bend (overshoot / loop), a second peak a bump.
void connector_checks(QualityReport& rep, const Profile2D& prof, int k, const QualityOptions& opt)
{
    const std::string edge = k == 0 ? "LE" : "TE";
    const double ta = prof.breaks()[static_cast<std::size_t>(k)];
    const double tb = prof.breaks()[static_cast<std::size_t>(k + 1)];
    const int n = 4 * opt.n_edge + 1;
    std::vector<double> kap(static_cast<std::size_t>(n));
    std::vector<double> tt(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double t = ta + (tb - ta) * i / (n - 1);
        tt[static_cast<std::size_t>(i)] = t;
        kap[static_cast<std::size_t>(i)] = prof.curvature(t, i == n - 1);
    }
    const double kend = std::min(kap.front(), kap.back());
    const auto imin = static_cast<std::size_t>(std::distance(kap.begin(), std::min_element(kap.begin(), kap.end())));
    const double kpeak = *std::max_element(kap.begin(), kap.end());
    const double range = kpeak - kend;
    auto add = [&](std::size_t i, Severity sv, const std::string& code, const std::string& msg) {
        const double t = tt[i];
        const bool ss_half = (k == 0) ? t > prof.t_le_point() : t < prof.t_te_point();
        rep.issues.push_back({sv, code, msg, prof.point(t), ss_half ? 0 : 1, 0.0});
    };
    if (kap[imin] < kend - 0.02 * range) {
        if (kap[imin] < 0.0)
            add(imin, Severity::Error, "edge_inflection",
                edge + " connector is S-shaped (curvature dips to " + std::to_string(kap[imin]) +
                    "): the edge overshoots or loops; change the " + edge + " radius, stretch, wedge angles or tension");
        else
            add(imin, Severity::Warning, "edge_curvature_dip",
                edge + " connector curvature dips below its junction values (flat spot on the edge)");
    }
    const auto ex = extrema(kap, 0.05 * range);
    int maxima = 0;
    for (std::size_t e : ex) if (e > 0 && e + 1 < kap.size() && kap[e] >= kap[e - 1] && kap[e] >= kap[e + 1]) ++maxima;
    if (maxima > 1)
        add(ex.front(), Severity::Warning, "edge_double_peak",
            edge + " connector curvature has " + std::to_string(maxima) + " peaks (bumpy edge)");
}

}  // namespace

int count_self_intersections(const std::vector<Vec2>& closed)
{
    const std::size_t n = closed.size() < 2 ? 0 : closed.size() - 1;
    int cnt = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;
            if (segments_cross(closed[i], closed[i + 1], closed[j], closed[j + 1])) ++cnt;
        }
    return cnt;
}

QualityReport check_quality(const Profile2D& prof, const QualityOptions& opt)
{
    QualityReport rep;
    rep.ss = prof.sample_side(Side::Suction, opt.n_edge, opt.n_body);
    rep.ps = prof.sample_side(Side::Pressure, opt.n_edge, opt.n_body);
    rep.junctions = prof.junctions();

    double kmax = 0.0, le_kmax = 0.0;
    for (const auto* c : {&rep.ss, &rep.ps})
        for (std::size_t i = 0; i < c->kappa.size(); ++i) {
            kmax = std::max(kmax, std::abs(c->kappa[i]));
            if (c->piece[i] == static_cast<int>(Piece::LeadingEdge)) le_kmax = std::max(le_kmax, c->kappa[i]);
        }
    rep.max_curvature = kmax;

    // continuity at the four junctions
    const bool c2 = prof.params().connector == ConnectorMode::C2Parametric;
    for (const auto& j : rep.junctions) {
        const Vec2 where = prof.point(j.t);
        if (j.gap > 1e-9)
            rep.issues.push_back({Severity::Error, "junction_gap", j.name + ": position gap " + std::to_string(j.gap), where});
        if (j.tangent_angle > opt.junction_angle)
            rep.issues.push_back({Severity::Error, "junction_kink",
                                  j.name + ": tangent kink " + std::to_string(Math::rad2deg(j.tangent_angle)) + " deg", where});
        if (j.curvature_jump > opt.junction_kappa_rel * std::max(kmax, 1.0))
            rep.issues.push_back({Severity::Error, "junction_curvature_jump",
                                  j.name + ": curvature jump " + std::to_string(j.curvature_jump) + " (not C2)", where});
        if (c2 && (j.d1_jump > 1e-7 || j.d2_jump > 1e-6))
            rep.issues.push_back({Severity::Error, "junction_parametric",
                                  j.name + ": parametric C2 violated (d1 " + std::to_string(j.d1_jump) +
                                      ", d2 " + std::to_string(j.d2_jump) + ")", where});
    }

    connector_checks(rep, prof, 0, opt);
    connector_checks(rep, prof, 2, opt);
    side_checks(rep, rep.ss, 0, opt, le_kmax);
    side_checks(rep, rep.ps, 1, opt, le_kmax);

    // closed loop: orientation and self-intersection
    const auto loop = prof.sample_loop(opt.n_edge / 2 + 2, opt.n_body / 2 + 2);
    double area = 0.0;
    for (std::size_t i = 0; i + 1 < loop.xy.size(); ++i) area += Math::cross2(loop.xy[i], loop.xy[i + 1]);
    area *= 0.5;
    if (!(area < 0.0))
        rep.issues.push_back({Severity::Error, "orientation",
                              "profile loop is not clockwise: the pressure side lies above the suction side", Vec2(0.5, 0.0)});
    rep.self_intersections = count_self_intersections(loop.xy);
    if (rep.self_intersections > 0) {
        Vec2 where(0.5, 0.0);
        const std::size_t n = loop.xy.size() - 1;
        for (std::size_t i = 0; i < n && where == Vec2(0.5, 0.0); ++i)
            for (std::size_t j = i + 2; j < n; ++j)
                if (!(i == 0 && j == n - 1) && segments_cross(loop.xy[i], loop.xy[i + 1], loop.xy[j], loop.xy[j + 1])) {
                    where = loop.xy[i];
                    break;
                }
        rep.issues.push_back({Severity::Error, "self_intersection",
                              "profile self-intersects (" + std::to_string(rep.self_intersections) + " crossings)", where});
    }

    // thickness
    const auto th = prof.thickness(200);
    rep.min_thickness = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < th.xi.size(); ++i) {
        if (th.xi[i] < 0.03 || th.xi[i] > 0.97) continue;
        const double t = th.measured[i];
        if (!std::isfinite(t) || t <= opt.min_thickness_frac) {
            rep.issues.push_back({Severity::Error, "negative_thickness",
                                  "pressure side touches or crosses the suction side near m' = " + std::to_string(th.mprime[i]),
                                  prof.suction_side().eval(static_cast<double>(i) / (static_cast<double>(th.xi.size()) - 1.0)), 0, th.xi[i]});
            break;
        }
        rep.min_thickness = std::min(rep.min_thickness, t);
    }
    for (std::size_t i = 0; i < prof.xi_ctrl().size(); ++i)
        if (thickness_law(prof.xi_ctrl()[i], prof.params().thickness, prof.xi_max()) <= 0.0)
            rep.issues.push_back({Severity::Error, "thickness_law",
                                  "thickness law is <= 0 at PS station " + std::to_string(i), prof.ps_free_points()[i]});

    // m' must increase along both sides (no fold back)
    for (int sd = 0; sd < 2; ++sd) {
        const auto& c = sd == 0 ? prof.suction_side() : prof.pressure_side();
        double prev = c.eval(0.0).x();
        for (int i = 1; i <= 400; ++i) {
            const Vec2 q = c.eval(i / 400.0);
            if (q.x() < prev - 1e-12) {
                rep.issues.push_back({Severity::Warning, "mprime_fold",
                                      std::string(sd == 0 ? "suction" : "pressure") + " side folds back in m'", q, sd});
                break;
            }
            prev = q.x();
        }
    }
    std::stable_sort(rep.issues.begin(), rep.issues.end(),
                     [](const Issue& a, const Issue& b) { return a.severity > b.severity; });
    return rep;
}

}  // namespace BladeGen::Profile
