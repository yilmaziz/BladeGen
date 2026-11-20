#include "fit/SectionFit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace BladeGen::Fit {

Blade::SectionPlacement placement_from_frame(const Normalized& f, const Blade::MeridionalCurve& sl)
{
    Blade::SectionPlacement p;
    const double L = sl.length();
    p.m_le = sl.m_of_mprime(f.le.x()) / L;
    p.m_te = sl.m_of_mprime(f.le.x() + f.chord * std::cos(f.stagger)) / L;
    p.stagger_deg = Math::rad2deg(f.stagger);
    p.theta_le_deg = Math::rad2deg(f.le.y());
    p.mirror = f.mirrored;
    return p;
}

std::vector<double> polyline_distances(const std::vector<Vec3>& pts, const std::vector<Vec3>& poly)
{
    std::vector<double> d(pts.size(), std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t k = 0; k + 1 < poly.size(); ++k) {
            const Vec3 a = poly[k], ab = poly[k + 1] - a;
            const double L2 = ab.squaredNorm();
            const double t = L2 > 0.0 ? std::clamp((pts[i] - a).dot(ab) / L2, 0.0, 1.0) : 0.0;
            d[i] = std::min(d[i], (pts[i] - (a + t * ab)).norm());
        }
    return d;
}

SectionFitResult fit_section(const std::vector<Vec2>& loop, const Blade::MeridionalCurve& sl, double span,
                             const Profile::ProfileParams& start, const FitOptions& opt)
{
    SectionFitResult r;
    r.frame = normalize_contour(loop);
    r.fit = fit_profile(r.frame.ss, r.frame.ps, start, opt);
    r.station.span = span;
    r.station.profile = r.fit.params;
    r.station.placement = placement_from_frame(r.frame, sl);
    r.station.streamline = sl;
    for (const auto& q : loop) {
        const Vec2 xr = sl.at_m(sl.m_of_mprime(q.x()));
        r.published_xyz.push_back(Blade::cartesian(xr.x(), xr.y(), q.y()));
    }
    const auto sec = Blade::build_section(r.fit.params, r.station.placement, sl, span, 120, 400);
    r.fitted_xyz = sec.xyz;
    const auto d = polyline_distances(r.published_xyz, r.fitted_xyz);
    double s2 = 0.0;
    for (double v : d) {
        s2 += v * v;
        r.max_3d = std::max(r.max_3d, v);
    }
    r.rms_3d = std::sqrt(s2 / static_cast<double>(d.size()));
    return r;
}

}  // namespace BladeGen::Fit
