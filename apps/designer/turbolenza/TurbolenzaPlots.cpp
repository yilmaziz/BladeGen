#include "TurbolenzaPlots.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace Designer::Tl {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::size_t sz(int k) { return static_cast<std::size_t>(k); }

QPointF qp(const T::Point& p, double dy = 0.0) { return {p.m, p.y + dy}; }

Series line(const QString& name, QVector<QPointF> pts, const QColor& c, double width = 1.4, Qt::PenStyle st = Qt::SolidLine,
            bool legend = true)
{
    Series s;
    s.name = name;
    s.pts = std::move(pts);
    s.color = c;
    s.width = width;
    s.style = st;
    s.legend = legend;
    return s;
}

struct V3
{
    double x, y, v;
};

/// Sutherland-Hodgman on a scalar: the part of the polygon with v >= level (above) or v <= level.
void clip(const std::vector<V3>& in, double level, bool above, std::vector<V3>& out)
{
    out.clear();
    const std::size_t n = in.size();
    for (std::size_t k = 0; k < n; ++k) {
        const V3& a = in[k];
        const V3& b = in[(k + 1) % n];
        const bool ia = above ? a.v >= level : a.v <= level;
        const bool ib = above ? b.v >= level : b.v <= level;
        if (ia) out.push_back(a);
        if (ia != ib) {
            const double t = (level - a.v) / (b.v - a.v);
            out.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), level});
        }
    }
}

}  // namespace

QColor bandColor(int b, int bands) { return turbo((b + 0.5) / std::max(bands, 1)); }

QPolygonF bladeOutline(const T::PassageBoundaries& b, int k)
{
    const double sx = b.periodic_translation.m, sy = b.periodic_translation.y;
    QPolygonF out;
    for (auto it = b.lower_blade_boundary.rbegin(); it != b.lower_blade_boundary.rend(); ++it)
        out << QPointF(it->m + k * sx, it->y + k * sy);
    for (std::size_t q = 1; q < b.upper_blade_boundary.size(); ++q)
        out << QPointF(b.upper_blade_boundary[q].m - sx + k * sx, b.upper_blade_boundary[q].y - sy + k * sy);
    return out;
}

void drawProfile(PlotWidget* w, const BladeGen::Design::Section2D& s, const QString& title)
{
    QVector<QPointF> ss, ps, cam;
    for (const auto& p : s.quality.ss.xy) ss << QPointF(p.x(), p.y());
    for (const auto& p : s.quality.ps.xy) ps << QPointF(p.x(), p.y());
    for (const auto& p : s.camber.point) cam << QPointF(p.x(), p.y());
    w->setTitle(title);
    w->setLabels("x / c", "y / c");
    w->proportional = true;
    w->setFills({});
    w->setMarkers({});
    w->setSeries({line("suction side", ss, kSuction), line("pressure side", ps, kPressure),
                  line("camber line", cam, QColor(120, 120, 120), 1.0, Qt::DashLine)});
    w->commit(true);
}

void drawCascade(PlotWidget* w, const T::Cascade& c, double inlet_angle_deg)
{
    const auto& a = c.active;
    QVector<Series> series;
    QVector<Fill> fills;
    for (int k = -1; k <= 1; ++k) {
        const double dx = k * a.periodic_translation.m, dy = k * a.periodic_translation.y;
        QVector<QPointF> ss, ps;
        QPolygonF poly;
        for (const auto& p : a.suction_surface) ss << QPointF(p.m + dx, p.y + dy);
        for (const auto& p : a.pressure_surface) ps << QPointF(p.m + dx, p.y + dy);
        poly << ss;
        for (auto it = ps.rbegin(); it != ps.rend(); ++it) poly << *it;
        fills.push_back({poly, QColor(225, 225, 225)});
        series << line("suction side", ss, kSuction, 1.4, Qt::SolidLine, k == 0) << line("pressure side", ps, kPressure, 1.4, Qt::SolidLine, k == 0);
    }
    series << line(QString("throat o = %1 c").arg(a.throat_opening, 0, 'f', 4), {qp(a.throat_point_reference), qp(a.throat_point_neighbor)},
                   QColor(0, 150, 0), 2.2);
    const double b1 = inlet_angle_deg * std::numbers::pi / 180.0;
    const QPointF le = qp(a.leading_edge_point), dir(std::cos(b1), std::sin(b1));
    series << line("inlet flow direction", {le - 0.75 * dir, le - 0.12 * dir}, QColor(90, 90, 90), 1.6, Qt::DashLine);
    PointMarker m;
    m.p = le - 0.12 * dir;
    m.shape = 'o';
    m.color = QColor(90, 90, 90);
    m.size = 3;
    w->setTitle(QString("Cascade: stagger %1 deg, solidity %2, throat o/s %3")
                    .arg(c.stagger * 180.0 / std::numbers::pi, 0, 'f', 2)
                    .arg(c.solidity, 0, 'f', 3)
                    .arg(c.throat_over_pitch, 0, 'f', 4));
    w->setLabels("m / c (axial)", "y / c (tangential)");
    w->proportional = true;
    w->setFills(fills);
    w->setMarkers({m});
    w->setSeries(series);
    w->commit(true);
}

void drawGrid(PlotWidget* w, const T::Grid& g)
{
    const int nm = g.n_m(), ne = g.n_eta();
    QVector<QPointF> st, eta, edge;
    for (int i = 0; i < nm; ++i) {
        auto& dst = (i == g.first_blade_i || i == g.last_blade_i) ? edge : st;
        for (int j = 0; j < ne; ++j) dst << QPointF(g.node_m(i, j), g.node_y(i, j));
        dst << QPointF(kNaN, kNaN);
    }
    for (int j = 0; j < ne; ++j) {
        for (int i = 0; i < nm; ++i) eta << QPointF(g.node_m(i, j), g.node_y(i, j));
        eta << QPointF(kNaN, kNaN);
    }
    QVector<Fill> fills;
    QVector<Series> series = {line("stations (constant m)", st, QColor(70, 110, 170), 0.8),
                              line("eta lines", eta, QColor(40, 40, 40), 0.8),
                              line("first / last blade station", edge, QColor(220, 0, 0), 2.0)};
    for (int k = 0; k <= 1; ++k) {
        const QPolygonF b = bladeOutline(g.boundaries, k);
        fills.push_back({b, QColor(215, 215, 215)});
        QVector<QPointF> o(b.begin(), b.end());
        o << b.front();
        series << line("blade", o, QColor(50, 50, 50), 1.0, Qt::SolidLine, k == 0);
    }
    w->setTitle(QString("Grid %1 x %2 (stations x eta), blade stations %3 .. %4").arg(nm).arg(ne).arg(g.first_blade_i).arg(g.last_blade_i));
    w->setLabels("m / c (axial)", "y / c (tangential)");
    w->proportional = true;
    w->setFills(fills);
    w->setMarkers({});
    w->setSeries(series);
    w->commit(true);
}

void drawResiduals(PlotWidget* w, const std::vector<T::PassRecord>& records, double tol_rho, double tol_w, double tol_k,
                   const QString& status)
{
    QVector<QPointF> rho, vel, kut, lin;
    double x_max = 1.0;
    for (const auto& r : records) {
        const double x = r.index;
        x_max = std::max(x_max, x);
        rho << QPointF(x, r.density_change);
        vel << QPointF(x, r.velocity_change);
        kut << QPointF(x, r.kutta_normalized);
        lin << QPointF(x, r.linear_residual);
    }
    const QColor cr(31, 119, 180), cw(255, 127, 14), ck(44, 160, 44), cl(148, 103, 189);
    QVector<Series> s = {line("ε_ρ density change", rho, cr, 1.5), line("ε_W velocity change", vel, cw, 1.5),
                         line("ε_K Kutta", kut, ck, 1.5), line("linear residual |FΨ - Q|∞", lin, cl, 1.2)};
    auto tol = [&](double v, const QColor& c) {
        if (std::isfinite(v) && v > 0.0) s << line({}, {QPointF(0.0, v), QPointF(x_max, v)}, c, 1.0, Qt::DashLine, false);
    };
    tol(tol_rho, cr);
    tol(tol_w, cw);
    tol(tol_k, ck);
    w->logY = true;
    w->setTitle("Residuals");
    w->setLabels("pass", "change / residual");
    w->setSeries(s);
    w->setBadge({status, 0});
    w->commit(false);
}

ContourGeometry contourGeometry(const T::Field& x, const T::Field& y, const T::Field& v, double lo, double hi, int bands,
                                double dy)
{
    ContourGeometry out;
    if (!(hi > lo) || bands < 1) hi = lo + 1.0;
    auto level = [&](int b) { return lo + (hi - lo) * b / bands; };
    auto band_of = [&](double val) { return std::clamp(static_cast<int>(std::floor((val - lo) / (hi - lo) * bands)), 0, bands - 1); };
    std::vector<V3> quad(4), t1, t2;
    for (int i = 0; i + 1 < v.n_m; ++i)
        for (int j = 0; j + 1 < v.n_eta; ++j) {
            const int ci[4] = {i, i + 1, i + 1, i}, cj[4] = {j, j, j + 1, j + 1};
            double vmin = 1e300, vmax = -1e300;
            for (int c = 0; c < 4; ++c) {
                quad[sz(c)] = {x(ci[c], cj[c]), y(ci[c], cj[c]) + dy, v(ci[c], cj[c])};
                vmin = std::min(vmin, quad[sz(c)].v);
                vmax = std::max(vmax, quad[sz(c)].v);
            }
            if (!std::isfinite(vmin) || !std::isfinite(vmax)) continue;
            const int b0 = band_of(vmin), b1 = band_of(vmax);
            for (int b = b0; b <= b1; ++b) {
                const std::vector<V3>* poly = &quad;
                if (b0 != b1) {
                    if (b > 0) {
                        clip(*poly, level(b), true, t1);
                        poly = &t1;
                    }
                    if (b < bands - 1 && poly->size() >= 3) {
                        clip(*poly, level(b + 1), false, t2);
                        poly = &t2;
                    }
                }
                if (poly->size() < 3) continue;
                QPolygonF pg;
                for (const auto& p : *poly) pg << QPointF(p.x, p.y);
                out.fills.push_back({pg, bandColor(b, bands)});
            }
            // iso-lines: marching squares at the band edges inside [vmin, vmax]
            for (int b = std::max(1, b0); b <= std::min(bands - 1, b1 + 1); ++b) {
                const double L = level(b);
                if (L < vmin || L > vmax) continue;
                QPointF p[4];
                bool has[4] = {false, false, false, false};
                int n = 0;
                for (int e = 0; e < 4; ++e) {
                    const V3& a = quad[sz(e)];
                    const V3& c = quad[sz((e + 1) % 4)];
                    if ((a.v >= L) != (c.v >= L)) {
                        const double t = (L - a.v) / (c.v - a.v);
                        p[e] = {a.x + t * (c.x - a.x), a.y + t * (c.y - a.y)};
                        has[e] = true;
                        ++n;
                    }
                }
                if (n == 2) {
                    for (int e = 0; e < 4; ++e)
                        if (has[e]) out.segments << p[e];
                } else if (n == 4) {
                    // saddle: the cell centre value decides which corners the iso-line cuts off
                    const double centre = 0.25 * (quad[0].v + quad[1].v + quad[2].v + quad[3].v);
                    if ((centre >= L) == (quad[0].v >= L)) out.segments << p[0] << p[1] << p[2] << p[3];
                    else out.segments << p[3] << p[0] << p[1] << p[2];
                }
            }
        }
    return out;
}

void drawContours(PlotWidget* w, const T::SolutionData& s, ContourField f, int bands, bool isolines, const QString& title)
{
    const T::Field& v = f == ContourField::Mach ? s.mach : s.pressure;
    const auto& g = s.grid;
    double lo = 1e300, hi = -1e300;
    for (double q : v.v)
        if (std::isfinite(q)) {
            lo = std::min(lo, q);
            hi = std::max(hi, q);
        }
    if (!(hi > lo)) hi = lo + 1e-9 * std::max(1.0, std::abs(lo));
    QVector<Fill> fills;
    QVector<QPointF> iso;
    const double pitch = g.boundaries.periodic_translation.y;
    for (int k = -1; k <= 1; ++k) {
        auto geo = contourGeometry(g.node_m, g.node_y, v, lo, hi, bands, k * pitch);
        fills += geo.fills;
        for (int q = 0; q + 1 < geo.segments.size(); q += 2) iso << geo.segments[q] << geo.segments[q + 1] << QPointF(kNaN, kNaN);
    }
    QVector<Series> series;
    if (isolines) series << line("iso-lines", iso, QColor(0, 0, 0, 130), 0.6, Qt::SolidLine, false);
    for (int k = -1; k <= 2; ++k) {  // copy k spans blade k (suction) to blade k + 1 (pressure)
        const QPolygonF b = bladeOutline(g.boundaries, k);
        fills.push_back({b, QColor(217, 217, 217)});
        QVector<QPointF> o(b.begin(), b.end());
        o << b.front();
        series << line("blade", o, QColor(51, 51, 51), 1.0, Qt::SolidLine, false);
    }
    const QString label = f == ContourField::Mach ? "relative Mach number" : "static pressure [Pa]";
    w->setTitle(title);
    w->setLabels("m / c (axial)", "y / c (tangential)");
    w->proportional = true;
    w->showLegend = false;
    w->setFills(fills);
    w->setMarkers({});
    w->setSeries(series);
    w->setColorBar(true, lo, hi, label, [bands](double t) {
        return bandColor(std::clamp(static_cast<int>(std::floor(t * bands)), 0, bands - 1), bands);
    });
    w->commit(true);
}

void drawLoading(PlotWidget* w, const T::SolutionData& s, LoadingQuantity q, const QString& title)
{
    auto pts = [&](const T::SurfaceData& d) {
        QVector<QPointF> out;
        const auto& val = q == LoadingQuantity::Mach ? d.mach : q == LoadingQuantity::Pressure ? d.pressure : d.dunavant_S;
        for (std::size_t k = 0; k < d.x_c.size() && k < val.size(); ++k) out << QPointF(d.x_c[k], val[k]);
        return out;
    };
    QVector<Series> series = {line("suction side", pts(s.suction), kSuction, 1.8), line("pressure side", pts(s.pressure_side), kPressure, 1.8)};
    for (auto& x : series) x.markers = true, x.marker_size = 2;
    if (q == LoadingQuantity::Mach) series << line("M = 1", {QPointF(0.0, 1.0), QPointF(1.0, 1.0)}, QColor(85, 85, 85), 0.9, Qt::DashLine);
    const QString yl = q == LoadingQuantity::Mach       ? "surface relative Mach number"
                       : q == LoadingQuantity::Pressure ? "surface static pressure p [Pa]"
                                                        : "Dunavant S = (Pt1 - p) / q1";
    w->setTitle(title);
    w->setLabels("x / c", yl);
    w->proportional = false;
    w->setFills({});
    w->setMarkers({});
    w->setSeries(series);
    w->commit(true);
}

}  // namespace Designer::Tl
