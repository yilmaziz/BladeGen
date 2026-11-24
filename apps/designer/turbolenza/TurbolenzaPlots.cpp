#include "TurbolenzaPlots.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <tuple>

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
    QVector<PointMarker> markers;
    if (g.touching()) {
        // the touching stations (EXPERIMENT): both wall nodes on the nose tip / tail end
        for (int i : {g.first_blade_i, g.last_blade_i})
            for (int j : {0, ne - 1}) {
                PointMarker m;
                m.p = QPointF(g.node_m(i, j), g.node_y(i, j));
                m.shape = 'o';
                m.color = QColor(230, 120, 0);
                m.size = 5;
                if (j == 0) m.label = i == g.first_blade_i ? "touching LE station" : "touching TE station";
                markers << m;
            }
    }
    w->setTitle(QString("Grid %1 x %2 (stations x eta), blade stations %3 .. %4%5")
                    .arg(nm)
                    .arg(ne)
                    .arg(g.first_blade_i)
                    .arg(g.last_blade_i)
                    .arg(g.touching() ? ", touching LE / TE stations (experiment)" : ""));
    w->setLabels("m / c (axial)", "y / c (tangential)");
    w->proportional = true;
    w->setFills(fills);
    w->setMarkers(markers);
    w->setSeries(series);
    w->commit(true);
}

void drawResiduals(PlotWidget* w, const std::vector<T::PassRecord>& records, double tol_rho, double tol_w, double tol_k,
                   const QString& status)
{
    QVector<QPointF> rho, vel, kut, lin;
    double x_max = 1.0, position = 0.0;
    for (const auto& r : records) {
        // the pass's position in the run: its index for an inviscid run, the passes of every boundary-layer
        // iteration one after the other for a viscous one (their indices restart at 0)
        const double x = position;
        position += 1.0;
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

void drawLoading(PlotWidget* w, const T::SolutionData& s, LoadingQuantity q, const QString& title, const T::ViscousData* inviscid)
{
    auto pts = [&](const T::SurfaceData& d) {
        QVector<QPointF> out;
        const auto& val = q == LoadingQuantity::Mach ? d.mach : q == LoadingQuantity::Pressure ? d.pressure : d.dunavant_S;
        for (std::size_t k = 0; k < d.x_c.size() && k < val.size(); ++k) out << QPointF(d.x_c[k], val[k]);
        return out;
    };
    const bool both = inviscid && !inviscid->inviscid_suction.x_c.empty();
    QVector<Series> series = {line(both ? "suction side (coupled)" : "suction side", pts(s.suction), kSuction, 1.8),
                              line(both ? "pressure side (coupled)" : "pressure side", pts(s.pressure_side), kPressure, 1.8)};
    for (auto& x : series) x.markers = true, x.marker_size = 2;
    if (both)
        series << line("suction side (inviscid, k = 0)", pts(inviscid->inviscid_suction), kSuction, 1.1, Qt::DashLine)
               << line("pressure side (inviscid, k = 0)", pts(inviscid->inviscid_pressure_side), kPressure, 1.1, Qt::DashLine);
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

namespace {

// Aungier (2003) Fig. 5-15, digitized (tests/turbolenza/validation/aungier2003_fig5_15/H_*.csv)
const QPointF kFigSuction[] = {
    {0.0000, 2.9547}, {0.0426, 2.8819}, {0.0863, 2.9972}, {0.1297, 1.8969}, {0.1662, 1.9048},
    {0.2579, 1.9090}, {0.3061, 1.9006}, {0.3943, 1.9006}, {0.4319, 1.9133}, {0.4778, 1.9578},
    {0.5224, 2.0320}, {0.5660, 2.1167}, {0.6059, 2.2100}, {0.6506, 2.3604}, {0.6965, 2.6042},
    {0.7376, 2.6656}, {0.7459, 2.6635}, {0.7517, 2.6571}, {1.0000, 2.5123},
};
const QPointF kFigPressure[] = {
    {0.0000, 2.8802}, {0.0429, 3.2135}, {0.1296, 2.7996}, {0.1732, 1.6553}, {0.1873, 1.6590},
    {0.2191, 1.6674}, {0.3884, 1.6717}, {0.4154, 1.6844}, {0.4319, 1.6844}, {0.4789, 1.7119},
    {0.7470, 1.9196}, {0.8329, 1.9726}, {0.8705, 1.9874}, {0.9093, 1.9663}, {0.9140, 1.9578},
    {0.9587, 1.7819}, {0.9975, 1.5827}, {1.0000, 1.5766},
};
const QPointF kFigTransition[] = {{0.1297, 1.8969}, {0.1732, 1.6553}};  // suction, pressure
const QPointF kFigSeparation(0.7459, 2.6635);                           // suction: the dashed part starts

}  // namespace

void drawBoundaryLayer(PlotWidget* w, const T::ViscousData& v, double chord, BlQuantity q, bool overlay, const QString& title)
{
    auto value = [&](const T::BoundaryLayerSurfaceData& d, std::size_t k) {
        switch (q) {
        case BlQuantity::H: return d.H[k];
        case BlQuantity::ThetaOverC: return d.theta[k] / chord;
        case BlQuantity::DeltaStarOverC: return d.delta_star[k] / chord;
        case BlQuantity::Cf: return d.c_f[k];
        case BlQuantity::ReTheta: return d.Re_theta[k];
        }
        return kNaN;
    };
    QVector<Series> series;
    QVector<PointMarker> markers;
    auto marker = [](QPointF p, char shape, const QColor& c, const QString& label) {
        PointMarker m;
        m.p = p;
        m.shape = shape;
        m.color = c;
        m.size = 6.0;
        m.label = label;
        return m;
    };
    for (const auto& [d, c, name] : {std::tuple{&v.suction, kSuction, QString("suction surface")},
                                     std::tuple{&v.pressure_side, kPressure, QString("pressure surface")}}) {
        QVector<QPointF> p;
        for (std::size_t k = 0; k < d->x_c.size(); ++k)
            if (std::isfinite(value(*d, k))) p << QPointF(d->x_c[k], value(*d, k));
        auto s = line(name, p, c, 1.8);
        s.markers = true;
        s.marker_size = 2;
        series << s;
        auto mark = [&](double x, char shape, const QString& label) {
            if (!std::isfinite(x)) return;
            for (std::size_t k = 0; k < d->x_c.size(); ++k)
                if (d->x_c[k] == x && std::isfinite(value(*d, k))) markers << marker(QPointF(x, value(*d, k)), shape, c, label);
        };
        mark(d->transition_x_c, 'o', "transition");
        mark(d->separation_x_c, 'x', "separation");
        mark(d->laminar_separation_x_c, 's', "laminar separation");
    }
    if (overlay && q == BlQuantity::H) {
        QVector<QPointF> ss(std::begin(kFigSuction), std::end(kFigSuction)), ps(std::begin(kFigPressure), std::end(kFigPressure));
        auto a = line("Aungier 2003 Fig. 5-15, suction", ss, kSuction.darker(150), 0.9, Qt::DotLine);
        auto b = line("Aungier 2003 Fig. 5-15, pressure", ps, kPressure.darker(150), 0.9, Qt::DotLine);
        a.markers = b.markers = true;
        a.marker_size = b.marker_size = 4;
        series << a << b;
        for (const QPointF& t : kFigTransition) markers << marker(t, '^', QColor(60, 60, 60), "Fig. 5-15 transition");
        markers << marker(kFigSeparation, 'x', QColor(60, 60, 60), "Fig. 5-15 separation");
    }
    const QString yl = q == BlQuantity::H              ? "shape factor H = δ*/θ"
                       : q == BlQuantity::ThetaOverC     ? "θ / c"
                       : q == BlQuantity::DeltaStarOverC ? "δ* / c"
                       : q == BlQuantity::Cf             ? "skin friction c_f"
                                                         : "Re_θ";
    w->setTitle(title);
    w->setLabels("x / c", yl);
    w->logY = false;
    w->proportional = false;
    w->setFills({});
    w->setMarkers(markers);
    w->setSeries(series);
    if (v.suction.x_c.empty())
        w->setBadge(Badge{"no boundary layer in this run (it failed before one was computed)", 1});
    else if (!v.layers_of_last_solve)
        w->setBadge(Badge{QString("the run stopped in BL iteration %1: this is the boundary layer of iteration %2, the last one "
                                  "computed")
                              .arg(v.layers_iteration + 1)
                              .arg(v.layers_iteration),
                          1});
    else
        w->setBadge(Badge{});
    w->commit(true);
}

void drawBlChange(PlotWidget* w, const T::ViscousData& v, double tolerance)
{
    QVector<QPointF> p;
    for (const auto& h : v.history) p << QPointF(h.k, h.max_change_m);
    auto s = line("max |m_BL - m_used| / max |m_BL|", p, QColor(31, 119, 180), 1.6);
    s.markers = true;
    QVector<Series> series{s};
    if (std::isfinite(tolerance) && tolerance > 0.0 && !p.isEmpty())
        series << line("tolerance", {QPointF(0.0, tolerance), QPointF(std::max(1.0, p.back().x()), tolerance)}, QColor(85, 85, 85),
                       1.0, Qt::DashLine);
    w->logY = true;
    w->setTitle(QString("Coupling history: change of the mass deficit (%1)").arg(v.converged ? "converged" : "NOT converged"));
    w->setLabels("BL iteration k", "max relative change of m");
    w->proportional = false;
    w->setFills({});
    w->setMarkers({});
    w->setSeries(series);
    w->setBadge({});
    w->commit(true);
}

void drawBlBetaOut(PlotWidget* w, const T::ViscousData& v)
{
    QVector<QPointF> p;
    for (const auto& h : v.history) p << QPointF(h.k, h.beta_out * 180.0 / std::numbers::pi);
    auto s = line("β_out of the inviscid solve of iteration k", p, QColor(214, 39, 40), 1.6);
    s.markers = true;
    w->logY = false;
    w->setTitle("Coupling history: β_out");
    w->setLabels("BL iteration k", "β_out [deg]");
    w->proportional = false;
    w->setFills({});
    w->setMarkers({});
    w->setSeries({s});
    w->setBadge({});
    w->commit(true);
}

namespace {

/// One fixed colour per standard level (the order of kGridLevels); other names grey.
QColor levelColor(const std::string& level)
{
    static const QColor colors[] = {QColor(148, 103, 189), QColor(31, 119, 180), QColor(44, 160, 44), QColor(255, 127, 14), QColor(214, 39, 40)};
    for (std::size_t k = 0; k < T::kGridLevels.size(); ++k)
        if (level == T::kGridLevels[k]) return colors[k];
    return QColor(110, 110, 110);
}

}  // namespace

void drawStudyBetaOut(PlotWidget* w, const std::vector<T::MeshStudyEntry>& levels, bool x_nodes)
{
    QVector<QPointF> ok, bad;
    QVector<PointMarker> markers;
    QStringList not_converged;
    for (const auto& e : levels) {
        if (!e.converged) not_converged << QString::fromStdString(e.level);
        const double x = x_nodes ? static_cast<double>(e.total_nodes) : e.dm_over_chord;
        if (!std::isfinite(x) || !std::isfinite(e.beta_out_deg) || (x_nodes && e.total_nodes <= 0)) continue;
        const QPointF p(x, e.beta_out_deg);
        if (e.converged) ok << p;
        else bad << p << QPointF(kNaN, kNaN);  // markers only, no line between them
        PointMarker m;
        m.p = p;
        m.shape = e.converged ? 'o' : 'x';
        m.color = e.converged ? kSuction : kPressure;
        m.size = e.converged ? 4.0 : 6.0;
        m.label = QString::fromStdString(e.level);
        markers << m;
    }
    std::sort(ok.begin(), ok.end(), [](const QPointF& a, const QPointF& b) { return a.x() < b.x(); });
    QVector<Series> series;
    if (!ok.isEmpty()) series << line("converged levels", ok, kSuction, 1.2);
    if (!bad.isEmpty()) {
        auto s = line("not converged (last pass)", bad, kPressure, 1.0);
        s.markers = true;
        s.marker_size = 3;
        series << s;
    }
    QString title = "Mesh study: β_out against the grid";
    if (!not_converged.isEmpty()) title += "   [NOT CONVERGED: " + not_converged.join(", ") + "]";
    w->setTitle(title);
    w->setLabels(x_nodes ? "total nodes" : "dm / c", "β_out [deg]");
    w->proportional = false;
    w->setFills({});
    w->setMarkers(markers);
    w->setSeries(series);
    w->setBadge(series.isEmpty() ? Badge{"no level with a solution", 1} : Badge{});
    w->commit(true);
}

void drawStudyLoading(PlotWidget* w, const std::vector<T::MeshStudyEntry>& levels, LoadingQuantity q)
{
    auto pts = [&](const T::SurfaceData& d) {
        QVector<QPointF> out;
        const auto& val = q == LoadingQuantity::Mach ? d.mach : q == LoadingQuantity::Pressure ? d.pressure : d.dunavant_S;
        for (std::size_t k = 0; k < d.x_c.size() && k < val.size(); ++k) out << QPointF(d.x_c[k], val[k]);
        return out;
    };
    QVector<Series> series;
    for (const auto& e : levels) {
        if (!e.result || !e.result->solution) continue;
        const auto& s = *e.result->solution;
        if (s.suction.x_c.empty()) continue;
        const QColor c = levelColor(e.level);
        const QString name = QString::fromStdString(e.level) + (e.converged ? "" : " (not converged)");
        series << line(name, pts(s.suction), c, 1.6) << line(name + " PS", pts(s.pressure_side), c, 1.6, Qt::DashLine, false);
    }
    if (q == LoadingQuantity::Mach && !series.isEmpty())
        series << line("M = 1", {QPointF(0.0, 1.0), QPointF(1.0, 1.0)}, QColor(85, 85, 85), 0.9, Qt::DotLine);
    const QString yl = q == LoadingQuantity::Mach       ? "surface relative Mach number"
                       : q == LoadingQuantity::Pressure ? "surface static pressure p [Pa]"
                                                        : "Dunavant S = (Pt1 - p) / q1";
    w->setTitle("Mesh study: blade loading (solid: suction side, dashed: pressure side)");
    w->setLabels("x / c", yl);
    w->proportional = false;
    w->setFills({});
    w->setMarkers({});
    w->setSeries(series);
    w->setBadge(series.isEmpty() ? Badge{"no level with blade surfaces", 1} : Badge{});
    w->commit(true);
}

}  // namespace Designer::Tl
