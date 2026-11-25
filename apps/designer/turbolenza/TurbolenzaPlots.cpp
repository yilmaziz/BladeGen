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

void drawStreamSheet(PlotWidget* w, const T::StreamSheet& sheet, double m_lo, double m_hi, double m_te)
{
    QVector<QPointF> r, b;
    const int n = 200;
    const double r_le = sheet.radius(0.0), b_lo = sheet.thickness(m_lo);
    double lo = 1.0, hi = 1.0;
    for (int k = 0; k <= n; ++k) {
        const double m = m_lo + (m_hi - m_lo) * k / n;
        r << QPointF(m, sheet.radius(m) / r_le);
        b << QPointF(m, sheet.thickness(m) / b_lo);
        lo = std::min({lo, r.back().y(), b.back().y()});
        hi = std::max({hi, r.back().y(), b.back().y()});
    }
    const double pad = 0.05 * std::max(hi - lo, 0.02);
    QVector<Series> series{line("r / r_LE", r, kSuction, 1.6), line("b / b_in", b, kPressure, 1.6)};
    series << line("LE, TE", {{0.0, lo - pad}, {0.0, hi + pad}}, QColor(120, 120, 120), 1.0, Qt::DashLine)
           << line("", {{m_te, lo - pad}, {m_te, hi + pad}}, QColor(120, 120, 120), 1.0, Qt::DashLine, false);
    w->setTitle(sheet.cylinder() ? "Stream sheet: a cylinder (r, b constant)" : "Stream sheet r(m), b(m)");
    w->setLabels("m / c (from the LE)", "ratio");
    w->proportional = false;
    w->setFills({});
    w->setMarkers({});
    w->setSeries(series);
    w->commit(true);
}

void drawVelocityTriangle(PlotWidget* w, const T::VelocityTriangle& t)
{
    // the vectors from the origin: C and W, and U from the tip of W to the tip of C (C = W + U)
    auto arrow = [](const QString& name, QPointF a, QPointF b, const QColor& c) {
        const QPointF d = b - a;
        const double L = std::hypot(d.x(), d.y());
        QVector<QPointF> pts{a, b};
        if (L > 0.0) {
            const QPointF u = d / L, nrm(-u.y(), u.x());
            const double h = 0.08 * L;
            pts << QPointF(kNaN, kNaN) << b - h * u + 0.5 * h * nrm << b << b - h * u - 0.5 * h * nrm;
        }
        return line(name, pts, c, 1.8);
    };
    const QPointF o(0.0, 0.0), C(t.C_m, t.C_theta), W(t.C_m, t.W_theta);
    QVector<Series> series{arrow(QString("C %1 m/s (alpha %2 deg)").arg(t.C, 0, 'f', 1).arg(t.alpha * 180.0 / std::numbers::pi, 0, 'f', 2), o, C,
                                 QColor(31, 119, 180)),
                           arrow(QString("W %1 m/s (beta %2 deg)").arg(t.W, 0, 'f', 1).arg(t.beta * 180.0 / std::numbers::pi, 0, 'f', 2), o, W,
                                 QColor(214, 39, 40)),
                           arrow(QString("U = omega r %1 m/s").arg(t.U, 0, 'f', 1), W, C, QColor(0, 150, 0))};
    w->setTitle(QString("Inlet velocity triangle: rpm %1, r_in %2 m, M1 %3 abs / %4 rel")
                    .arg(t.rpm, 0, 'g', 7)
                    .arg(t.radius, 0, 'g', 5)
                    .arg(t.mach_absolute, 0, 'f', 3)
                    .arg(t.mach_relative, 0, 'f', 3));
    w->setLabels("meridional [m/s]", "tangential, + toward +theta [m/s]");
    w->proportional = true;
    w->setFills({});
    w->setMarkers({});
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

namespace {
// Dunavant (1955), measured blade loading of the NACA 65-(12A10)10 cascade (x/c, Mach number, S), as supplied
// by the user; kept equal to tests/turbolenza/validation/dunavant1955/loading.csv
struct DunavantPoint
{
    double x_c, mach, S;
};
const DunavantPoint kDunavant1955Suction[] = {
    {0.0211, 0.6318, 1.0586},
    {0.0518, 0.75, 1.3987},
    {0.1025, 0.8204, 1.6047},
    {0.1986, 0.8863, 1.7962},
    {0.302, 0.8697, 1.7481},
    {0.4019, 0.8531, 1.7},
    {0.5035, 0.8185, 1.5991},
    {0.6016, 0.7513, 1.4025},
    {0.707, 0.6801, 1.1962},
    {0.7996, 0.595, 0.9564},
    {0.9049, 0.5561, 0.8508}};
const DunavantPoint kDunavant1955Pressure[] = {
    {0.0211, 0.6436, 1.0921},
    {0.0593, 0.5762, 0.9051},
    {0.1029, 0.5513, 0.838},
    {0.1955, 0.5348, 0.7947},
    {0.3025, 0.5256, 0.7705},
    {0.4024, 0.5087, 0.7272},
    {0.5022, 0.48, 0.6551},
    {0.6002, 0.4522, 0.5878},
    {0.7019, 0.4357, 0.5493},
    {0.7999, 0.4189, 0.5108}};
// Aungier (2003) Figs. 5-6 and 5-7, the potential-flow curves (solid lines), digitized
// (tools/turbolenza_reference/digitize_aungier_fig5_6_5_7.py; tests/turbolenza/validation/aungier2003_fig5_6_5_7)
const QPointF kFig56Suction[] = {{0.0000, 0.7587}, {0.0435, 0.8692}, {0.0870, 0.9036}, {0.1304, 0.9145}, {0.1739, 0.9188}, {0.2174, 0.9176}, {0.2609, 0.9089}, {0.3043, 0.8983}, {0.3478, 0.8875}, {0.3913, 0.8744}, {0.4348, 0.8503}, {0.4783, 0.8191}, {0.5217, 0.7839}, {0.5652, 0.7502}, {0.6087, 0.7166}, {0.6522, 0.6872}, {0.6957, 0.6553}, {0.7391, 0.6268}, {0.7826, 0.5991}, {0.8261, 0.5694}, {0.8696, 0.5370}, {0.9130, 0.4981}, {0.9565, 0.4529}, {1.0000, 0.3979}};
const QPointF kFig56Pressure[] = {{0.0000, 0.6580}, {0.0435, 0.5606}, {0.0870, 0.5341}, {0.1304, 0.5301}, {0.1739, 0.5211}, {0.2174, 0.5162}, {0.2609, 0.5088}, {0.3043, 0.5000}, {0.3478, 0.4930}, {0.3913, 0.4843}, {0.4348, 0.4737}, {0.4783, 0.4613}, {0.5217, 0.4492}, {0.5652, 0.4368}, {0.6087, 0.4244}, {0.6522, 0.4122}, {0.6957, 0.4015}, {0.7391, 0.3926}, {0.7826, 0.3820}, {0.8261, 0.3751}, {0.8696, 0.3678}, {0.9130, 0.3631}, {0.9565, 0.3676}, {1.0000, 0.3905}};
const QPointF kFig56Separation{0.778, 0.587};  // the predicted separation point (filled triangle)
const QPointF kFig57Suction[] = {{0.0000, 0.9001}, {0.0200, 0.9480}, {0.0400, 0.9959}, {0.0600, 1.0438}, {0.0800, 1.0725}, {0.1000, 1.0873}, {0.1200, 1.1039}, {0.1400, 1.1230}, {0.1600, 1.1283}, {0.1800, 1.1283}, {0.2000, 1.1300}, {0.2200, 1.1283}, {0.2400, 1.1204}, {0.2600, 1.1109}, {0.2800, 1.1013}, {0.3000, 1.0908}, {0.3200, 1.0821}, {0.3400, 1.0778}, {0.3600, 1.0586}, {0.3800, 1.0394}, {0.4000, 1.0203}, {0.4200, 1.0185}, {0.4400, 0.9994}, {0.4600, 0.9767}, {0.4800, 0.9523}, {0.5000, 0.9297}, {0.5200, 0.9071}, {0.5400, 0.8844}, {0.5600, 0.8635}, {0.5800, 0.8435}, {0.6000, 0.8243}, {0.6200, 0.8060}, {0.6400, 0.7869}, {0.6600, 0.7677}, {0.6800, 0.7503}, {0.7000, 0.7337}, {0.7200, 0.7172}, {0.7400, 0.6998}, {0.7600, 0.6841}, {0.7800, 0.6675}, {0.8000, 0.6519}, {0.8200, 0.6362}, {0.8400, 0.6188}, {0.8600, 0.6014}, {0.8800, 0.5822}, {0.9000, 0.5622}, {0.9200, 0.5413}, {0.9400, 0.5151}, {0.9600, 0.4908}, {0.9800, 0.4611}, {1.0000, 0.4368}};
const QPointF kFig57Pressure[] = {{0.0000, 0.8008}, {0.0200, 0.7355}, {0.0400, 0.6702}, {0.0600, 0.6527}, {0.0800, 0.6379}, {0.1000, 0.6310}, {0.1200, 0.6249}, {0.1400, 0.6188}, {0.1600, 0.6118}, {0.1800, 0.6048}, {0.2000, 0.5996}, {0.2200, 0.5944}, {0.2400, 0.5883}, {0.2600, 0.5822}, {0.2800, 0.5770}, {0.3000, 0.5709}, {0.3200, 0.5657}, {0.3400, 0.5596}, {0.3600, 0.5543}, {0.3800, 0.5482}, {0.4000, 0.5430}, {0.4200, 0.5369}, {0.4400, 0.5299}, {0.4600, 0.5238}, {0.4800, 0.5160}, {0.5000, 0.5090}, {0.5200, 0.5003}, {0.5400, 0.4934}, {0.5600, 0.4873}, {0.5800, 0.4803}, {0.6000, 0.4742}, {0.6200, 0.4672}, {0.6400, 0.4611}, {0.6600, 0.4559}, {0.6800, 0.4498}, {0.7000, 0.4446}, {0.7200, 0.4394}, {0.7400, 0.4341}, {0.7600, 0.4307}, {0.7800, 0.4246}, {0.8000, 0.4211}, {0.8200, 0.4159}, {0.8400, 0.4141}, {0.8600, 0.4106}, {0.8800, 0.4071}, {0.9000, 0.4037}, {0.9200, 0.4037}, {0.9400, 0.4071}, {0.9600, 0.4106}, {0.9800, 0.4237}, {1.0000, 0.4368}};


}  // namespace

void drawLoading(PlotWidget* w, const T::SolutionData& s, LoadingQuantity q, const QString& title, const LoadingOverlays& o)
{
    const T::ViscousData* inviscid = o.inviscid;
    const bool dunavant1955 = o.dunavant1955;
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
    if (dunavant1955 && q != LoadingQuantity::Pressure) {  // measured points only (no p in Pa)
        auto measured = [&](const char* name, const auto& data, const QColor& c, char shape) {
            QVector<QPointF> out;
            for (const auto& d : data) out << QPointF(d.x_c, q == LoadingQuantity::Mach ? d.mach : d.S);
            Series m = line(name, out, c, 1.4, Qt::NoPen);
            m.markers = true;
            m.marker_size = 4;
            m.marker_shape = shape;
            return m;
        };
        series << measured("suction side, Dunavant (1955) measured", kDunavant1955Suction, Qt::black, 'o')
               << measured("pressure side, Dunavant (1955) measured", kDunavant1955Pressure, Qt::black, '^');
    }
    QVector<PointMarker> marks;
    if (q == LoadingQuantity::Mach) {
        auto curve = [](const auto& data) {
            QVector<QPointF> out;
            for (const auto& pt : data) out << pt;
            return out;
        };
        if (o.aungier_fig5_6) {
            const QString n = "Aungier (2003) Fig. 5-6, potential flow (M_in 0.65)";
            series << line(n, curve(kFig56Suction), Qt::black, 1.2) << line(n, curve(kFig56Pressure), Qt::black, 1.2);
            marks << PointMarker{kFig56Separation, Qt::black, '^', "Fig. 5-6 separation", 6.0, 2.0};
        }
        if (o.aungier_fig5_7) {
            const QString n = "Aungier (2003) Fig. 5-7, potential flow (M_in 0.757)";
            series << line(n, curve(kFig57Suction), Qt::black, 1.2, Qt::DashDotLine)
                   << line(n, curve(kFig57Pressure), Qt::black, 1.2, Qt::DashDotLine);
        }
    }
    if (o.layers) {  // the boundary layer's separation points on the reported loading
        auto at = [&](const T::SurfaceData& d, double xc) {
            const auto& val = q == LoadingQuantity::Mach ? d.mach : q == LoadingQuantity::Pressure ? d.pressure : d.dunavant_S;
            for (std::size_t k = 1; k < d.x_c.size() && k < val.size(); ++k)
                if ((d.x_c[k - 1] - xc) * (d.x_c[k] - xc) <= 0.0 && d.x_c[k] != d.x_c[k - 1])
                    return val[k - 1] + (val[k] - val[k - 1]) * (xc - d.x_c[k - 1]) / (d.x_c[k] - d.x_c[k - 1]);
            return std::numeric_limits<double>::quiet_NaN();
        };
        const T::BoundaryLayerSurfaceData* layer[2] = {&o.layers->suction, &o.layers->pressure_side};
        const T::SurfaceData* loading[2] = {&s.suction, &s.pressure_side};
        for (int k = 0; k < 2; ++k) {
            const QColor c = k == 0 ? kSuction : kPressure;
            const double turbulent = layer[k]->separation_x_c, laminar = layer[k]->laminar_separation_x_c;
            for (const double xc : {turbulent, laminar}) {
                if (!std::isfinite(xc)) continue;
                const double y = at(*loading[k], xc);
                if (std::isfinite(y))
                    marks << PointMarker{QPointF(xc, y), c, 'x', xc == turbulent ? "separation (BL)" : "laminar separation (BL)", 7.0, 2.0};
            }
        }
    }
    if (q == LoadingQuantity::Mach) series << line("M = 1", {QPointF(0.0, 1.0), QPointF(1.0, 1.0)}, QColor(85, 85, 85), 0.9, Qt::DashLine);
    const QString yl = q == LoadingQuantity::Mach       ? "surface relative Mach number"
                       : q == LoadingQuantity::Pressure ? "surface static pressure p [Pa]"
                                                        : "Dunavant S = (Pt1 - p) / q1";
    w->setTitle(title);
    w->setLabels("x / c", yl);
    w->proportional = false;
    w->setFills({});
    w->setMarkers(marks);
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

void drawBoundaryLayer(PlotWidget* w, const T::ViscousData& v, double chord, BlQuantity q, bool overlay, const QString& title,
                       BlAbscissa abscissa)
{
    const bool on_f = abscissa == BlAbscissa::StationFraction;
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
        const std::vector<double>& xs = on_f ? d->f : d->x_c;
        const std::size_t n = std::min(xs.size(), d->H.size());
        QVector<QPointF> p;
        for (std::size_t k = 0; k < n; ++k)
            if (std::isfinite(value(*d, k))) p << QPointF(xs[k], value(*d, k));
        auto s = line(name, p, c, 1.8);
        s.markers = true;
        s.marker_size = 2;
        series << s;
        auto mark = [&](double x, char shape, const QString& label) {
            if (!std::isfinite(x)) return;
            for (std::size_t k = 0; k < n; ++k)
                if (xs[k] == x && std::isfinite(value(*d, k))) markers << marker(QPointF(x, value(*d, k)), shape, c, label);
        };
        mark(on_f ? d->transition_f : d->transition_x_c, 'o', "transition");
        mark(on_f ? d->separation_f : d->separation_x_c, 'x', "separation");
        mark(on_f ? d->laminar_separation_f : d->laminar_separation_x_c, 's', "laminar separation");
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
    w->setLabels(on_f ? "blade-station fraction f = (m − m_first) / (m_last − m_first)" : "x / c", yl);
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
