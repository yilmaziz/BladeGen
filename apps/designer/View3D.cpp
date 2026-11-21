#include "View3D.hpp"

#include <QContextMenuEvent>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

namespace Designer {

namespace {
constexpr double kDeg = 3.14159265358979323846 / 180.0;

struct Item
{
    double depth;
    int kind;  // 0 quad, 1 segment, 2 point
    QPointF q[4];
    QColor color;
    Qt::PenStyle style;
    double width;
};
}  // namespace

View3D::View3D(QWidget* parent) : QWidget(parent)
{
    setMinimumSize(260, 200);
    setMouseTracking(false);
}

void View3D::setScene(std::vector<Grid3D> grids, std::vector<Line3D> lines, bool refit)
{
    grids_ = std::move(grids);
    lines_ = std::move(lines);
    if (refit || !fitted_) zoomFit();
    else update();
}

void View3D::zoomFit()
{
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    bool any = false;
    for (const auto& g : grids_)
        for (const auto& p : g.p)
            for (int k = 0; k < 3; ++k) {
                lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
                hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
                any = true;
            }
    if (!any) return;
    for (int k = 0; k < 3; ++k) center_[static_cast<std::size_t>(k)] = 0.5 * (lo[k] + hi[k]);
    extent_ = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-9});
    zoom_ = 1.0;
    pan_ = {0, 0};
    fitted_ = true;
    update();
}

void View3D::setView(double yaw_deg, double pitch_deg)
{
    yaw_ = yaw_deg;
    pitch_ = pitch_deg;
    update();
}

std::array<double, 3> View3D::view(const std::array<double, 3>& p) const
{
    const double qx = p[0] - center_[0], qy = p[1] - center_[1], qz = p[2] - center_[2];
    // yaw about the world z axis (radial direction at theta = 0, drawn upward)
    const double cy = std::cos(yaw_ * kDeg), sy = std::sin(yaw_ * kDeg);
    const double ax = cy * qx - sy * qy, ay = sy * qx + cy * qy, az = qz;
    // view frame: x axial to the right, z up, -y toward the viewer; then pitch about the view x
    const double X = ax, Y = az, Z = -ay;
    const double cp = std::cos(pitch_ * kDeg), sp = std::sin(pitch_ * kDeg);
    return {X, cp * Y - sp * Z, sp * Y + cp * Z};
}

QPointF View3D::screen(const std::array<double, 3>& v) const
{
    const double s = 0.85 * zoom_ * std::min(width(), height()) / extent_;
    return {0.5 * width() + pan_.x() + s * v[0], 0.5 * height() + pan_.y() - s * v[1]};
}

void View3D::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    p.setRenderHint(QPainter::Antialiasing, true);
    std::vector<Item> items;
    const std::array<double, 3> light = {0.35, 0.55, 0.76};
    for (const auto& g : grids_) {
        if (g.ni < 1 || g.nj < 1) continue;
        std::vector<std::array<double, 3>> v(g.p.size());
        std::vector<QPointF> s(g.p.size());
        for (std::size_t k = 0; k < g.p.size(); ++k) {
            v[k] = view(g.p[k]);
            s[k] = screen(v[k]);
        }
        auto idx = [&](int i, int j) { return static_cast<std::size_t>(i + g.ni * j); };
        const bool quads = shaded && g.surface;
        if (quads)
            for (int j = 0; j + 1 < g.nj; ++j)
                for (int i = 0; i + 1 < g.ni; ++i) {
                    const auto &a = v[idx(i, j)], &b = v[idx(i + 1, j)], &c = v[idx(i + 1, j + 1)], &d = v[idx(i, j + 1)];
                    const double ux = c[0] - a[0], uy = c[1] - a[1], uz = c[2] - a[2];
                    const double wx = d[0] - b[0], wy = d[1] - b[1], wz = d[2] - b[2];
                    double nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
                    const double nn = std::sqrt(nx * nx + ny * ny + nz * nz);
                    if (nn > 0) {
                        nx /= nn;
                        ny /= nn;
                        nz /= nn;
                    }
                    const double lam = std::abs(nx * light[0] + ny * light[1] + nz * light[2]);
                    const double f = 0.35 + 0.65 * lam;
                    const QColor col(std::clamp(static_cast<int>(g.color.red() * f + 40 * lam), 0, 255),
                                     std::clamp(static_cast<int>(g.color.green() * f + 40 * lam), 0, 255),
                                     std::clamp(static_cast<int>(g.color.blue() * f + 40 * lam), 0, 255));
                    Item it{0.25 * (a[2] + b[2] + c[2] + d[2]), 0, {s[idx(i, j)], s[idx(i + 1, j)], s[idx(i + 1, j + 1)], s[idx(i, j + 1)]},
                            col, Qt::SolidLine, 1.0};
                    items.push_back(it);
                }
        // grid lines (in shaded mode: darker edges of the surface; walls always as lines)
        const QColor lc = quads ? g.color.darker(170) : g.color;
        for (int j = 0; j < g.nj; ++j)
            for (int i = 0; i + 1 < g.ni; ++i) {
                const auto &a = v[idx(i, j)], &b = v[idx(i + 1, j)];
                items.push_back({0.5 * (a[2] + b[2]) + (quads ? 1e-6 * extent_ : 0.0), 1, {s[idx(i, j)], s[idx(i + 1, j)]}, lc, Qt::SolidLine, 1.0});
            }
        for (int i = 0; i < g.ni; ++i)
            for (int j = 0; j + 1 < g.nj; ++j) {
                const auto &a = v[idx(i, j)], &b = v[idx(i, j + 1)];
                items.push_back({0.5 * (a[2] + b[2]) + (quads ? 1e-6 * extent_ : 0.0), 1, {s[idx(i, j)], s[idx(i, j + 1)]}, lc, Qt::SolidLine, 1.0});
            }
        if (showPoints || g.points)
            for (std::size_t k = 0; k < v.size(); ++k) items.push_back({v[k][2], 2, {s[k]}, g.color.darker(150), Qt::SolidLine, 1.0});
    }
    for (const auto& l : lines_)
        for (std::size_t k = 0; k + 1 < l.p.size(); ++k) {
            const auto a = view(l.p[k]), b = view(l.p[k + 1]);
            items.push_back({shaded ? 0.5 * (a[2] + b[2]) : 1e300, 1, {screen(a), screen(b)}, l.color, l.style, l.width});
        }
    if (shaded) std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.depth < b.depth; });
    for (const auto& it : items) {
        if (it.kind == 0) {
            p.setPen(QPen(it.color, 0));
            p.setBrush(it.color);
            p.drawPolygon(it.q, 4);
        } else if (it.kind == 1) {
            QPen pen(it.color, it.width, it.style);
            pen.setCosmetic(true);
            p.setPen(pen);
            p.drawLine(it.q[0], it.q[1]);
        } else {
            p.setPen(Qt::NoPen);
            p.setBrush(it.color);
            p.drawEllipse(it.q[0], 1.4, 1.4);
        }
    }
    // axis triad (bottom left): x axial red, y green, z blue
    const QPointF o(46, height() - 46);
    const std::array<std::array<double, 3>, 3> axes = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const QColor ac[3] = {QColor(200, 0, 0), QColor(0, 150, 0), QColor(0, 0, 220)};
    const char* an[3] = {"x", "y", "z"};
    for (int k = 0; k < 3; ++k) {
        const auto c0 = view(center_);
        auto e = axes[static_cast<std::size_t>(k)];
        for (int q = 0; q < 3; ++q) e[static_cast<std::size_t>(q)] = center_[static_cast<std::size_t>(q)] + e[static_cast<std::size_t>(q)];
        const auto c1 = view(e);
        const QPointF d(30 * (c1[0] - c0[0]), -30 * (c1[1] - c0[1]));
        QPen pen(ac[k], 1.5);
        pen.setCosmetic(true);
        p.setPen(pen);
        p.drawLine(o, o + d);
        p.drawText(o + d * 1.15 + QPointF(-3, 4), an[k]);
    }
    p.setPen(Qt::darkGray);
    p.drawText(8, 16, QString("yaw %1°  pitch %2°   drag: rotate · shift/middle drag: pan · wheel: zoom").arg(yaw_, 0, 'f', 0).arg(pitch_, 0, 'f', 0));
}

void View3D::mousePressEvent(QMouseEvent* e)
{
    last_ = e->position();
    if (e->button() == Qt::MiddleButton || (e->button() == Qt::LeftButton && (e->modifiers() & Qt::ShiftModifier))) drag_ = 2;
    else if (e->button() == Qt::LeftButton) drag_ = 1;
}

void View3D::mouseMoveEvent(QMouseEvent* e)
{
    const QPointF d = e->position() - last_;
    last_ = e->position();
    if (drag_ == 1) {
        yaw_ += 0.5 * d.x();
        pitch_ = std::clamp(pitch_ + 0.5 * d.y(), -89.0, 89.0);
        update();
    } else if (drag_ == 2) {
        pan_ += d;
        update();
    }
}

void View3D::mouseReleaseEvent(QMouseEvent*) { drag_ = 0; }
void View3D::mouseDoubleClickEvent(QMouseEvent*) { zoomFit(); }

void View3D::wheelEvent(QWheelEvent* e)
{
    zoom_ *= std::pow(1.15, e->angleDelta().y() / 120.0);
    update();
}

void View3D::contextMenuEvent(QContextMenuEvent* e) { emit contextMenuWanted(e->globalPos()); }

QImage View3D::image()
{
    QImage img(size() * devicePixelRatioF(), QImage::Format_ARGB32);
    img.setDevicePixelRatio(devicePixelRatioF());
    img.fill(Qt::white);
    render(&img);
    return img;
}

}  // namespace Designer
