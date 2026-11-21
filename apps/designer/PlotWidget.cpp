#include "PlotWidget.hpp"

#include <QContextMenuEvent>
#include <QImage>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

namespace Designer {

namespace {

double symlog(double v) { return (v < 0 ? -1.0 : 1.0) * std::log10(1.0 + std::abs(v)); }
double symexp(double t) { return (t < 0 ? -1.0 : 1.0) * (std::pow(10.0, std::abs(t)) - 1.0); }

double niceStep(double range, int n)
{
    if (!(range > 0)) return 1.0;
    const double raw = range / std::max(n, 1);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / mag;
    const double nf = f < 1.5 ? 1.0 : f < 3.0 ? 2.0 : f < 7.0 ? 5.0 : 10.0;
    return nf * mag;
}

QString fmt(double v, double step)
{
    if (std::abs(v) < 1e-12 * std::max(1.0, std::abs(step))) v = 0.0;
    const int dec = std::clamp(static_cast<int>(-std::floor(std::log10(std::abs(step)) + 1e-9)), 0, 8);
    if (std::abs(v) >= 1e5 || (std::abs(v) < 1e-4 && v != 0.0)) return QString::number(v, 'g', 4);
    return QString::number(v, 'f', dec);
}

/// tick values (transformed coordinates) and labels for one axis
struct Ticks
{
    QVector<double> t;
    QStringList label;
};

Ticks makeTicks(double a, double b, bool log, int n)
{
    Ticks k;
    if (a > b) std::swap(a, b);
    if (log) {
        QVector<double> vals = {0.0};
        for (double d = 1; d <= 1e9; d *= 10) {
            vals.push_back(d);
            vals.push_back(-d);
        }
        std::sort(vals.begin(), vals.end());
        for (double v : vals) {
            const double t = symlog(v);
            if (t >= a - 1e-12 && t <= b + 1e-12) {
                k.t.push_back(t);
                k.label << (std::abs(v) >= 1e4 ? QString::number(v, 'g', 2) : QString::number(v, 'f', 0));
            }
        }
        if (k.t.size() >= 2) return k;
        k = {};  // too few decades visible: plain ticks
    }
    const double st = niceStep(b - a, n);
    for (double v = std::ceil(a / st - 1e-9) * st; v <= b + 1e-9 * st; v += st) {
        const double val = log ? symexp(v) : v;
        k.t.push_back(v);
        k.label << fmt(val, log ? std::max(st, 1e-3) : st);
    }
    return k;
}

}  // namespace

QColor turbo(double t)
{
    t = std::clamp(t, 0.0, 1.0);
    // Mikhailov (Google, 2019) polynomial fit of the Turbo colour map
    const double r = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))));
    const double g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))));
    const double b = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))));
    auto c = [](double v) { return static_cast<int>(std::clamp(v, 0.0, 1.0) * 255.0 + 0.5); };
    return QColor(c(r), c(g), c(b));
}

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent)
{
    setMouseTracking(true);
    setMinimumSize(260, 180);
    setAutoFillBackground(false);
}

void PlotWidget::setColorBar(bool on, double lo, double hi, const QString& label, std::function<QColor(double)> map)
{
    colorbar_ = on;
    cb_lo_ = lo;
    cb_hi_ = hi;
    cb_label_ = label;
    cb_map_ = std::move(map);
}

QPointF PlotWidget::toData(const QPointF& p) const
{
    const double y = symlogY ? symlog(p.y()) : p.y();
    return swapAxes ? QPointF(y, p.x()) : QPointF(p.x(), y);
}

double PlotWidget::untransformY(double y) const { return symlogY ? symexp(y) : y; }

QPointF PlotWidget::dataAt(const QPointF& screen) const
{
    const QPointF t = fromScreen(screen);
    double x = t.x(), y = t.y();
    if (swapAxes) std::swap(x, y);
    return {x, untransformY(y)};
}

int PlotWidget::handleAt(const QPointF& screen) const
{
    int best = -1;
    double bd = 8.0;  // pick radius [px]
    for (int i = 0; i < handles_.size(); ++i) {
        if (handles_[i].fixed) continue;
        const double d = QLineF(toScreen(toData(handles_[i].p)), screen).length();
        if (d < bd) {
            bd = d;
            best = i;
        }
    }
    return best;
}

int PlotWidget::badgeHeight() const
{
    const QFontMetrics fm(font());
    const int w = std::max(60, width() - 62 - (colorbar_ ? 92 : 16) - 12);
    return fm.boundingRect(QRect(0, 0, w, 400), Qt::TextWordWrap, badge_.text).height() + 8;
}

QRect PlotWidget::plotRect() const
{
    const int left = 62, bottom = 40, right = colorbar_ ? 92 : 16;
    int top = title_.isEmpty() ? 12 : 26;
    if (showBadge && !badge_.text.isEmpty()) top += badgeHeight() + 6;  // the status strip sits above the plot
    return QRect(left, top, std::max(10, width() - left - right), std::max(10, height() - top - bottom));
}

QPointF PlotWidget::toScreen(const QPointF& t) const
{
    const QRectF r = plotRect();
    double fx = (t.x() - x0_) / (x1_ - x0_), fy = (t.y() - y0_) / (y1_ - y0_);
    if (flipX) fx = 1.0 - fx;
    if (flipY) fy = 1.0 - fy;
    return {r.left() + fx * r.width(), r.bottom() - fy * r.height()};
}

QPointF PlotWidget::fromScreen(const QPointF& s) const
{
    const QRectF r = plotRect();
    double fx = (s.x() - r.left()) / r.width(), fy = (r.bottom() - s.y()) / r.height();
    if (flipX) fx = 1.0 - fx;
    if (flipY) fy = 1.0 - fy;
    return {x0_ + fx * (x1_ - x0_), y0_ + fy * (y1_ - y0_)};
}

QRectF PlotWidget::dataBounds() const
{
    double a = std::numeric_limits<double>::infinity(), b = -a, c = a, d = -a;
    auto add = [&](const QPointF& q) {
        const QPointF t = toData(q);
        if (!std::isfinite(t.x()) || !std::isfinite(t.y())) return;
        a = std::min(a, t.x());
        b = std::max(b, t.x());
        c = std::min(c, t.y());
        d = std::max(d, t.y());
    };
    for (const auto& s : series_)
        for (const auto& q : s.pts) add(q);
    for (const auto& f : fills_)
        for (const auto& q : f.poly) add(q);
    if (!std::isfinite(a)) return QRectF(0, 0, 1, 1);
    return QRectF(QPointF(a, c), QPointF(b, d));
}

void PlotWidget::applyProportional()
{
    if (!proportional) return;
    const QRect r = plotRect();
    const double sx = (x1_ - x0_) / r.width(), sy = (y1_ - y0_) / r.height();
    if (sx > sy) {
        const double c = 0.5 * (y0_ + y1_), h = 0.5 * sx * r.height();
        y0_ = c - h;
        y1_ = c + h;
    } else {
        const double c = 0.5 * (x0_ + x1_), h = 0.5 * sy * r.width();
        x0_ = c - h;
        x1_ = c + h;
    }
}

void PlotWidget::zoomFit()
{
    QRectF b = dataBounds();
    double w = b.width(), h = b.height();
    if (!(w > 0)) w = std::max(1e-6, std::abs(b.left()) * 0.1 + 1e-3);
    if (!(h > 0)) h = std::max(1e-6, std::abs(b.top()) * 0.1 + 1e-3);
    x0_ = b.left() - 0.04 * w;
    x1_ = b.left() + 1.04 * w;
    y0_ = b.top() - 0.06 * h;
    y1_ = b.top() + 1.06 * h;
    applyProportional();
    fitted_ = true;
    user_view_ = false;
    update();
}

void PlotWidget::setHandles(QVector<DragHandle> h)
{
    if (drag_handle_ >= 0 && drag_handle_ < handles_.size()) {
        const auto held = handles_[drag_handle_];
        drag_handle_ = -1;
        for (int i = 0; i < h.size(); ++i)
            if (h[i].id == held.id) {
                h[i].p = held.p;
                drag_handle_ = i;
            }
    }
    handles_ = std::move(h);
}

void PlotWidget::commit(bool refit)
{
    if (drag_handle_ >= 0) {  // never rescale under the user's hand
        update();
        return;
    }
    if (refit || !fitted_ || !user_view_) {
        zoomFit();
        return;
    }
    const QRectF b = dataBounds();
    const QRectF v(QPointF(std::min(x0_, x1_), std::min(y0_, y1_)), QPointF(std::max(x0_, x1_), std::max(y0_, y1_)));
    if (!v.intersects(b)) zoomFit();
    else update();
}

void PlotWidget::resizeEvent(QResizeEvent*)
{
    if (proportional && fitted_) applyProportional();
}

void PlotWidget::drawAxes(QPainter& p, const QRect& r)
{
    const bool xlog = symlogY && swapAxes, ylog = symlogY && !swapAxes;
    const auto tx = makeTicks(x0_, x1_, xlog, std::max(3, r.width() / 90));
    const auto ty = makeTicks(y0_, y1_, ylog, std::max(3, r.height() / 45));
    QPen grid(QColor(190, 190, 190), 0, Qt::DotLine);
    grid.setCosmetic(true);
    QPen axis(QColor(60, 60, 60), 0);
    axis.setCosmetic(true);
    const QFontMetrics fm(font());
    p.setPen(grid);
    for (int i = 0; i < tx.t.size(); ++i) {
        const double sx = toScreen({tx.t[i], y0_}).x();
        p.drawLine(QPointF(sx, r.top()), QPointF(sx, r.bottom()));
    }
    for (int i = 0; i < ty.t.size(); ++i) {
        const double sy = toScreen({x0_, ty.t[i]}).y();
        p.drawLine(QPointF(r.left(), sy), QPointF(r.right(), sy));
    }
    // zero lines
    QPen zero(QColor(20, 20, 20), 0);
    zero.setCosmetic(true);
    p.setPen(zero);
    if (x0_ < 0 && x1_ > 0) {
        const double sx = toScreen({0, 0}).x();
        p.drawLine(QPointF(sx, r.top()), QPointF(sx, r.bottom()));
    }
    if (y0_ < 0 && y1_ > 0) {
        const double sy = toScreen({0, 0}).y();
        p.drawLine(QPointF(r.left(), sy), QPointF(r.right(), sy));
    }
    p.setPen(axis);
    p.drawRect(r);
    for (int i = 0; i < tx.t.size(); ++i) {
        const double sx = toScreen({tx.t[i], y0_}).x();
        p.drawLine(QPointF(sx, r.bottom()), QPointF(sx, r.bottom() + 4));
        const int w = fm.horizontalAdvance(tx.label[i]);
        p.drawText(QPointF(sx - w / 2.0, r.bottom() + 6 + fm.ascent()), tx.label[i]);
    }
    for (int i = 0; i < ty.t.size(); ++i) {
        const double sy = toScreen({x0_, ty.t[i]}).y();
        p.drawLine(QPointF(r.left() - 4, sy), QPointF(r.left(), sy));
        const int w = fm.horizontalAdvance(ty.label[i]);
        p.drawText(QPointF(r.left() - 6 - w, sy + fm.ascent() / 2.0 - 1), ty.label[i]);
    }
    const QString xl = swapAxes ? ylabel_ : xlabel_, yl = swapAxes ? xlabel_ : ylabel_;
    p.drawText(QRect(r.left(), r.bottom() + 20, r.width(), 18), Qt::AlignHCenter | Qt::AlignTop, xl);
    p.save();
    p.translate(12, r.center().y());
    p.rotate(-90);
    p.drawText(QRect(-r.height() / 2, -8, r.height(), 18), Qt::AlignCenter, yl);
    p.restore();
}

void PlotWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    if (!fitted_) zoomFit();
    const QRect r = plotRect();
    p.setRenderHint(QPainter::Antialiasing, true);
    if (!title_.isEmpty()) {
        QFont f = font();
        f.setBold(true);
        p.setFont(f);
        p.setPen(Qt::black);
        p.drawText(QRect(0, 4, width(), 20), Qt::AlignHCenter, title_);
        p.setFont(font());
    }
    // fills first (contours), clipped
    p.save();
    p.setClipRect(r);
    p.setRenderHint(QPainter::Antialiasing, false);  // no seams between neighbouring cells
    for (const auto& f : fills_) {
        p.setPen(QPen(f.color, 0));
        QPolygonF s;
        s.reserve(f.poly.size());
        for (const auto& q : f.poly) s << toScreen(toData(q));
        p.setBrush(f.color);
        p.drawPolygon(s);
    }
    p.restore();
    p.setRenderHint(QPainter::Antialiasing, true);
    drawAxes(p, r);
    p.save();
    p.setClipRect(r.adjusted(1, 1, -1, -1));
    for (const auto& s : series_) {
        QPen pen(s.color, s.width, s.style, Qt::RoundCap, Qt::RoundJoin);
        pen.setCosmetic(true);
        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
        QPolygonF line;
        auto flush = [&] {
            if (line.size() >= 2) p.drawPolyline(line);
            line.clear();
        };
        for (const auto& q : s.pts) {
            const QPointF t = toData(q);
            if (!std::isfinite(t.x()) || !std::isfinite(t.y())) {
                flush();
                continue;
            }
            line << toScreen(t);
        }
        flush();
        if (showPoints || s.markers) {
            const double m = s.marker_size;
            for (const auto& q : s.pts) {
                const QPointF t = toData(q);
                if (!std::isfinite(t.x()) || !std::isfinite(t.y())) continue;
                const QPointF c = toScreen(t);
                p.drawRect(QRectF(c.x() - m, c.y() - m, 2 * m, 2 * m));
            }
        }
    }
    for (const auto& m : markers_) {
        const QPointF c = toScreen(toData(m.p));
        QPen pen(m.color, m.width);
        pen.setCosmetic(true);
        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
        const double z = m.size;
        if (m.shape == 'x') {
            p.drawLine(c + QPointF(-z, -z), c + QPointF(z, z));
            p.drawLine(c + QPointF(-z, z), c + QPointF(z, -z));
        } else if (m.shape == '^') {
            p.drawPolygon(QPolygonF({c + QPointF(0, -1.15 * z), c + QPointF(z, 0.8 * z), c + QPointF(-z, 0.8 * z)}));
        } else if (m.shape == 's') {
            p.drawRect(QRectF(c.x() - 0.7 * z, c.y() - 0.7 * z, 1.4 * z, 1.4 * z));
        } else {
            p.drawEllipse(c, 0.8 * z, 0.8 * z);
        }
        if (!m.label.isEmpty()) p.drawText(c + QPointF(8, -6), m.label);
    }
    for (int i = 0; i < handles_.size(); ++i) {
        const auto& h = handles_[i];
        const QPointF c = toScreen(toData(h.p));
        QPen pen(h.color, i == drag_handle_ ? 2.0 : 1.2);
        pen.setCosmetic(true);
        p.setPen(pen);
        p.setBrush(h.fixed ? QBrush(Qt::white) : QBrush(h.color.lighter(170)));
        if (h.fixed) p.drawRect(QRectF(c.x() - 2.5, c.y() - 2.5, 5, 5));
        else p.drawEllipse(c, 4.5, 4.5);
    }
    p.restore();

    const QFontMetrics fm(font());
    // legend (top right)
    if (showLegend) {
        QStringList names;
        QVector<const Series*> items;
        for (const auto& s : series_)
            if (s.legend && !s.name.isEmpty() && !names.contains(s.name)) {
                names << s.name;
                items << &s;
            }
        if (!items.isEmpty()) {
            int w = 0;
            for (const auto& n : names) w = std::max(w, fm.horizontalAdvance(n));
            const int lh = fm.height(), bw = w + 40, bh = lh * static_cast<int>(items.size()) + 8;
            const QRect box(r.right() - bw - 6, r.top() + 6, bw, bh);
            p.setPen(QPen(QColor(150, 150, 150), 0));
            p.setBrush(QColor(255, 255, 255, 225));
            p.drawRect(box);
            for (int i = 0; i < items.size(); ++i) {
                const int y = box.top() + 4 + lh * i + lh / 2;
                QPen pen(items[i]->color, items[i]->width, items[i]->style);
                pen.setCosmetic(true);
                p.setPen(pen);
                p.drawLine(box.left() + 6, y, box.left() + 30, y);
                p.setPen(Qt::black);
                p.drawText(box.left() + 36, y + fm.ascent() / 2 - 1, names[i]);
            }
        }
    }
    // status strip (above the plot, never over the data or the legend)
    if (showBadge && !badge_.text.isEmpty()) {
        const QColor bg = badge_.severity == 2 ? QColor(255, 220, 220) : badge_.severity == 1 ? QColor(255, 242, 214) : QColor(222, 245, 222);
        const QColor bd = badge_.severity == 2 ? QColor(200, 0, 0) : badge_.severity == 1 ? QColor(220, 130, 0) : QColor(0, 140, 0);
        const int bh = badgeHeight();
        const QRect box(r.left(), r.top() - bh - 6, r.width(), bh);
        p.setPen(QPen(bd, 1.2));
        p.setBrush(bg);
        p.drawRoundedRect(box, 4, 4);
        p.setPen(Qt::black);
        p.drawText(box.adjusted(6, 4, -6, -4), Qt::TextWordWrap, badge_.text);
    }
    // colour bar
    if (colorbar_ && cb_map_) {
        const QRect cb(r.right() + 14, r.top() + 4, 16, r.height() - 8);
        for (int i = 0; i < cb.height(); ++i) {
            const double t = 1.0 - static_cast<double>(i) / std::max(1, cb.height() - 1);
            p.setPen(cb_map_(t));
            p.drawLine(cb.left(), cb.top() + i, cb.right(), cb.top() + i);
        }
        p.setPen(QPen(Qt::black, 0));
        p.setBrush(Qt::NoBrush);
        p.drawRect(cb);
        const auto tk = makeTicks(cb_lo_, cb_hi_, false, 6);
        for (int i = 0; i < tk.t.size(); ++i) {
            const double f = (tk.t[i] - cb_lo_) / std::max(cb_hi_ - cb_lo_, 1e-300);
            const int y = static_cast<int>(cb.bottom() - f * cb.height());
            p.drawLine(cb.right(), y, cb.right() + 3, y);
            p.drawText(cb.right() + 5, y + fm.ascent() / 2 - 1, tk.label[i]);
        }
        p.save();
        p.translate(width() - 8, cb.center().y());
        p.rotate(-90);
        p.drawText(QRect(-cb.height() / 2, -12, cb.height(), 14), Qt::AlignCenter, cb_label_);
        p.restore();
    }
    // cursor position (bottom left)
    if (r.contains(cursor_.toPoint())) {
        const QPointF t = fromScreen(cursor_);
        double x = t.x(), y = t.y();
        if (swapAxes) std::swap(x, y);
        y = untransformY(y);
        p.setPen(Qt::black);
        p.drawText(4, height() - 4, QString("(%1, %2)").arg(x, 0, 'g', 6).arg(y, 0, 'g', 6));
    }
}

void PlotWidget::wheelEvent(QWheelEvent* e)
{
    const double f = std::pow(0.85, e->angleDelta().y() / 120.0);
    const QPointF c = fromScreen(e->position());
    x0_ = c.x() + (x0_ - c.x()) * f;
    x1_ = c.x() + (x1_ - c.x()) * f;
    y0_ = c.y() + (y0_ - c.y()) * f;
    y1_ = c.y() + (y1_ - c.y()) * f;
    user_view_ = true;
    update();
}

void PlotWidget::mousePressEvent(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton) {
        drag_handle_ = handleAt(e->position());
        if (drag_handle_ >= 0) {
            setCursor(Qt::ClosedHandCursor);
            update();
            return;
        }
    }
    if (e->button() == Qt::LeftButton || e->button() == Qt::MiddleButton) {
        dragging_ = true;
        drag_start_ = e->position();
        dx0_ = x0_;
        dx1_ = x1_;
        dy0_ = y0_;
        dy1_ = y1_;
    }
}

void PlotWidget::mouseMoveEvent(QMouseEvent* e)
{
    cursor_ = e->position();
    if (drag_handle_ >= 0) {
        handles_[drag_handle_].p = dataAt(e->position());
        emit handleMoved(handles_[drag_handle_].id, handles_[drag_handle_].p, false);
        update();
        return;
    }
    if (!(e->buttons() & (Qt::LeftButton | Qt::MiddleButton)))
        setCursor(handleAt(e->position()) >= 0 ? Qt::OpenHandCursor : Qt::CrossCursor);
    if (dragging_) {
        const QRectF r = plotRect();
        double ddx = (e->position().x() - drag_start_.x()) / r.width() * (dx1_ - dx0_);
        double ddy = (e->position().y() - drag_start_.y()) / r.height() * (dy1_ - dy0_);
        if (flipX) ddx = -ddx;
        if (flipY) ddy = -ddy;
        x0_ = dx0_ - ddx;
        x1_ = dx1_ - ddx;
        y0_ = dy0_ + ddy;
        y1_ = dy1_ + ddy;
        user_view_ = true;
    }
    update();
}

void PlotWidget::mouseReleaseEvent(QMouseEvent* e)
{
    if (drag_handle_ >= 0) {
        const int i = drag_handle_;
        drag_handle_ = -1;
        handles_[i].p = dataAt(e->position());
        emit handleMoved(handles_[i].id, handles_[i].p, true);
        setCursor(Qt::OpenHandCursor);
        update();
    }
    dragging_ = false;
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent*) { zoomFit(); }

void PlotWidget::contextMenuEvent(QContextMenuEvent* e) { emit contextMenuWanted(e->globalPos()); }

QString PlotWidget::csv() const
{
    QString out = QString("# %1\nseries,%2,%3\n").arg(title_, xlabel_, ylabel_);
    for (const auto& s : series_)
        for (const auto& q : s.pts)
            if (std::isfinite(q.x()) && std::isfinite(q.y()))
                out += QString("%1,%2,%3\n").arg(s.name).arg(q.x(), 0, 'g', 12).arg(q.y(), 0, 'g', 12);
    return out;
}

QImage PlotWidget::image()
{
    QImage img(size() * devicePixelRatioF(), QImage::Format_ARGB32);
    img.setDevicePixelRatio(devicePixelRatioF());
    img.fill(Qt::white);
    render(&img);
    return img;
}

}  // namespace Designer
