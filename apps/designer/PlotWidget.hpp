// PlotWidget: a light 2-D plot drawn with QPainter (thin cosmetic lines, antialiased).
// Series are polylines (a NaN point breaks a line), optional point markers, issue markers,
// filled polygons (contour plots) with a colour bar, a status badge, legend, axis labels.
// Interaction: wheel = zoom at the cursor, left drag = pan, double click = zoom fit; the
// cursor position is shown bottom left.  View options: legend, points, swap axes,
// proportional axes, flip x / y, symmetric-log y (curvature).
#pragma once

#include <QColor>
#include <QPointF>
#include <QPolygonF>
#include <QString>
#include <QVector>
#include <QWidget>

#include <functional>

namespace Designer {

struct Series
{
    QString name;
    QVector<QPointF> pts;
    QColor color = Qt::black;
    double width = 1.0;
    Qt::PenStyle style = Qt::SolidLine;
    bool legend = true;
    bool markers = false;      ///< draw point markers even when "show points" is off
    int marker_size = 3;
};

struct PointMarker
{
    QPointF p;
    QColor color = Qt::red;
    char shape = 'x';          ///< 'x', '^', 'o', 's'
    QString label;
};

struct Fill
{
    QPolygonF poly;
    QColor color;
};

struct Badge
{
    QString text;
    int severity = 0;          ///< 0 ok, 1 warning, 2 error
};

class PlotWidget : public QWidget
{
    Q_OBJECT
public:
    explicit PlotWidget(QWidget* parent = nullptr);

    void setTitle(const QString& t) { title_ = t; }
    void setLabels(const QString& x, const QString& y) { xlabel_ = x; ylabel_ = y; }
    void setSeries(QVector<Series> s) { series_ = std::move(s); }
    void setMarkers(QVector<PointMarker> m) { markers_ = std::move(m); }
    void setFills(QVector<Fill> f) { fills_ = std::move(f); }
    void setColorBar(bool on, double lo = 0, double hi = 1, const QString& label = {}, std::function<QColor(double)> map = {});
    void setBadge(const Badge& b) { badge_ = b; }
    /// Repaint with the new data; the view is refitted when `refit` or when the data
    /// left the current view entirely.
    void commit(bool refit = false);

    // view options
    bool showLegend = true, showPoints = false, swapAxes = false, proportional = false, flipX = false, flipY = false;
    bool symlogY = false;      ///< y' = sign(y) log10(1 + |y|) (linear below 1)
    bool showBadge = true;     ///< status box (top left)

    void zoomFit();
    [[nodiscard]] QString csv() const;  ///< every series as name,x,y rows
    [[nodiscard]] QImage image();

signals:
    void contextMenuWanted(const QPoint& global);

protected:
    void paintEvent(QPaintEvent*) override;
    void wheelEvent(QWheelEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;
    void contextMenuEvent(QContextMenuEvent*) override;
    void resizeEvent(QResizeEvent*) override;

private:
    [[nodiscard]] QPointF toData(const QPointF& p) const;      // data (after swap / symlog) -> transformed
    [[nodiscard]] QPointF toScreen(const QPointF& t) const;    // transformed -> pixels
    [[nodiscard]] QPointF fromScreen(const QPointF& s) const;  // pixels -> transformed
    [[nodiscard]] double untransformY(double y) const;
    [[nodiscard]] QRectF dataBounds() const;
    [[nodiscard]] QRect plotRect() const;
    [[nodiscard]] int badgeHeight() const;
    void applyProportional();
    void drawAxes(QPainter& p, const QRect& r);

    QString title_, xlabel_, ylabel_;
    QVector<Series> series_;
    QVector<PointMarker> markers_;
    QVector<Fill> fills_;
    Badge badge_;
    bool colorbar_ = false;
    double cb_lo_ = 0, cb_hi_ = 1;
    QString cb_label_;
    std::function<QColor(double)> cb_map_;
    double x0_ = 0, x1_ = 1, y0_ = 0, y1_ = 1;   // visible range in transformed coordinates
    bool fitted_ = false;
    bool user_view_ = false;  // the user zoomed / panned: keep the view on data updates
    bool dragging_ = false;
    QPointF drag_start_;
    double dx0_ = 0, dx1_ = 0, dy0_ = 0, dy1_ = 0;
    QPointF cursor_{-1, -1};
};

/// Perceptual "turbo"-like colour map on [0, 1] (Mikhailov 2019 polynomial approximation).
[[nodiscard]] QColor turbo(double t);

}  // namespace Designer
