// View3D: the meshed 3-D blade, hub and casing (QPainter, orthographic, no OpenGL).
// Shaded mode sorts the quads and grid lines back to front (painter's algorithm) with
// Lambert shading; wireframe mode draws the grid lines only.  Mouse: left drag rotates,
// middle drag (or shift + left) pans, wheel zooms, double click fits.
#pragma once

#include <QColor>
#include <QString>
#include <QVector>
#include <QWidget>

#include <array>
#include <vector>

namespace Designer {

struct Grid3D
{
    QString name;
    int ni = 0, nj = 0;
    std::vector<std::array<double, 3>> p;  ///< i fastest
    QColor color;
    bool surface = false;                  ///< shaded quads in shaded mode (else lines only)
    bool points = false;                   ///< node dots
};

struct Line3D
{
    std::vector<std::array<double, 3>> p;
    QColor color;
    Qt::PenStyle style = Qt::SolidLine;
    double width = 1.0;
};

class View3D : public QWidget
{
    Q_OBJECT
public:
    explicit View3D(QWidget* parent = nullptr);

    void setScene(std::vector<Grid3D> grids, std::vector<Line3D> lines, bool refit);
    bool shaded = true;
    bool showPoints = false;
    void zoomFit();
    void setView(double yaw_deg, double pitch_deg);
    [[nodiscard]] QImage image();

signals:
    void contextMenuWanted(const QPoint& global);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;
    void wheelEvent(QWheelEvent*) override;
    void contextMenuEvent(QContextMenuEvent*) override;

private:
    [[nodiscard]] std::array<double, 3> view(const std::array<double, 3>& p) const;  // view coords (x right, y up, z toward viewer)
    [[nodiscard]] QPointF screen(const std::array<double, 3>& v) const;

    std::vector<Grid3D> grids_;
    std::vector<Line3D> lines_;
    std::array<double, 3> center_{0, 0, 0};
    double extent_ = 1.0;
    double yaw_ = -35.0, pitch_ = 25.0, zoom_ = 1.0;
    QPointF pan_{0, 0};
    bool fitted_ = false;
    int drag_ = 0;  // 1 rotate, 2 pan
    QPointF last_;
};

}  // namespace Designer
