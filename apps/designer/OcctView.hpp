// OcctView: the blade, hub and casing in OpenCASCADE's own viewer (V3d / AIS, OpenGL).  Two
// sources: the fast display mesh of the design rows (Design::SceneMesh, our own triangles drawn
// as OpenGL primitive arrays: no OCCT modelling, updated with every design change) or the exact
// OCCT solid (Design::Scene: trimmed B-spline solid, meshed by BRepMesh; slow).
// The view renders offscreen (V3d_View::ToPixMap into a frame buffer of the widget's size) and
// paints the image, so it behaves like any other widget inside the MDI area and its image can
// be saved.  Mouse: left drag rotates, middle drag (or shift + left) pans, wheel zooms, double
// click fits, right click opens the view menu.
//
// Each part (blade, hub, casing) has a look: hidden, wireframe (edges and iso lines), solid
// (opaque, lit, with its face boundaries) or shaded (lit, smooth, see-through).  Without the OCCT
// visualization toolkits, or without an OpenGL context (a headless test machine), the view shows
// why instead.  BLADEGEN_NO_OPENGL (set to anything) switches the view off.
#pragma once

#include "design/DesignSolid.hpp"

#include <QWidget>

#include <memory>

namespace Designer {

class OcctView : public QWidget
{
    Q_OBJECT
public:
    enum class Part { Blade, Hub, Casing };
    enum class Look { Hidden, Wireframe, Solid, Shaded };
    enum class Direction { Iso, Axial, Meridional, Radial };

    explicit OcctView(QWidget* parent = nullptr);
    ~OcctView() override;

    /// The viewer works (toolkits built in and an OpenGL context created).
    [[nodiscard]] bool available() const;
    /// Why it does not (empty when it does).
    [[nodiscard]] QString unavailableReason() const;
    /// Meshes the scene for display: call it off the GUI thread, before setScene.
    static void prepare(const BladeGen::Design::Scene& s);
    void setScene(std::shared_ptr<const BladeGen::Design::Scene> s, bool refit);
    /// The fast display mesh (cheap: call it on the GUI thread with every result).
    void setMeshScene(std::shared_ptr<const BladeGen::Design::SceneMesh> m, bool refit);
    /// Text drawn over the view (progress, errors); empty: none.
    void setMessage(const QString& m);
    void setLook(Part p, Look l);
    [[nodiscard]] Look look(Part p) const;
    void setDirection(Direction d);
    void fitAll();
    [[nodiscard]] QImage image();

signals:
    void contextMenuWanted(const QPoint& global);

protected:
    void paintEvent(QPaintEvent*) override;
    void resizeEvent(QResizeEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;
    void wheelEvent(QWheelEvent*) override;
    void contextMenuEvent(QContextMenuEvent*) override;

private:
    struct Impl;
    std::unique_ptr<Impl> d_;
    QString message_;
    QPoint press_;
    bool moved_ = false;
};

}  // namespace Designer
