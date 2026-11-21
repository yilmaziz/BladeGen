// ContentView: one MDI child.  Its content (3-D view, blade-to-blade view, meridional
// contour, angle / thickness / curvature plots, information tables) is chosen from the
// right-click menu, together with the layer, the horizontal axis of the angle / thickness
// plots and the plot options.
#pragma once

#include "Engine.hpp"

#include <QStackedWidget>
#include <QWidget>

class QTableWidget;
class QMenu;

namespace Designer {

class PlotWidget;
class View3D;

enum class Content
{
    View3D,
    BladeToBlade,
    ProfileSection,
    MeridionalContour,
    MeridionalBeta,
    MeridionalThickness,
    MeridionalCurvature,
    ProfileCurvature,
    LeTeTheta,
    LeTeBeta,
    LeTeParameters,
    AirfoilArea,
    InfoCurrent,
    InfoAll,
    AnglesCurrent,
    BetaAll,
    ThetaAll,
    BetaVsTheta,
    ThicknessCurrent,
    ThicknessAll,
};

[[nodiscard]] QString contentName(Content c);
[[nodiscard]] QList<QList<Content>> contentGroups();  ///< for the menus (separators between groups)

enum class HAxis { MeridionalFraction, MPrime, Meridional, ChordFraction };
enum class B2BCoords { MPrimeTheta, MRTheta };
enum class CurvatureRange { Full, LeadingEdge, TrailingEdge };

class ContentView : public QWidget
{
    Q_OBJECT
public:
    ContentView(Session& s, Content c, QWidget* parent = nullptr);

    [[nodiscard]] Content content() const { return content_; }
    void setContent(Content c);
    /// Redraw from the session (refit: reset the zoom).
    void refresh(bool refit = false);
    [[nodiscard]] QImage image();
    /// Profile Section drag targets: kSsHandle + SS control index, kPsHandle + PS control index.
    static constexpr int kSsHandle = 1000, kPsHandle = 2000;
    /// As if the user dragged control point `id` to p (self-test).
    void simulateDrag(int id, QPointF p) { dragControlPoint(id, p, true); }
    int layer = -1;  ///< -1: follow the active layer
    HAxis haxis = HAxis::MeridionalFraction;
    B2BCoords b2b = B2BCoords::MPrimeTheta;
    CurvatureRange crange = CurvatureRange::Full;

signals:
    void newViewWanted(Content c);
    void meshPointsWanted();
    void titleChanged(const QString& t);
    /// The view changed the design (a dragged control point); `done` on release.
    void designEdited(bool done);

private:
    void showMenu(const QPoint& global);
    void dragControlPoint(int id, QPointF p, bool done);
    [[nodiscard]] int layerIndex() const;
    void buildPlot(bool refit);
    void build3D(bool refit);
    void buildTable();
    void saveData();

    Session& s_;
    Content content_;
    QStackedWidget* stack_;
    PlotWidget* plot_;
    View3D* view3d_;
    QTableWidget* table_;
    bool show_walls_ = true;
};

}  // namespace Designer
