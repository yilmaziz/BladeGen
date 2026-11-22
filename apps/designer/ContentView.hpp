// ContentView: one MDI child.  Its content (3-D view, OpenCASCADE solid view, blade-to-blade
// view, meridional contour, angle / thickness / curvature plots, information tables) is chosen
// from the right-click menu, together with the layer, the horizontal axis of the angle /
// thickness plots and the plot options.  The solid view builds its OCCT shapes on a worker
// thread whenever the design result changes (one build at a time; a newer result waits).
#pragma once

#include "Engine.hpp"
#include "design/DesignSolid.hpp"

#include <QFutureWatcher>
#include <QStackedWidget>
#include <QWidget>

#include <memory>

class QTableWidget;
class QMenu;

namespace Designer {

class OcctView;
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
    Solid3D,  ///< the OCCT blade solid, hub and casing (single blade or whole annulus)
};

[[nodiscard]] QString contentName(Content c);
[[nodiscard]] QList<QList<Content>> contentGroups();  ///< for the menus (separators between groups)

enum class HAxis { MeridionalFraction, MPrime, Meridional, ChordFraction };
/// Blade-to-blade view coordinates: conformal (m', theta); (m, r theta); (z, r theta); the axial view
/// (x, y) = (r sin theta, r cos theta), the section seen along the machine axis from downstream.
enum class B2BCoords { MPrimeTheta, MRTheta, ZRTheta, XYAxial };
enum class CurvatureRange { Full, LeadingEdge, TrailingEdge };

class ContentView : public QWidget
{
    Q_OBJECT
public:
    ContentView(Session& s, Content c, QWidget* parent = nullptr);
    ~ContentView() override;

    [[nodiscard]] Content content() const { return content_; }
    void setContent(Content c);
    /// Redraw from the session (refit: reset the zoom).
    void refresh(bool refit = false);
    [[nodiscard]] QImage image();
    /// Profile Section drag targets: kSsHandle + SS control index, kPsHandle + PS control index.
    static constexpr int kSsHandle = 1000, kPsHandle = 2000;
    /// As if the user dragged control point `id` to p (self-test).
    void simulateDrag(int id, QPointF p) { dragControlPoint(id, p, true); }
    /// 3D Solid: every blade around the annulus (else one blade on one pitch of the walls).
    void setAnnulus(bool on);
    /// 3D Solid: the exact OCCT solid (trimmed B-spline solid, built on a worker thread in about
    /// 2 s) instead of the fast display mesh of the design rows (updated with every change).
    void setExact(bool on);
    [[nodiscard]] bool exact() const { return exact_; }
    [[nodiscard]] bool annulus() const { return annulus_; }
    [[nodiscard]] OcctView* occtView() const { return occt_; }
    /// 3D Solid: waits for the solid being built (tests); its error, empty when it was built.
    QString waitForSolid();
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
    void buildSolid();
    void startSolid();
    void solidDone();

    Session& s_;
    Content content_;
    QStackedWidget* stack_;
    PlotWidget* plot_;
    View3D* view3d_;
    QTableWidget* table_;
    OcctView* occt_;
    bool show_walls_ = true;
    // 3D Solid
    struct SolidJob
    {
        std::shared_ptr<const BladeGen::Design::Scene> scene;
        QString error;
        double seconds = 0.0;
        bool annulus = false;
    };
    QFutureWatcher<SolidJob>* solid_watch_;
    std::unique_ptr<BladeGen::Design::Evaluator> solid_ev_;  ///< used by the worker only
    std::shared_ptr<const BladeGen::Design::DesignResult> solid_for_;
    QString solid_error_;
    bool exact_ = false;
    bool annulus_ = false, solid_annulus_ = false, solid_pending_ = false, solid_applied_ = true, solid_shown_ = false, solid_refit_ = false;
};

}  // namespace Designer
