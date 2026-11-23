#include "ContentView.hpp"

#include "OcctView.hpp"
#include "PlotWidget.hpp"
#include "View3D.hpp"
#include "design/DesignSolid.hpp"

#include <Standard_Failure.hxx>

#include <QActionGroup>
#include <QApplication>
#include <QClipboard>
#include <QFile>
#include <QFileDialog>
#include <QHeaderView>
#include <QMenu>
#include <QMessageBox>
#include <QTableWidget>
#include <QTextStream>
#include <QVBoxLayout>
#include <QtConcurrent/QtConcurrentRun>

#include <chrono>
#include <cmath>

namespace Designer {

using namespace BladeGen;
using Design::LayerResult;

namespace {

const QColor kSS(25, 80, 225), kPS(220, 30, 30), kEdge(0, 150, 60), kCamber(120, 120, 120);
const QColor kLayer[] = {QColor(30, 60, 235), QColor(0, 175, 195), QColor(0, 150, 0), QColor(120, 170, 20), QColor(225, 30, 30)};

QColor layerColor(std::size_t k) { return kLayer[k % 5]; }

QColor pieceColor(int piece, bool suction)
{
    if (piece == static_cast<int>(Profile::Piece::LeadingEdge) || piece == static_cast<int>(Profile::Piece::TrailingEdge)) return kEdge;
    return suction ? kSS : kPS;
}

/// Polyline split into runs of equal piece (shared end points), coloured by piece.
void addPieces(QVector<Series>& out, const std::vector<QPointF>& pts, const std::vector<int>& piece, bool suction_default,
               double width = 1.0, bool legend = true)
{
    if (pts.empty()) return;
    std::size_t a = 0;
    bool first_ss = true, first_ps = true, first_edge = true;
    // in a closed loop the side changes at the TE: pieces after the TE belong to the pressure side
    bool suction = suction_default;
    for (std::size_t i = 1; i <= pts.size(); ++i) {
        if (i < pts.size() && piece[i] == piece[a]) continue;
        Series s;
        const bool edge = piece[a] == static_cast<int>(Profile::Piece::LeadingEdge) || piece[a] == static_cast<int>(Profile::Piece::TrailingEdge);
        if (piece[a] == static_cast<int>(Profile::Piece::SuctionSide)) suction = true;
        if (piece[a] == static_cast<int>(Profile::Piece::PressureSide)) suction = false;
        s.color = pieceColor(piece[a], suction);
        s.width = width;
        s.name = edge ? "LE / TE closure" : suction ? "suction side" : "pressure side";
        s.legend = legend && (edge ? first_edge : suction ? first_ss : first_ps);
        (edge ? first_edge : suction ? first_ss : first_ps) = false;
        for (std::size_t k = a; k < std::min(i + 1, pts.size()); ++k) s.pts << pts[k];
        out << s;
        a = i;
    }
}

/// Geometry-check location: a small x, red for an error, black for a warning.
PointMarker issueMarker(const QPointF& p, Profile::Severity sev)
{
    PointMarker m;
    m.p = p;
    m.color = sev == Profile::Severity::Error ? QColor(220, 0, 0) : QColor(0, 0, 0);
    m.shape = 'x';
    m.size = 3.5;
    m.width = 1.3;
    return m;
}

/// The custom thickness law with the layer's numbers (ProfileParams.hpp, ThicknessParams).
QString thicknessFormula(const Profile::ProfileParams& p, double xi_max)
{
    const auto& t = p.thickness;
    return QString("thickness law along the SS arc fraction ξ (ξ_max = %1 at m′ = %2):\n"
                   "ξ ≤ ξ_max:  t = %3 + (%4 − %3) (ξ / ξ_max)^%5\n"
                   "ξ > ξ_max:  t = %6 + (%4 − %6) ((1 − ξ) / (1 − ξ_max))^%7%8")
        .arg(xi_max, 0, 'f', 4)
        .arg(t.mprime_max, 0, 'f', 4)
        .arg(t.t_le, 0, 'f', 5)
        .arg(t.t_max, 0, 'f', 5)
        .arg(t.shape_le, 0, 'f', 4)
        .arg(t.t_te, 0, 'f', 5)
        .arg(t.shape_te, 0, 'f', 4)
        .arg(t.enforce_c1 ? QString("   (shape_te set by the C1 option)") : QString());
}

QString num(double v, int prec = 4) { return QString::number(v, 'g', prec); }

struct Conf
{
    double c, cg, sg, m0, th, sign;
    [[nodiscard]] QPointF operator()(const Math::Vec2& q) const
    {
        return {m0 + c * (q.x() * cg - sign * q.y() * sg), th + c * (q.x() * sg + sign * q.y() * cg)};
    }
};

Conf confOf(const LayerResult& L, bool mirror)
{
    return {L.chord_mprime, std::cos(L.gamma), std::sin(L.gamma), L.mprime_le, L.theta_le, mirror ? -1.0 : 1.0};
}

Badge badgeFor(const Profile::QualityReport& q, const QString& head)
{
    const int ne = q.count(Profile::Severity::Error), nw = q.count(Profile::Severity::Warning);
    Badge b;
    b.severity = ne ? 2 : nw ? 1 : 0;
    b.text = head + "\n" + (ne ? QString("✗ %1 geometry error(s): not manufacturable as designed").arg(ne)
                            : nw ? QString("⚠ %1 warning(s): geometry built (list: Geometry checks)").arg(nw)
                                 : QString("✓ geometry OK: C2 / G2 junctions, no edge defects"));
    for (const auto& i : q.issues)  // the first error in full; the rest is in the Geometry checks dock
        if (i.severity == Profile::Severity::Error) {
            b.text += "\n" + QString::fromStdString(i.message).left(110);
            break;
        }
    return b;
}

double edgeWindow(const Profile::QualityReport& q, bool le)
{
    double w = 0.0;
    for (const auto* sd : {&q.ss, &q.ps}) {
        for (std::size_t i = 1; i < sd->piece.size(); ++i)
            if (sd->piece[i] != sd->piece[i - 1]) {
                const double f = sd->s[i] / sd->s.back();
                if (le) {
                    w = std::max(w, f);
                    break;
                }
            }
        if (!le)
            for (std::size_t i = sd->piece.size() - 1; i > 0; --i)
                if (sd->piece[i] != sd->piece[i - 1]) {
                    w = std::max(w, 1.0 - sd->s[i] / sd->s.back());
                    break;
                }
    }
    return std::clamp(3.0 * w, 0.004, 0.5);
}

}  // namespace

QString contentName(Content c)
{
    switch (c) {
    case Content::View3D: return "3D View";
    case Content::BladeToBlade: return "Blade-to-Blade View";
    case Content::ProfileSection: return "Profile Section (chord frame)";
    case Content::MeridionalContour: return "Meridional Contour View";
    case Content::MeridionalBeta: return "Meridional Metal Angle Contour";
    case Content::MeridionalThickness: return "Meridional Thickness Contour";
    case Content::MeridionalCurvature: return "Meridional Curvature";
    case Content::ProfileCurvature: return "Blade-to-Blade Curvature (SS / PS)";
    case Content::LeTeTheta: return "LE/TE Theta";
    case Content::LeTeBeta: return "LE/TE Beta";
    case Content::LeTeParameters: return "LE/TE Parameters";
    case Content::AirfoilArea: return "Airfoil Area";
    case Content::InfoCurrent: return "Blade Information Table, Current Layer";
    case Content::InfoAll: return "Blade Information Table, All Layers";
    case Content::AnglesCurrent: return "Blade Angles, Current Layer";
    case Content::BetaAll: return "Beta, All Layers";
    case Content::ThetaAll: return "Theta, All Layers";
    case Content::BetaVsTheta: return "Beta vs. Theta, All Layers";
    case Content::ThicknessCurrent: return "Thickness, Current Layer";
    case Content::ThicknessAll: return "Thickness, All Layers";
    case Content::Solid3D: return "3D Solid (OpenCASCADE)";
    case Content::BladeSections: return "Blade Sections (all layers)";
    }
    return "?";
}

QList<QList<Content>> contentGroups()
{
    return {{Content::View3D, Content::Solid3D, Content::BladeToBlade, Content::BladeSections, Content::ProfileSection, Content::MeridionalContour},
            {Content::MeridionalBeta, Content::MeridionalThickness},
            {Content::MeridionalCurvature, Content::ProfileCurvature},
            {Content::LeTeTheta, Content::LeTeBeta, Content::LeTeParameters},
            {Content::AirfoilArea},
            {Content::InfoCurrent, Content::InfoAll},
            {Content::AnglesCurrent, Content::BetaAll, Content::ThetaAll, Content::BetaVsTheta},
            {Content::ThicknessCurrent, Content::ThicknessAll}};
}

ContentView::ContentView(Session& s, Content c, QWidget* parent) : QWidget(parent), s_(s), content_(c)
{
    auto* lay = new QVBoxLayout(this);
    lay->setContentsMargins(0, 0, 0, 0);
    stack_ = new QStackedWidget(this);
    plot_ = new PlotWidget(stack_);
    view3d_ = new View3D(stack_);
    table_ = new QTableWidget(stack_);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setContextMenuPolicy(Qt::CustomContextMenu);
    table_->verticalHeader()->setDefaultSectionSize(20);
    occt_ = new OcctView(stack_);
    solid_watch_ = new QFutureWatcher<SolidJob>(this);
    solid_ev_ = std::make_unique<Design::Evaluator>();
    connect(solid_watch_, &QFutureWatcher<SolidJob>::finished, this, &ContentView::solidDone);
    stack_->addWidget(plot_);
    stack_->addWidget(view3d_);
    stack_->addWidget(table_);
    stack_->addWidget(occt_);
    lay->addWidget(stack_);
    connect(plot_, &PlotWidget::contextMenuWanted, this, &ContentView::showMenu);
    connect(view3d_, &View3D::contextMenuWanted, this, &ContentView::showMenu);
    connect(occt_, &OcctView::contextMenuWanted, this, &ContentView::showMenu);
    connect(plot_, &PlotWidget::handleMoved, this, &ContentView::dragControlPoint);
    connect(table_, &QTableWidget::customContextMenuRequested, this, [this](const QPoint& p) { showMenu(table_->viewport()->mapToGlobal(p)); });
    setContent(c);
}

ContentView::~ContentView()
{
    if (solid_watch_->isRunning()) solid_watch_->waitForFinished();  // the worker uses solid_ev_
}

void ContentView::setContent(Content c)
{
    content_ = c;
    const bool shape = c == Content::BladeToBlade || c == Content::BladeSections || c == Content::ProfileSection || c == Content::MeridionalContour ||
                       c == Content::MeridionalBeta || c == Content::MeridionalThickness;
    plot_->proportional = shape;
    plot_->swapAxes = plot_->flipX = plot_->flipY = false;
    plot_->symlogY = c == Content::ProfileCurvature;
    plot_->showLegend = true;
    refresh(true);
}

void ContentView::dragControlPoint(int id, QPointF q, bool done)
{
    const int k = layerIndex();
    if (k >= static_cast<int>(s_.design.layers.size())) return;
    auto& lay = s_.design.layers[static_cast<std::size_t>(k)];
    if (!lay.uses_custom()) return;
    auto& p = lay.custom;
    try {
        const Profile::Profile2D prof(p);
        const Math::Vec2 P(q.x(), q.y());
        if (id >= kSsHandle && id < kPsHandle) {
            const auto& c = prof.suction_side().ctrl();
            const int i = id - kSsHandle;
            const double span = (c.back() - c.front()).norm();
            if (i == 2) {
                p.ss.ctrl2_mprime = P.x();
                p.ss.ctrl2_theta = P.y();
            } else if (i == 1) {
                const Math::Vec2 v = Profile::unit_from_beta(p.le.beta_deg + p.le.d_beta_ss_deg, +1.0);
                p.ss.handle_le_frac = std::max(0.005, (P - c.front()).dot(v) / span);
            } else if (i == 3) {
                const Math::Vec2 v = Profile::unit_from_beta(p.te.beta_deg + p.te.d_beta_ss_deg, -1.0);
                p.ss.handle_te_frac = std::max(0.005, (P - c.back()).dot(v) / span);
            }
        } else if (id >= kPsHandle) {
            const auto& c = prof.pressure_side().ctrl();
            const int i = id - kPsHandle, n = static_cast<int>(c.size());
            const double span = (c.back() - c.front()).norm();
            if (i == 1) {
                const Math::Vec2 v = Profile::unit_from_beta(p.le.beta_deg - p.le.d_beta_ps_deg, +1.0);
                p.ps.handle_le_frac = std::max(0.005, (P - c.front()).dot(v) / span);
            } else if (i == n - 2) {
                const Math::Vec2 v = Profile::unit_from_beta(p.te.beta_deg - p.te.d_beta_ps_deg, -1.0);
                p.ps.handle_te_frac = std::max(0.005, (P - c.back()).dot(v) / span);
            } else if (i >= 2 && i <= n - 3) {
                // a free PS point sits at the law thickness below the SS point of its station:
                // dragging moves the station (keeping the point's x offset from it)
                auto& st = p.ps.mprime_ctrl;
                const auto j = static_cast<std::size_t>(i - 2);
                if (j < st.size()) {
                    const double off = c[static_cast<std::size_t>(i)].x() - st[j];
                    const double lo = j == 0 ? 0.02 : st[j - 1] + 0.01, hi = j + 1 == st.size() ? 0.98 : st[j + 1] - 0.01;
                    st[j] = std::clamp(P.x() - off, lo, hi);
                }
            }
        }
    } catch (const std::exception&) {
        return;  // an invalid intermediate profile: keep the last valid parameters
    }
    emit designEdited(done);
}

int ContentView::layerIndex() const
{
    const int n = static_cast<int>(s_.design.layers.size());
    return std::clamp(layer >= 0 ? layer : s_.active_layer, 0, std::max(0, n - 1));
}

QImage ContentView::image()
{
    if (stack_->currentWidget() == occt_) {
        waitForSolid();
        return occt_->image();
    }
    if (stack_->currentWidget() == plot_) return plot_->image();
    if (stack_->currentWidget() == view3d_) return view3d_->image();
    return table_->grab().toImage();
}

void ContentView::refresh(bool refit)
{
    const int k = layerIndex();
    QString t = contentName(content_);
    const bool per_layer = content_ == Content::BladeToBlade || content_ == Content::ProfileSection ||
                           content_ == Content::ProfileCurvature || content_ == Content::InfoCurrent ||
                           content_ == Content::AnglesCurrent || content_ == Content::ThicknessCurrent;
    if (content_ == Content::Solid3D) {
        t += annulus_ ? QString(" — whole annulus (%1 blades)").arg(s_.design.n_blades) : QString(" — single blade");
        t += exact_ ? QString(", exact OCCT solid") : QString(", design mesh");
    }
    if (per_layer && k < static_cast<int>(s_.design.layers.size()))
        t += QString(" — %1 (span %2)").arg(QString::fromStdString(s_.design.layers[static_cast<std::size_t>(k)].name)).arg(s_.design.layers[static_cast<std::size_t>(k)].span);
    if (content_ == Content::BladeSections)
        t += sections == SectionCoords::MRTheta ? " — M, R·Theta" : sections == SectionCoords::ZRTheta ? " — Z, R·Theta" : " — X, Y (axial view)";
    setWindowTitle(t);
    emit titleChanged(t);
    if (content_ == Content::View3D) {
        stack_->setCurrentWidget(view3d_);
        build3D(refit);
    } else if (content_ == Content::Solid3D) {
        stack_->setCurrentWidget(occt_);
        buildSolid();
    } else if (content_ == Content::InfoCurrent || content_ == Content::InfoAll) {
        stack_->setCurrentWidget(table_);
        buildTable();
    } else {
        stack_->setCurrentWidget(plot_);
        buildPlot(refit);
    }
}

// ─────────────────────────────────────────────────────────────────────────────

void ContentView::buildPlot(bool refit)
{
    QVector<Series> S;
    QVector<PointMarker> M;
    QVector<Fill> F;
    plot_->setColorBar(false);
    plot_->setBadge({});
    plot_->setHandles({});
    plot_->setTitle({});
    const auto R = s_.result;
    if (!R) {
        plot_->setSeries({});
        plot_->setMarkers({});
        plot_->setFills({});
        plot_->setBadge({"evaluating ...", 1});
        plot_->commit(refit);
        return;
    }
    const auto& D = s_.design;
    const std::size_t k = static_cast<std::size_t>(layerIndex());
    const bool have_k = k < R->layers.size() && R->layers[k].error.empty() && R->layers[k].section;
    if (!R->errors.empty()) {
        Badge b{"✗ design error: " + QString::fromStdString(R->errors.front()), 2};
        plot_->setBadge(b);
    } else if (k < R->layers.size() && !R->layers[k].error.empty()) {
        plot_->setBadge({"✗ layer error: " + QString::fromStdString(R->layers[k].error), 2});
    }
    auto hval = [&](const LayerResult& L, std::size_t i) {
        switch (haxis) {
        case HAxis::MeridionalFraction: return L.mfrac[i];
        case HAxis::MPrime: return L.mprime[i];
        case HAxis::Meridional: return L.m[i];
        case HAxis::ChordFraction: return L.xc[i];
        }
        return L.mfrac[i];
    };
    // horizontal-axis value at chord fraction x (interpolated along the camber points)
    auto hOfX = [&](const LayerResult& L, double x) {
        if (haxis == HAxis::ChordFraction || L.xc.size() < 2) return x;
        std::size_t i = 1;
        while (i + 1 < L.xc.size() && L.xc[i] < x) ++i;
        const double w = std::clamp((x - L.xc[i - 1]) / (L.xc[i] - L.xc[i - 1]), 0.0, 1.0);
        return hval(L, i - 1) + w * (hval(L, i) - hval(L, i - 1));
    };
    const QString hlabel = haxis == HAxis::MeridionalFraction ? "meridional chord fraction (LE 0 → TE 1)"
                           : haxis == HAxis::MPrime          ? "m' (from the layer LE)"
                           : haxis == HAxis::Meridional      ? "m [mm] (from the layer LE)"
                                                             : "chord fraction x/c";
    auto camberSeries = [&](std::size_t j, const std::vector<double>& y, const QString& name, Qt::PenStyle st = Qt::SolidLine) {
        const auto& L = R->layers[j];
        Series s;
        s.name = name;
        s.color = layerColor(j);
        s.style = st;
        for (std::size_t i = 0; i < y.size(); ++i) s.pts << QPointF(hval(L, i), y[i]);
        return s;
    };
    auto layerName = [&](std::size_t j) { return QString::fromStdString(D.layers[j].name); };

    switch (content_) {
    case Content::BladeToBlade: {
        if (!have_k) break;
        const auto& L = R->layers[k];
        const auto& loop = L.section->loop;
        const double pitch = 2.0 * Math::kPi / std::max(D.n_blades, 1);
        // (m', theta) on the layer's stream surface -> the chosen view coordinates
        auto map = [&](const Math::Vec2& mt, double dth) {
            const double th = mt.y() + dth;
            if (b2b == B2BCoords::MPrimeTheta) return QPointF(mt.x(), th);
            return QPointF(L.surface.m_of_mprime(mt.x()), L.surface.xr(mt.x()).y() * th);
        };
        for (double dth : {-pitch, pitch}) {
            Series s;
            s.name = "neighbour blades";
            s.color = QColor(170, 170, 170);
            s.legend = dth > 0;
            for (const auto& q : L.conformal) s.pts << map(q, dth);
            S << s;
        }
        std::vector<QPointF> pts;
        for (const auto& q : L.conformal) pts.push_back(map(q, 0.0));
        addPieces(S, pts, loop.piece, true);
        Series cam;
        cam.name = "camber line";
        cam.color = kCamber;
        cam.style = Qt::DashLine;
        for (std::size_t i = 0; i < L.mprime.size(); ++i) cam.pts << map({L.mprime_le + L.mprime[i], Math::deg2rad(L.theta_deg[i])}, 0.0);
        S << cam;
        const auto cf = confOf(L, D.mirror);
        for (const auto& is : L.section->quality.issues) {
            if (is.severity == Profile::Severity::Info) continue;
            const QPointF c = cf(is.where);
            M << issueMarker(map({c.x(), c.y()}, 0.0), is.severity);
        }
        M << PointMarker{map({L.mprime_le, L.theta_le}, 0.0), Qt::black, 'o', "LE"};
        M << PointMarker{map({L.mprime_le + L.mprime.back(), Math::deg2rad(L.theta_deg.back())}, 0.0), Qt::black, 'o', "TE"};
        M << PointMarker{map({L.cg_mprime, L.cg_theta}, 0.0), QColor(0, 0, 0), 's', "CG", 4.0, 1.6};
        if (b2b == B2BCoords::MPrimeTheta) plot_->setLabels("m'", "θ [rad]");
        else plot_->setLabels("m [mm]", "r θ [mm]");
        const QString rs = D.restagger.angle_deg != 0.0 ? QString(" (re-staggered %1%2°)").arg(D.restagger.angle_deg > 0 ? "+" : "").arg(D.restagger.angle_deg) : QString();
        plot_->setBadge(badgeFor(L.section->quality, QString("%1 · γ %2°%3 · %4")
                                                         .arg(QString::fromStdString(L.section->designation))
                                                         .arg(Math::rad2deg(L.gamma), 0, 'f', 3)
                                                         .arg(rs)
                                                         .arg(L.surface.cylindrical() ? "cylindrical wrap" : "surface of revolution")));
        break;
    }
    case Content::BladeSections: {
        // every layer's section in one plot, joined by the LE and TE lines; the axial view also shows
        // the hub and casing circles at the LE and TE radii (of the first and last layer)
        std::vector<QPointF> le, te, cgs;
        bool first = true;
        double tmin = 1e300, tmax = -1e300;
        for (std::size_t j = 0; j < R->layers.size(); ++j) {
            const auto& L = R->layers[j];
            if (!L.error.empty() || !L.section) continue;
            auto map = [&](const Math::Vec2& mt) {
                switch (sections) {
                case SectionCoords::MRTheta: return QPointF(L.surface.m_of_mprime(mt.x()), L.surface.xr(mt.x()).y() * mt.y());
                case SectionCoords::ZRTheta: {
                    const Math::Vec2 xr = L.surface.xr(mt.x());
                    return QPointF(xr.x(), xr.y() * mt.y());
                }
                case SectionCoords::XYAxial: {
                    const Math::Vec3 p = L.surface.cartesian(mt.x(), mt.y());  // (axial, r sin theta, r cos theta)
                    return QPointF(p.y(), p.z());
                }
                }
                return QPointF(mt.x(), mt.y());
            };
            std::vector<QPointF> pts;
            for (const auto& q : L.conformal) {
                pts.push_back(map(q));
                tmin = std::min(tmin, q.y());
                tmax = std::max(tmax, q.y());
            }
            addPieces(S, pts, L.section->loop.piece, true, 1.0, first);
            Series cam;
            cam.name = "camber line";
            cam.color = QColor(40, 40, 40);
            cam.width = 0.7;
            cam.legend = first;
            for (std::size_t i = 0; i < L.mprime.size(); ++i) cam.pts << map({L.mprime_le + L.mprime[i], Math::deg2rad(L.theta_deg[i])});
            S << cam;
            le.push_back(map({L.mprime_le, L.theta_le}));
            cgs.push_back(map({L.cg_mprime, L.cg_theta}));
            M << PointMarker{cgs.back(), QColor(0, 0, 0), 's', j == 0 ? QString("CG") : QString(), 4.0, 1.6};
            te.push_back(map({L.mprime_le + L.mprime.back(), Math::deg2rad(L.theta_deg.back())}));
            first = false;
        }
        for (const auto* line : {&le, &te}) {
            Series s;
            s.name = "LE / TE lines";
            s.color = kEdge;
            s.width = 1.3;
            s.legend = line == &le;
            for (const auto& p : *line) s.pts << p;
            S << s;
        }
        {
            Series s;
            s.name = "centres of gravity";
            s.color = QColor(0, 0, 0);
            s.style = Qt::DashLine;
            s.width = 1.2;
            for (const auto& p : cgs) s.pts << p;
            S << s;
        }
        if (R->blade_volume > 0.0) {
            const auto& g = R->blade_cg;
            plot_->setBadge({QString("blade centre of gravity (uniform density): Z %1 mm · R %2 mm · θ %3° · volume %4 mm³")
                                 .arg(g.x(), 0, 'f', 3)
                                 .arg(std::hypot(g.y(), g.z()), 0, 'f', 3)
                                 .arg(Math::rad2deg(std::atan2(g.y(), g.z())), 0, 'f', 3)
                                 .arg(R->blade_volume, 0, 'f', 1),
                             0});
        }
        if (sections == SectionCoords::XYAxial && !first) {
            const auto& h = R->layers.front();
            const auto& c = R->layers.back();
            const double pad = 0.25 * (tmax - tmin) + 0.02;
            bool leg = true;
            for (const double r : {h.r_le, h.r_te, c.r_le, c.r_te}) {
                Series s;
                s.name = "hub / casing circles at LE, TE";
                s.color = QColor(30, 30, 30);
                s.style = Qt::DashLine;
                s.legend = leg;
                leg = false;
                for (int i = 0; i <= 80; ++i) {
                    const double t = tmin - pad + (tmax - tmin + 2.0 * pad) * i / 80.0;
                    s.pts << QPointF(r * std::sin(t), r * std::cos(t));
                }
                S << s;
            }
        }
        switch (sections) {
        case SectionCoords::MRTheta: plot_->setLabels("M [mm] (from each layer's LE line)", "R·Theta [mm]"); break;
        case SectionCoords::ZRTheta: plot_->setLabels("Z (axial) [mm]", "R·Theta [mm]"); break;
        case SectionCoords::XYAxial: plot_->setLabels("X = r sin θ [mm]   (axial view from downstream)", "Y = r cos θ [mm]"); break;
        }
        break;
    }
    case Content::ProfileSection: {
        if (!have_k) break;
        const auto& L = R->layers[k];
        const auto& q = L.section->quality;
        for (const auto* sd : {&q.ss, &q.ps}) {
            std::vector<QPointF> pts;
            for (const auto& p : sd->xy) pts.emplace_back(p.x(), p.y());
            addPieces(S, pts, sd->piece, sd == &q.ss);
        }
        Series cam;
        cam.name = "camber line (inscribed circles)";
        cam.color = kCamber;
        cam.style = Qt::DashLine;
        for (const auto& p : L.section->camber.point) cam.pts << QPointF(p.x(), p.y());
        S << cam;
        for (const auto& is : q.issues) {
            if (is.severity == Profile::Severity::Info) continue;
            M << issueMarker(QPointF(is.where.x(), is.where.y()), is.severity);
        }
        if (L.section->control) {
            const auto& C = *L.section->control;
            auto polygon = [&](const std::vector<Math::Vec2>& c, const QString& name, QColor col, bool legend) {
                Series s;
                s.name = name;
                s.color = col;
                s.style = Qt::DotLine;
                s.legend = legend;
                for (const auto& v : c) s.pts << QPointF(v.x(), v.y());
                S << s;
            };
            polygon(C.ss, "control polygon (drag the dots)", QColor(25, 80, 225, 150), true);
            polygon(C.ps, "control polygon (drag the dots)", QColor(220, 30, 30, 150), false);
            polygon(C.le, "connector control points (fixed)", QColor(0, 150, 60, 150), true);
            polygon(C.te, "connector control points (fixed)", QColor(0, 150, 60, 150), false);
            QVector<DragHandle> H;
            for (std::size_t i = 0; i < C.ss.size(); ++i) {
                const bool free = i >= 1 && i + 1 < C.ss.size();
                H << DragHandle{kSsHandle + static_cast<int>(i), QPointF(C.ss[i].x(), C.ss[i].y()), kSS, !free, {}};
            }
            for (std::size_t i = 0; i < C.ps.size(); ++i) {
                const bool free = i >= 1 && i + 1 < C.ps.size();
                H << DragHandle{kPsHandle + static_cast<int>(i), QPointF(C.ps[i].x(), C.ps[i].y()), kPS, !free, {}};
            }
            for (const auto* c : {&C.le, &C.te})
                for (std::size_t i = 1; i + 1 < c->size(); ++i) H << DragHandle{-1, QPointF((*c)[i].x(), (*c)[i].y()), kEdge, true, {}};
            plot_->setHandles(H);
        }
        plot_->setLabels("x / c", "y / c");
        plot_->setBadge(badgeFor(q, QString::fromStdString(L.section->designation) +
                                        (L.section->control ? QString("  ·  drag a control point: SS handles / CP2, PS handles / stations") : QString())));
        break;
    }
    case Content::ProfileCurvature: {
        if (!have_k) break;
        const auto& q = R->layers[k].section->quality;
        const double wle = edgeWindow(q, true), wte = edgeWindow(q, false);
        auto inRange = [&](double f) {
            return crange == CurvatureRange::Full || (crange == CurvatureRange::LeadingEdge && f <= wle) ||
                   (crange == CurvatureRange::TrailingEdge && f >= 1.0 - wte);
        };
        for (const auto* sd : {&q.ss, &q.ps}) {
            std::vector<QPointF> pts;
            std::vector<int> pc;
            for (std::size_t i = 0; i < sd->xy.size(); ++i) {
                const double f = sd->s[i] / sd->s.back();
                if (!inRange(f)) continue;
                pts.emplace_back(f, sd->kappa[i]);
                pc.push_back(sd->piece[i]);
            }
            addPieces(S, pts, pc, sd == &q.ss);
        }
        for (const auto& is : q.issues) {
            if (is.severity == Profile::Severity::Info || is.side < 0 || !inRange(is.s)) continue;
            const auto& sd = is.side == 0 ? q.ss : q.ps;
            std::size_t best = 0;
            for (std::size_t i = 0; i < sd.s.size(); ++i)
                if (std::abs(sd.s[i] / sd.s.back() - is.s) < std::abs(sd.s[best] / sd.s.back() - is.s)) best = i;
            M << issueMarker(QPointF(is.s, sd.kappa[best]), is.severity);
        }
        plot_->setLabels("arc length / side length (LE point → TE point)", "κ [1/chord]  (symlog: linear below 1)");
        plot_->setBadge(badgeFor(q, crange == CurvatureRange::Full ? "curvature = 1 / R, signed (> 0 convex)"
                                    : crange == CurvatureRange::LeadingEdge ? "leading-edge region" : "trailing-edge region"));
        break;
    }
    case Content::MeridionalContour:
    case Content::MeridionalBeta:
    case Content::MeridionalThickness: {
        if (content_ != Content::MeridionalContour) {
            const auto& fld = content_ == Content::MeridionalBeta ? R->beta : R->thickness;
            if (fld.ni > 1 && fld.nj > 1) {
                const int bands = 16;
                const double lo = fld.lo, hi = fld.hi > fld.lo ? fld.hi : fld.lo + 1.0;
                auto cmap = [&](double t) { return turbo(std::clamp((std::floor(t * bands) + 0.5) / bands, 0.0, 1.0)); };
                for (int j = 0; j + 1 < fld.nj; ++j)
                    for (int i = 0; i + 1 < fld.ni; ++i) {
                        auto id = [&](int a, int b) { return static_cast<std::size_t>(a + fld.ni * b); };
                        const double v = 0.25 * (fld.value[id(i, j)] + fld.value[id(i + 1, j)] + fld.value[id(i + 1, j + 1)] + fld.value[id(i, j + 1)]);
                        QPolygonF poly;
                        for (const auto& [a, b] : {std::pair{i, j}, {i + 1, j}, {i + 1, j + 1}, {i, j + 1}}) poly << QPointF(fld.xr[id(a, b)].x(), fld.xr[id(a, b)].y());
                        F << Fill{poly, cmap((v - lo) / (hi - lo))};
                    }
                plot_->setColorBar(true, lo, hi, content_ == Content::MeridionalBeta ? "metal angle β [deg]" : "thickness t / c", cmap);
            }
        } else if (R->ok()) {
            QPolygonF blade;
            const auto& h = R->layers.front(), &c = R->layers.back();
            for (const auto& p : h.camber_xr) blade << QPointF(p.x(), p.y());
            for (auto it = c.camber_xr.rbegin(); it != c.camber_xr.rend(); ++it) blade << QPointF(it->x(), it->y());
            F << Fill{blade, QColor(225, 232, 250)};
        }
        auto poly = [&](const std::vector<Math::Vec2>& p, const QString& name, QColor col, double w, Qt::PenStyle st, bool legend) {
            Series s;
            s.name = name;
            s.color = col;
            s.width = w;
            s.style = st;
            s.legend = legend;
            for (const auto& q : p) s.pts << QPointF(q.x(), q.y());
            S << s;
        };
        poly(R->hub, "hub", Qt::black, 1.4, Qt::SolidLine, true);
        poly(R->casing, "casing", QColor(60, 60, 60), 1.4, Qt::SolidLine, true);
        for (std::size_t j = 1; j + 1 < R->streamlines.size(); ++j)
            poly(R->streamlines[j], layerName(j) + " layer", layerColor(j), 1.0, Qt::DashLine, content_ == Content::MeridionalContour);
        poly(R->le_line, "LE / TE lines", QColor(0, 0, 200), 1.2, Qt::SolidLine, true);
        poly(R->te_line, "LE / TE lines", QColor(0, 0, 200), 1.2, Qt::SolidLine, false);
        const auto& m = D.meridional;
        for (const auto* w : {&m.hub, &m.casing}) {
            M << PointMarker{QPointF(w->z_le, w->r_le), Qt::black, 's', QString("Z %1  R %2").arg(w->z_le).arg(w->r_le)};
            M << PointMarker{QPointF(w->z_te, w->r_te), Qt::black, 's', QString("Z %1  R %2").arg(w->z_te).arg(w->r_te)};
        }
        plot_->setLabels("Z (axial) [mm]", "R [mm]");
        if (content_ == Content::MeridionalContour && R->errors.empty())
            plot_->setBadge({QString("hub: %1 · casing: %2").arg(m.hub_law == Blade::WallLaw::BellMehta ? "Bell & Mehta" : "linear",
                                                               m.casing_law == Blade::WallLaw::BellMehta ? "Bell & Mehta" : "linear"), 0});
        break;
    }
    case Content::MeridionalCurvature: {
        if (R->layers.empty()) break;
        for (std::size_t j = 0; j < R->layers.size(); ++j) {
            const auto& L = R->layers[j];
            if (!L.error.empty()) continue;
            Series s;
            s.name = layerName(j);
            s.color = layerColor(j);
            const auto& Sf = L.surface;
            for (int i = 0; i <= 240; ++i) {
                const double mp = Sf.mprime_inlet() + (Sf.mprime_outlet() - Sf.mprime_inlet()) * i / 240.0;
                s.pts << QPointF(Sf.xr(mp).x(), Sf.curvature(mp));
            }
            S << s;
        }
        plot_->setLabels("Z (axial) [mm]", "streamline curvature in (z, R) [1/mm]");
        break;
    }
    case Content::LeTeTheta:
    case Content::LeTeBeta:
    case Content::LeTeParameters:
    case Content::AirfoilArea: {
        Series a, b, c;
        a.markers = b.markers = c.markers = true;
        a.color = QColor(30, 60, 235);
        b.color = QColor(220, 30, 30);
        c.color = QColor(0, 150, 0);
        for (std::size_t j = 0; j < R->layers.size(); ++j) {
            const auto& L = R->layers[j];
            if (!L.error.empty() || !L.section) continue;
            const double sp = D.layers[j].span;
            if (content_ == Content::LeTeTheta) {
                a.pts << QPointF(Math::rad2deg(L.theta_le), sp);
                b.pts << QPointF(L.theta_te_deg, sp);
            } else if (content_ == Content::LeTeBeta) {
                a.pts << QPointF(L.beta_le, sp);
                b.pts << QPointF(L.beta_te, sp);
            } else if (content_ == Content::LeTeParameters) {
                a.pts << QPointF(L.section->le_radius, sp);
                b.pts << QPointF(L.section->te_radius, sp);
                c.pts << QPointF(L.section->camber.max_thickness, sp);
            } else {
                a.pts << QPointF(L.section->area * L.chord * L.chord, sp);
            }
        }
        if (content_ == Content::LeTeTheta) {
            a.name = "θ LE";
            b.name = "θ TE";
            plot_->setLabels("θ [deg]", "span fraction");
        } else if (content_ == Content::LeTeBeta) {
            a.name = "β LE (metal)";
            b.name = "β TE (metal)";
            plot_->setLabels("β [deg]", "span fraction");
        } else if (content_ == Content::LeTeParameters) {
            a.name = "LE radius / c";
            b.name = "TE radius / c";
            c.name = "t_max / c";
            plot_->setLabels("fraction of chord", "span fraction");
        } else {
            a.name = "airfoil area";
            plot_->setLabels("section area [mm²]", "span fraction");
        }
        S << a;
        if (!b.pts.isEmpty()) S << b;
        if (!c.pts.isEmpty()) S << c;
        break;
    }
    case Content::AnglesCurrent: {
        if (!have_k) break;
        auto sb = camberSeries(k, R->layers[k].beta_deg, "β (metal angle)");
        auto st = camberSeries(k, R->layers[k].theta_deg, "θ (camber)", Qt::DashLine);
        sb.color = QColor(30, 60, 235);
        st.color = QColor(220, 30, 30);
        S << sb << st;
        plot_->setLabels(hlabel, "angle [deg]");
        break;
    }
    case Content::BetaAll:
    case Content::ThetaAll:
    case Content::ThicknessAll:
        for (std::size_t j = 0; j < R->layers.size(); ++j) {
            if (!R->layers[j].error.empty()) continue;
            const auto& L = R->layers[j];
            const auto& y = content_ == Content::BetaAll ? L.beta_deg : content_ == Content::ThetaAll ? L.theta_deg : L.thickness;
            const QString what = content_ == Content::BetaAll ? "Beta" : content_ == Content::ThetaAll ? "Theta" : "t/c";
            S << camberSeries(j, y, QString("%1 @ Span: %2").arg(what).arg(D.layers[j].span, 0, 'f', 4));
        }
        plot_->setLabels(hlabel, content_ == Content::BetaAll ? "Beta [deg]" : content_ == Content::ThetaAll ? "Theta [deg]" : "thickness t / c");
        if (content_ == Content::ThicknessAll && have_k && R->layers[k].section->thickness_law)
            plot_->setBadge({QString::fromStdString(D.layers[k].name) + ": " + thicknessFormula(D.layers[k].custom, R->layers[k].section->xi_max), 0});
        break;
    case Content::BetaVsTheta:
        for (std::size_t j = 0; j < R->layers.size(); ++j) {
            const auto& L = R->layers[j];
            if (!L.error.empty()) continue;
            Series s;
            s.name = QString("Beta @ Span: %1").arg(D.layers[j].span, 0, 'f', 4);
            s.color = layerColor(j);
            for (std::size_t i = 0; i < L.beta_deg.size(); ++i) s.pts << QPointF(L.theta_deg[i], L.beta_deg[i]);
            S << s;
        }
        plot_->setLabels("Theta [deg]", "Beta [deg]");
        break;
    case Content::ThicknessCurrent: {
        if (!have_k) break;
        auto s = camberSeries(k, R->layers[k].thickness, "t / c (inscribed circle)");
        s.color = QColor(30, 60, 235);
        S << s;
        const auto& sec = *R->layers[k].section;
        if (sec.thickness_law) {
            Series law;
            law.name = "law t(ξ), normal to the SS, at the SS point's x";
            law.color = QColor(220, 30, 30);
            law.style = Qt::DashLine;
            const auto& L = R->layers[k];
            for (std::size_t i = 0; i < sec.thickness_law->law.size(); ++i)
                law.pts << QPointF(hOfX(L, sec.thickness_law->mprime[i]), sec.thickness_law->law[i]);
            S << law;
            plot_->setBadge({thicknessFormula(D.layers[k].custom, sec.xi_max), 0});
        }
        plot_->setLabels(hlabel, "thickness t / c");
        break;
    }
    default: break;
    }
    plot_->setSeries(S);
    plot_->setMarkers(M);
    plot_->setFills(F);
    plot_->commit(refit);
}

void ContentView::build3D(bool refit)
{
    const auto R = s_.result;
    if (!R || R->blade.p.empty()) {
        view3d_->setScene({}, {}, false);
        return;
    }
    auto grid = [](const Design::Grid& g, const QString& name, QColor c, bool surf, bool pts) {
        Grid3D o;
        o.name = name;
        o.ni = g.ni;
        o.nj = g.nj;
        o.color = c;
        o.surface = surf;
        o.points = pts;
        for (const auto& p : g.p) o.p.push_back({p.x(), p.y(), p.z()});
        return o;
    };
    std::vector<Grid3D> G;
    G.push_back(grid(R->blade, "blade", QColor(170, 20, 20), true, false));
    if (show_walls_) {
        G.push_back(grid(R->hub_wall, "hub", QColor(0, 150, 0), false, true));
        G.push_back(grid(R->casing_wall, "casing", QColor(40, 40, 40), false, true));
    }
    std::vector<Line3D> L;
    // tip cap: chordwise ribs between suction and pressure side points of the last row
    const int j = R->blade.nj - 1, n = R->blade.ni - 1;
    for (int i = 1; i < n / 2; i += 2) {
        const auto& a = R->blade.at(i, j);
        const auto& b = R->blade.at(n - i, j);
        L.push_back({{{a.x(), a.y(), a.z()}, {b.x(), b.y(), b.z()}}, QColor(200, 0, 200), Qt::SolidLine, 1.0});
    }
    Line3D tip;
    tip.color = QColor(200, 0, 200);
    tip.width = 1.5;
    for (int i = 0; i <= n; ++i) tip.p.push_back({R->blade.at(i, j).x(), R->blade.at(i, j).y(), R->blade.at(i, j).z()});
    L.push_back(tip);
    // centres of gravity: the sections' (dashed line) and the blade's (cross)
    const QColor cgc(0, 0, 0);
    Line3D cgl;
    cgl.color = cgc;
    cgl.style = Qt::DashLine;
    cgl.width = 1.5;
    for (const auto& lay : R->layers)
        if (lay.error.empty()) cgl.p.push_back({lay.cg_xyz.x(), lay.cg_xyz.y(), lay.cg_xyz.z()});
    if (cgl.p.size() > 1) L.push_back(cgl);
    if (R->blade_volume > 0.0) {
        const auto& g = R->blade_cg;
        const double h = 2.0;  // mm
        for (const Math::Vec3& e : {Math::Vec3(h, 0, 0), Math::Vec3(0, h, 0), Math::Vec3(0, 0, h)})
            L.push_back({{{g.x() - e.x(), g.y() - e.y(), g.z() - e.z()}, {g.x() + e.x(), g.y() + e.y(), g.z() + e.z()}}, cgc, Qt::SolidLine, 2.5});
    }
    // rotation axis (dash-dot)
    double x0 = 1e300, x1 = -1e300;
    for (const auto& p : R->hub_wall.p) {
        x0 = std::min(x0, p.x());
        x1 = std::max(x1, p.x());
    }
    L.push_back({{{x0, 0, 0}, {x1, 0, 0}}, QColor(40, 40, 220), Qt::DashDotLine, 1.0});
    view3d_->setScene(std::move(G), std::move(L), refit);
}

// ── 3D Solid ─────────────────────────────────────────────────────────────────

void ContentView::setAnnulus(bool on)
{
    if (on == annulus_) return;
    annulus_ = on;
    solid_refit_ = true;
    refresh();
}

void ContentView::setExact(bool on)
{
    if (on == exact_) return;
    exact_ = on;
    solid_for_.reset();  // the exact solid is built anew for the current result
    refresh();
}

void ContentView::buildSolid()
{
    const auto R = s_.result;
    if (!R) {
        occt_->setMessage("evaluating ...");
        return;
    }
    if (!R->ok()) {
        occt_->setMessage("✗ the design has errors: no solid (see Geometry checks)");
        return;
    }
    if (!exact_) {  // the design rows as triangles: immediate
        try {
            const auto t0 = std::chrono::steady_clock::now();
            auto m = std::make_shared<const Design::SceneMesh>(Design::build_scene_mesh(s_.design, *R, annulus_));
            occt_->setMeshScene(std::move(m), !solid_shown_ || solid_refit_);  // includes the OpenGL arrays
            solid_shown_ = true;
            solid_refit_ = false;
            solid_error_.clear();
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            occt_->setMessage(QString("%1 · design mesh · updated in %2 ms")
                                  .arg(annulus_ ? QString("%1 blades").arg(s_.design.n_blades) : QString("1 of %1 blades").arg(s_.design.n_blades))
                                  .arg(ms, 0, 'f', 0));
        } catch (const std::exception& e) {
            solid_error_ = QString::fromUtf8(e.what());
            occt_->setMessage("✗ no display mesh: " + solid_error_);
        }
        return;
    }
    if (R == solid_for_ && annulus_ == solid_annulus_) return;  // built (or being built) for this result
    if (solid_watch_->isRunning()) {
        solid_pending_ = true;
        return;
    }
    startSolid();
}

void ContentView::startSolid()
{
    solid_for_ = s_.result;
    solid_annulus_ = annulus_;
    solid_pending_ = false;
    solid_applied_ = false;
    const auto d = s_.design;
    const bool ann = annulus_;
    occt_->setMessage(ann ? QString("building the solid (%1 blades) ...").arg(d.n_blades) : QString("building the solid ..."));
    auto* ev = solid_ev_.get();
    solid_watch_->setFuture(QtConcurrent::run([d, ann, ev] {
        SolidJob j;
        j.annulus = ann;
        const auto t0 = std::chrono::steady_clock::now();
        try {
            auto sc = std::make_shared<Design::Scene>(Design::build_scene(d, *ev, {ann, 11}));
            OcctView::prepare(*sc);
            j.scene = std::move(sc);
        } catch (const Standard_Failure& e) {
            j.error = QString("OCCT: %1").arg(e.GetMessageString());
        } catch (const std::exception& e) {
            j.error = QString::fromUtf8(e.what());
        }
        j.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return j;
    }));
}

void ContentView::solidDone()
{
    if (solid_applied_ || !solid_watch_->future().isFinished()) return;
    solid_applied_ = true;
    if (!exact_) return;  // switched back to the display mesh meanwhile
    const SolidJob j = solid_watch_->result();
    if (!j.error.isEmpty()) {
        solid_error_ = j.error;
        occt_->setMessage("✗ solid not built: " + j.error);
    } else {
        solid_error_.clear();
        occt_->setScene(j.scene, !solid_shown_ || solid_refit_);
        solid_shown_ = true;
        solid_refit_ = false;
        occt_->setMessage(QString("%1 · exact OCCT solid · built in %2 s")
                              .arg(j.annulus ? QString("%1 blades").arg(j.scene->copies) : QString("1 of %1 blades").arg(j.scene->n_blades))
                              .arg(j.seconds, 0, 'f', 2));
    }
    // a newer result (or a blade-count switch) arrived meanwhile
    if (exact_ && content_ == Content::Solid3D && (solid_pending_ || s_.result != solid_for_ || annulus_ != solid_annulus_)) buildSolid();
}

QString ContentView::waitForSolid()
{
    if (content_ != Content::Solid3D) return "not a solid view";
    if (!exact_) return solid_shown_ ? solid_error_ : QString("no display mesh");
    if (solid_watch_->isRunning()) solid_watch_->waitForFinished();
    solidDone();
    if (solid_watch_->isRunning()) {  // solidDone started the newer result: wait for it too
        solid_watch_->waitForFinished();
        solidDone();
    }
    if (!solid_shown_ && solid_error_.isEmpty()) return "no solid built";
    return solid_error_;
}

void ContentView::buildTable()
{
    table_->clear();
    const auto R = s_.result;
    if (!R) return;
    const auto& D = s_.design;
    QStringList props = {"construction", "span fraction", "stagger [deg]", "β LE [deg]", "β TE [deg]", "camber β LE − β TE [deg]",
                         "θ LE [deg]", "θ TE [deg]", "mapping", "R LE [mm]", "R TE [mm]", "conformal chord c'", "meridional length LE→TE [mm]",
                         "chord [mm]", "pitch at LE [mm]", "solidity (3-D, n_blades)", "solidity (2-D input)", "t_max / c", "x/c of t_max",
                         "LE radius / c", "TE radius / c", "area [mm²]", "centre of gravity Z [mm]", "centre of gravity R [mm]",
                         "centre of gravity θ [deg]", "quality", "fit RMS / deviation [% c]",
                         "NACA 65 equivalent arc: γ + χ LE / γ − χ TE [deg]"};
    auto column = [&](std::size_t j) {
        QStringList v;
        const auto& L = R->layers[j];
        const auto& l = D.layers[j];
        if (!L.error.empty() || !L.section) {
            v << "error: " + QString::fromStdString(L.error);
            while (v.size() < props.size()) v << "";
            return v;
        }
        const auto& sec = *L.section;
        const int ne = sec.quality.count(Profile::Severity::Error), nw = sec.quality.count(Profile::Severity::Warning);
        v << QString::fromStdString(sec.designation) << num(l.span) << num(l.stagger_deg) << num(L.beta_le, 5) << num(L.beta_te, 5)
          << num(L.beta_le - L.beta_te, 5) << num(Math::rad2deg(L.theta_le), 5) << num(L.theta_te_deg, 5)
          << (L.surface.cylindrical() ? "cylindrical wrap" : "surface of revolution") << num(L.r_le, 6) << num(L.r_te, 6)
          << num(L.chord_mprime, 6) << num(L.surface.m_te(), 6) << num(L.chord, 6) << num(L.pitch_le, 6) << num(L.solidity_3d, 4)
          << num(l.solidity, 4) << num(sec.camber.max_thickness, 4) << num(sec.camber.max_thickness_x, 4) << num(sec.le_radius, 4)
          << num(sec.te_radius, 4) << num(L.section_area, 5) << num(L.cg_xr.x(), 6) << num(L.cg_xr.y(), 6)
          << num(Math::rad2deg(L.cg_theta), 6)
          << (ne ? QString("%1 error(s)").arg(ne) : nw ? QString("%1 warning(s)").arg(nw) : QString("OK"))
          << (sec.deviation_rms ? num(100.0 * *sec.deviation_rms, 3) : l.fitted_to ? num(100.0 * l.fit_rms, 3) : QString("—"))
          << (sec.equivalent_arc_deg ? QString("%1 / %2").arg(num(l.stagger_deg + (*sec.equivalent_arc_deg)[0], 5), num(l.stagger_deg - (*sec.equivalent_arc_deg)[1], 5))
                                     : QString("—"));
        return v;
    };
    if (content_ == Content::InfoCurrent) {
        const auto k = static_cast<std::size_t>(layerIndex());
        if (k >= R->layers.size()) return;
        table_->setColumnCount(2);
        table_->setRowCount(static_cast<int>(props.size()));
        table_->setHorizontalHeaderLabels({"property", QString::fromStdString(D.layers[k].name)});
        const auto v = column(k);
        for (int i = 0; i < props.size(); ++i) {
            table_->setItem(i, 0, new QTableWidgetItem(props[i]));
            table_->setItem(i, 1, new QTableWidgetItem(v.value(i)));
        }
    } else {
        table_->setColumnCount(static_cast<int>(R->layers.size()) + 1);
        table_->setRowCount(static_cast<int>(props.size()));
        QStringList h = {"property"};
        for (const auto& l : D.layers) h << QString::fromStdString(l.name);
        table_->setHorizontalHeaderLabels(h);
        for (int i = 0; i < props.size(); ++i) table_->setItem(i, 0, new QTableWidgetItem(props[i]));
        for (std::size_t j = 0; j < R->layers.size() && j < D.layers.size(); ++j) {
            const auto v = column(j);
            for (int i = 0; i < props.size(); ++i) table_->setItem(i, static_cast<int>(j) + 1, new QTableWidgetItem(v.value(i)));
        }
    }
    table_->resizeColumnsToContents();
}

// ─────────────────────────────────────────────────────────────────────────────

void ContentView::saveData()
{
    if (content_ == Content::View3D) {
        const QString f = QFileDialog::getSaveFileName(this, "Save mesh", "blade_mesh.obj", "Wavefront OBJ (*.obj)");
        if (f.isEmpty() || !s_.result) return;
        try {
            Design::write_obj(*s_.result, f.toStdString());
        } catch (const std::exception& e) {
            QMessageBox::warning(this, "Save mesh", e.what());
        }
        return;
    }
    const QString f = QFileDialog::getSaveFileName(this, "Save graph data", "graph.csv", "CSV (*.csv)");
    if (f.isEmpty()) return;
    QFile file(f);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
    QTextStream out(&file);
    if (stack_->currentWidget() == table_) {
        for (int r = -1; r < table_->rowCount(); ++r) {
            QStringList row;
            for (int c = 0; c < table_->columnCount(); ++c)
                row << (r < 0 ? table_->horizontalHeaderItem(c)->text() : table_->item(r, c) ? table_->item(r, c)->text() : QString());
            out << row.join(",") << "\n";
        }
    } else {
        out << plot_->csv();
    }
}

void ContentView::showMenu(const QPoint& global)
{
    QMenu m(this);
    // layer
    auto* lm = m.addMenu("Layer");
    auto* lg = new QActionGroup(&m);
    auto* follow = lm->addAction("Follow the active layer");
    follow->setCheckable(true);
    follow->setChecked(layer < 0);
    lg->addAction(follow);
    connect(follow, &QAction::triggered, this, [this] { layer = -1; refresh(); });
    lm->addSeparator();
    for (int j = 0; j < static_cast<int>(s_.design.layers.size()); ++j) {
        const auto& l = s_.design.layers[static_cast<std::size_t>(j)];
        auto* a = lm->addAction(QString("%1 (span %2)").arg(QString::fromStdString(l.name)).arg(l.span));
        a->setCheckable(true);
        a->setChecked(layer == j);
        lg->addAction(a);
        connect(a, &QAction::triggered, this, [this, j] { layer = j; refresh(); });
    }
    m.addSeparator();
    auto* nv = m.addMenu("New Auxiliary View");
    auto* cm = m.addMenu("Auxiliary View Content");
    for (const auto& grp : contentGroups()) {
        for (Content c : grp) {
            nv->addAction(contentName(c), this, [this, c] { emit newViewWanted(c); });
            auto* a = cm->addAction(contentName(c), this, [this, c] { setContent(c); });
            a->setCheckable(true);
            a->setChecked(c == content_);
        }
        nv->addSeparator();
        cm->addSeparator();
    }
    auto* ha = m.addMenu("Ang/Thk Plot Horizontal Axis");
    auto* hg = new QActionGroup(&m);
    const std::pair<HAxis, const char*> axes[] = {{HAxis::MeridionalFraction, "Meridional chord fraction"},
                                                  {HAxis::MPrime, "M-prime (m')"},
                                                  {HAxis::Meridional, "Meridional length m"},
                                                  {HAxis::ChordFraction, "Chord fraction x/c"}};
    for (const auto& [h, name] : axes) {
        auto* a = ha->addAction(name, this, [this, hv = h] { haxis = hv; refresh(true); });
        a->setCheckable(true);
        a->setChecked(haxis == h);
        hg->addAction(a);
    }
    auto* bc = m.addMenu("Blade-to-Blade Coordinates");
    auto* bg = new QActionGroup(&m);
    for (const auto& [c, name] : {std::pair{B2BCoords::MPrimeTheta, "M-prime, Theta (conformal)"}, std::pair{B2BCoords::MRTheta, "M, R·Theta [mm]"}}) {
        auto* a = bc->addAction(name, this, [this, cv = c] { b2b = cv; refresh(true); });
        a->setCheckable(true);
        a->setChecked(b2b == c);
        bg->addAction(a);
    }
    auto* sc = m.addMenu("Blade Sections Coordinates");
    auto* sg = new QActionGroup(&m);
    for (const auto& [c, name] : {std::pair{SectionCoords::MRTheta, "M, RTheta Plot"}, std::pair{SectionCoords::ZRTheta, "Z, RTheta Plot"},
                                  std::pair{SectionCoords::XYAxial, "X, Y Plot (axial view)"}}) {
        auto* a = sc->addAction(name, this, [this, cv = c] { sections = cv; refresh(true); });
        a->setCheckable(true);
        a->setChecked(sections == c);
        sg->addAction(a);
    }
    auto* cr = m.addMenu("Curvature Range");
    auto* cg = new QActionGroup(&m);
    for (const auto& [c, name] : {std::pair{CurvatureRange::Full, "LE point → TE point"}, std::pair{CurvatureRange::LeadingEdge, "Leading-edge zoom"},
                                  std::pair{CurvatureRange::TrailingEdge, "Trailing-edge zoom"}}) {
        auto* a = cr->addAction(name, this, [this, cv = c] { crange = cv; refresh(true); });
        a->setCheckable(true);
        a->setChecked(crange == c);
        cg->addAction(a);
    }
    m.addSeparator();
    auto toggle = [&](const QString& name, bool& flag, bool refit, bool enabled = true) {
        auto* a = m.addAction(name);
        a->setCheckable(true);
        a->setChecked(flag);
        a->setEnabled(enabled);
        connect(a, &QAction::toggled, this, [this, &flag, refit](bool on) { flag = on; refresh(refit); });
    };
    const bool is_plot = stack_->currentWidget() == plot_, is3d = content_ == Content::View3D;
    toggle("Show Legend", plot_->showLegend, false, is_plot);
    toggle("Show Status Box", plot_->showBadge, false, is_plot);
    toggle("Show Points", is3d ? view3d_->showPoints : plot_->showPoints, false, is_plot || is3d);
    if (is3d) {
        toggle("Shaded", view3d_->shaded, false);
        toggle("Show Hub and Casing", show_walls_, false);
    }
    if (content_ == Content::Solid3D) {
        auto* gm = m.addMenu("Geometry");
        auto* ggp = new QActionGroup(&m);
        for (const bool ex : {false, true}) {
            auto* a = gm->addAction(ex ? QString("Exact OCCT solid (trimmed B-spline solid, slow: ~2 s per change)")
                                       : QString("Design mesh (fast: follows every change)"),
                                    this, [this, ex] { setExact(ex); });
            a->setCheckable(true);
            a->setChecked(exact_ == ex);
            ggp->addAction(a);
        }
        auto* bm = m.addMenu("Blades");
        auto* bgp = new QActionGroup(&m);
        for (const bool ann : {false, true}) {
            auto* a = bm->addAction(ann ? QString("Whole annulus (%1 blades)").arg(s_.design.n_blades) : QString("Single blade (one pitch of hub and casing)"),
                                    this, [this, ann] { setAnnulus(ann); });
            a->setCheckable(true);
            a->setChecked(annulus_ == ann);
            bgp->addAction(a);
        }
        const std::pair<OcctView::Part, const char*> parts[] = {
            {OcctView::Part::Blade, "Blade Display"}, {OcctView::Part::Hub, "Hub Display"}, {OcctView::Part::Casing, "Casing Display"}};
        const std::pair<OcctView::Look, const char*> looks[] = {{OcctView::Look::Hidden, "Hidden"},
                                                                {OcctView::Look::Wireframe, "Wireframe"},
                                                                {OcctView::Look::Solid, "Solid (opaque, with edges)"},
                                                                {OcctView::Look::Shaded, "Shaded (smooth, see-through)"}};
        for (const auto& [part, name] : parts) {
            auto* pm = m.addMenu(name);
            auto* pg = new QActionGroup(&m);
            for (const auto& [lk, lname] : looks) {
                if (part == OcctView::Part::Blade && lk == OcctView::Look::Hidden) continue;
                auto* a = pm->addAction(lname, this, [this, p = part, l = lk] { occt_->setLook(p, l); });
                a->setCheckable(true);
                a->setChecked(occt_->look(part) == lk);
                pg->addAction(a);
            }
        }
        auto* vm = m.addMenu("Standard View");
        const std::pair<OcctView::Direction, const char*> dirs[] = {{OcctView::Direction::Iso, "Isometric"},
                                                                    {OcctView::Direction::Axial, "Axial (along the machine axis)"},
                                                                    {OcctView::Direction::Meridional, "Meridional (from the side)"},
                                                                    {OcctView::Direction::Radial, "Radial (onto the blade-to-blade plane)"}};
        for (const auto& [dir, name] : dirs) vm->addAction(name, this, [this, dv = dir] { occt_->setDirection(dv); });
    }
    m.addSeparator();
    toggle("Swap Axes", plot_->swapAxes, true, is_plot);
    toggle("Proportional Axes", plot_->proportional, true, is_plot);
    m.addSeparator();
    toggle("Flip X Axis Direction", plot_->flipX, false, is_plot);
    toggle("Flip Y Axis Direction", plot_->flipY, false, is_plot);
    toggle("Symmetric Log Y Scale", plot_->symlogY, true, is_plot);
    m.addSeparator();
    m.addAction("Number of Points...", this, [this] { emit meshPointsWanted(); });
    m.addAction("Save Graph Data...", this, [this] { saveData(); });
    m.addAction("Copy View Image to Clipboard", this, [this] { QApplication::clipboard()->setImage(image()); });
    m.addSeparator();
    m.addAction("Zoom Fit", this, [this] {
        if (content_ == Content::View3D) view3d_->zoomFit();
        else if (content_ == Content::Solid3D) occt_->fitAll();
        else plot_->zoomFit();
    });
    m.exec(global);
}

}  // namespace Designer
