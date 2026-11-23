// windows.h first: OCCT's Standard_Macro.hxx defines NOUSER & co. for any windows.h after it,
// and the hidden context window needs user32 (DefWindowProcW, CS_OWNDC, WS_POPUP)
#if defined(_WIN32) && defined(BLADEGEN_HAVE_OCCT_VIEWER)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#include "OcctView.hpp"

#include <QContextMenuEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include <array>
#include <cmath>
#include <cstring>

#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
#include <AIS_InteractiveContext.hxx>
#include <AIS_InteractiveObject.hxx>
#include <AIS_Shape.hxx>
#include <AIS_ViewController.hxx>
#include <Aspect_DisplayConnection.hxx>
#include <Aspect_NeutralWindow.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Graphic3d_ArrayOfPolylines.hxx>
#include <Graphic3d_ArrayOfTriangles.hxx>
#include <Graphic3d_AspectFillArea3d.hxx>
#include <Graphic3d_AspectLine3d.hxx>
#include <Prs3d_Presentation.hxx>
#include <PrsMgr_PresentationManager.hxx>
#include <SelectMgr_Selection.hxx>
#include <Image_PixMap.hxx>
#include <OpenGl_Context.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <Prs3d_Drawer.hxx>
#include <Prs3d_IsoAspect.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Standard_Failure.hxx>
#include <V3d_Trihedron.hxx>
#include <V3d_View.hxx>
#include <V3d_Viewer.hxx>
#if defined(_WIN32)
#include <WNT_WClass.hxx>
#include <WNT_Window.hxx>
#elif defined(__APPLE__)
#include <Cocoa_Window.hxx>
#else
#include <Xw_Window.hxx>
#endif
#endif

namespace Designer {

namespace {

#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
// Tessellation of the shown shapes (metres): the drawers ask for the same, so AIS reuses it.
// 0.05 mm on the blade (about 1/1000 of the chord) and 0.2 mm on the walls look smooth on screen;
// the angular limit (about 11 deg) keeps the LE / TE rounded.
constexpr double kBladeDeflection = 5e-5, kBladeAngle = 0.2;
constexpr double kWallDeflection = 2e-4, kWallAngle = 0.2;
#endif

}  // namespace

#if defined(BLADEGEN_HAVE_OCCT_VIEWER)

namespace {

/// One part of the fast display (blade with its copies, hub or casing) as OpenGL primitive
/// arrays: mode 1 shaded triangles with vertex normals (+ outline polylines for the solid look),
/// mode 0 wireframe polylines.  Not selectable.
class MeshObject : public AIS_InteractiveObject
{
    DEFINE_STANDARD_RTTI_INLINE(MeshObject, AIS_InteractiveObject)
public:
    MeshObject(std::shared_ptr<const BladeGen::Design::SceneMesh> scene, std::size_t part, const Quantity_Color& color,
               Graphic3d_NameOfMaterial material)
        : scene_(std::move(scene)), part_(part), color_(color), material_(material)
    {
        setLook(false, 0.0);
    }

    void setLook(bool edges, double transparency)
    {
        edges_ = edges;
        Graphic3d_MaterialAspect mat(material_);
        mat.SetColor(color_);
        mat.SetTransparency(static_cast<float>(transparency));
        fill_ = new Graphic3d_AspectFillArea3d(Aspect_IS_SOLID, color_, Quantity_NOC_BLACK, Aspect_TOL_SOLID, 1.0, mat, mat);
        fill_->SetDistinguishOff();
        Quantity_Color dark = color_;
        dark.ChangeIntensity(-40.0);
        wire_ = new Graphic3d_AspectLine3d(dark, Aspect_TOL_SOLID, 1.0);
        edge_ = new Graphic3d_AspectLine3d(Quantity_Color(0.15, 0.15, 0.15, Quantity_TOC_sRGB), Aspect_TOL_SOLID, 1.0);
    }

    Standard_Boolean AcceptDisplayMode(const Standard_Integer mode) const override { return mode == 0 || mode == 1; }

protected:
    void Compute(const Handle(PrsMgr_PresentationManager)&, const Handle(Prs3d_Presentation)& prs, const Standard_Integer mode) override
    {
        const auto& sc = *scene_;
        const auto& m = part_ == 0 ? sc.blade : part_ == 1 ? sc.hub : sc.casing;
        const auto& lines = part_ == 0 ? sc.blade_lines : part_ == 1 ? sc.hub_lines : sc.casing_lines;
        const auto& edges = part_ == 0 ? sc.blade_edges : part_ == 1 ? sc.hub_edges : sc.casing_edges;
        const int copies = part_ == 0 ? std::max(1, sc.copies) : 1;
        std::vector<gp_Trsf> turn;
        for (int k = 0; k < copies; ++k) turn.push_back(BladeGen::Design::theta_rotation(k * sc.pitch));
        auto pnt = [](const BladeGen::Math::Vec3& p, const gp_Trsf& t) { return gp_Pnt(p.x(), p.y(), p.z()).Transformed(t); };
        auto polylines = [&](const std::vector<std::vector<BladeGen::Math::Vec3>>& ls, const Handle(Graphic3d_AspectLine3d)& aspect) {
            int nv = 0;
            for (const auto& l : ls) nv += static_cast<int>(l.size());
            if (nv == 0) return;
            Handle(Graphic3d_ArrayOfPolylines) arr = new Graphic3d_ArrayOfPolylines(nv * copies, static_cast<int>(ls.size()) * copies);
            for (const auto& t : turn)
                for (const auto& l : ls) {
                    arr->AddBound(static_cast<int>(l.size()));
                    for (const auto& p : l) arr->AddVertex(pnt(p, t));
                }
            Handle(Graphic3d_Group) g = prs->NewGroup();
            g->SetGroupPrimitivesAspect(aspect);
            g->AddPrimitiveArray(arr);
        };
        if (mode == 0) {
            polylines(lines, wire_);
            polylines(edges, wire_);
            return;
        }
        const int nv = static_cast<int>(m.p.size()), nt = static_cast<int>(m.tri.size());
        Handle(Graphic3d_ArrayOfTriangles) tri = new Graphic3d_ArrayOfTriangles(nv * copies, 3 * nt * copies, Graphic3d_ArrayFlags_VertexNormal);
        for (std::size_t k = 0; k < turn.size(); ++k) {
            const auto& t = turn[k];
            for (std::size_t v = 0; v < m.p.size(); ++v)
                tri->AddVertex(pnt(m.p[v], t), gp_Dir(m.n[v].x(), m.n[v].y(), m.n[v].z()).Transformed(t));
            const int base = static_cast<int>(k) * nv + 1;  // 1-based
            for (const auto& [a, b, c] : m.tri) tri->AddEdges(base + a, base + b, base + c);
        }
        Handle(Graphic3d_Group) g = prs->NewGroup();
        g->SetGroupPrimitivesAspect(fill_);
        g->AddPrimitiveArray(tri);
        if (edges_) polylines(edges, edge_);
    }

    void ComputeSelection(const Handle(SelectMgr_Selection)&, const Standard_Integer) override {}

private:
    std::shared_ptr<const BladeGen::Design::SceneMesh> scene_;
    std::size_t part_;
    Quantity_Color color_;
    Graphic3d_NameOfMaterial material_;
    bool edges_ = false;
    Handle(Graphic3d_AspectFillArea3d) fill_;
    Handle(Graphic3d_AspectLine3d) wire_, edge_;
};

/// One OpenGL driver (display connection, shared GL context) and, on Windows, one window class for
/// every view of the process.  Created on first use and kept to the end of the process: a window
/// class per view failed once views had been closed ("Unable to register window class").
struct SharedGl
{
    Handle(Aspect_DisplayConnection) disp;
    Handle(OpenGl_GraphicDriver) driver;
#if defined(_WIN32)
    Handle(WNT_WClass) cls;
#endif
    // the owner of the driver's shared GL context: a view on a hidden window that lives as long
    // as the process (the first view's context would die with that view's window)
    Handle(V3d_Viewer) viewer;
    Handle(V3d_View) view;
    Handle(Aspect_Window) native;
    Handle(Aspect_NeutralWindow) win;
    QString error;
};

/// A hidden native window (only its OpenGL context is used; the views render offscreen).
Handle(Aspect_Window) hiddenWindow(const SharedGl& gl)
{
#if defined(_WIN32)
    return new WNT_Window("BladeGen OCCT", gl.cls, WS_POPUP, 0, 0, 16, 16, Quantity_NOC_BLACK);
#elif defined(__APPLE__)
    (void)gl;
    return new Cocoa_Window("BladeGen OCCT", 0, 0, 16, 16);
#else
    return new Xw_Window(gl.disp, "BladeGen OCCT", 0, 0, 16, 16);
#endif
}

Handle(Aspect_NeutralWindow) neutralOn(const Handle(Aspect_Window)& native)
{
    Handle(Aspect_NeutralWindow) w = new Aspect_NeutralWindow();
    w->SetNativeHandles(native->NativeHandle(), 0, native->NativeFBConfig());
    w->SetSize(64, 64);
    w->SetVirtual(Standard_True);
    return w;
}

SharedGl& sharedGl()
{
    static SharedGl* g = [] {
        auto* s = new SharedGl;  // never deleted: no teardown order against Qt or the GL driver at exit
        if (qEnvironmentVariableIsSet("BLADEGEN_NO_OPENGL")) {  // broken drivers; tests of the fallback
            s->error = "switched off by the environment variable BLADEGEN_NO_OPENGL";
            return s;
        }
        try {
            s->disp = new Aspect_DisplayConnection();
            s->driver = new OpenGl_GraphicDriver(s->disp);
#if defined(_WIN32)
            s->cls = new WNT_WClass("BladeGenOcctView", reinterpret_cast<Standard_Address>(DefWindowProcW), CS_OWNDC);
#endif
            s->viewer = new V3d_Viewer(s->driver);
            s->view = s->viewer->CreateView();
            s->native = hiddenWindow(*s);
            s->win = neutralOn(s->native);
            s->view->SetWindow(s->win);
            // the views render offscreen into frame buffers: OpenGL 2.1 with FBOs.  A software
            // driver (Windows "GDI Generic" 1.1, e.g. a remote desktop or a CI machine) makes OCCT
            // fall back to the fixed pipeline, and the dumps come out blank: say so instead.
            const Handle(OpenGl_Context)& c = s->driver->GetSharedContext();
            if (c.IsNull()) s->error = "no OpenGL context";
            else if (!c->IsGlGreaterEqual(2, 1) || c->arbFBO == nullptr)
                s->error = QString("the OpenGL driver offers version %1.%2%3; the view needs OpenGL 2.1 with frame buffers "
                                   "(update the graphics driver)")
                               .arg(c->VersionMajor())
                               .arg(c->VersionMinor())
                               .arg(c->arbFBO == nullptr ? " without frame buffers" : "");
        } catch (const Standard_Failure& e) {
            s->error = QString::fromUtf8(e.GetMessageString());
        } catch (const std::exception& e) {
            s->error = QString::fromUtf8(e.what());
        }
        return s;
    }();
    return *g;
}

const Quantity_Color& partColor(std::size_t k)
{
    static const Quantity_Color colors[] = {Quantity_Color(0.74, 0.20, 0.17, Quantity_TOC_sRGB), Quantity_Color(0.36, 0.60, 0.40, Quantity_TOC_sRGB),
                                            Quantity_Color(0.62, 0.64, 0.70, Quantity_TOC_sRGB)};
    return colors[k];
}

}  // namespace

struct OcctView::Impl : public AIS_ViewController
{
    OcctView* w = nullptr;
    Handle(V3d_Viewer) viewer;
    Handle(V3d_View) view;
    Handle(AIS_InteractiveContext) ctx;
    Handle(Aspect_Window) native;  // hidden platform window: owns the OpenGL context
    Handle(Aspect_NeutralWindow) win;
    std::array<Handle(AIS_InteractiveObject), 3> parts;
    std::array<Look, 3> looks{Look::Solid, Look::Solid, Look::Shaded};
    std::shared_ptr<const BladeGen::Design::Scene> scene;
    std::shared_ptr<const BladeGen::Design::SceneMesh> mesh;
    QString reason;
    QImage frame;
    bool dirty = true;
    double theta_blade = 0.0;

    bool started = false;

    ~Impl() override
    {
        // the view's GL resources first, while its hidden window still exists; then the window
        try {
            if (!ctx.IsNull())
                for (auto& p : parts)
                    if (!p.IsNull()) ctx->Remove(p, Standard_False);
            parts = {};
            if (!view.IsNull()) view->Remove();
        } catch (const Standard_Failure&) {
        }
        ctx.Nullify();
        view.Nullify();
        viewer.Nullify();
        win.Nullify();
        native.Nullify();
    }

    bool ok() const { return reason.isEmpty() && !view.IsNull(); }

    /// Starts the viewer on first use (only views that show the 3D Solid get an OpenGL context).
    bool ensure()
    {
        if (!started) {
            started = true;
            init();
        }
        return ok();
    }

    void init()
    {
        auto& gl = sharedGl();
        if (!gl.error.isEmpty()) {
            reason = "OpenGL viewer could not start: " + gl.error;
            return;
        }
        try {
            viewer = new V3d_Viewer(gl.driver);
            viewer->SetDefaultLights();
            viewer->SetLightOn();
            ctx = new AIS_InteractiveContext(viewer);
            view = viewer->CreateView();
            view->SetImmediateUpdate(Standard_False);
            native = hiddenWindow(gl);
            win = neutralOn(native);
            view->SetWindow(win);
            view->SetBgGradientColors(Quantity_Color(0.98, 0.98, 0.99, Quantity_TOC_sRGB), Quantity_Color(0.78, 0.81, 0.87, Quantity_TOC_sRGB),
                                      Aspect_GradientFillMethod_Vertical);
            // machine convention, as the 3D View: internal x (axial) is Z, y (r sin theta) is X, z (r cos theta) is Y
            view->Trihedron()->SetLabels("Z", "X", "Y");
            view->Trihedron()->SetArrowsColor(Quantity_Color(0.0, 0.0, 0.86, Quantity_TOC_sRGB), Quantity_Color(0.78, 0.0, 0.0, Quantity_TOC_sRGB),
                                              Quantity_Color(0.0, 0.59, 0.0, Quantity_TOC_sRGB));
            view->TriedronDisplay(Aspect_TOTP_LEFT_LOWER, Quantity_NOC_GRAY20, 0.08, V3d_ZBUFFER);
            view->ChangeRenderingParams().IsAntialiasingEnabled = Standard_True;
            view->ChangeRenderingParams().NbMsaaSamples = 4;
            SetAllowHighlight(false);
            SetAllowDragging(false);
            ChangeMouseGestureMap().UnBind(Aspect_VKeyMouse_RightButton);  // right click: the view menu
            ChangeMouseGestureMap().UnBind(static_cast<unsigned int>(Aspect_VKeyMouse_LeftButton) | static_cast<unsigned int>(Aspect_VKeyFlags_ALT));
            Image_PixMap probe;  // offscreen rendering works (frame buffers, OpenGL >= 2.1)?
            if (!view->ToPixMap(probe, 32, 32, Graphic3d_BT_RGBA)) reason = "offscreen OpenGL rendering failed (no frame buffer support)";
        } catch (const Standard_Failure& e) {
            reason = QString("OpenGL viewer could not start: %1").arg(e.GetMessageString());
        } catch (const std::exception& e) {
            reason = QString("OpenGL viewer could not start: %1").arg(e.what());
        }
        if (!reason.isEmpty()) {
            view.Nullify();
            ctx.Nullify();
            viewer.Nullify();
        }
    }

    void handleViewRedraw(const Handle(AIS_InteractiveContext)&, const Handle(V3d_View)&) override
    {
        dirty = true;
        w->update();
    }

    void resize(int W, int H)
    {
        if (!ok() || W < 2 || H < 2) return;
        if (win->SetSize(W, H)) {
            view->MustBeResized();
            dirty = true;
        }
    }

    bool render(int W, int H)
    {
        if (!ok() || W < 2 || H < 2) return false;
        resize(W, H);
        try {
            Image_PixMap pix;
            pix.InitZero(Image_Format_RGBA, static_cast<Standard_Size>(W), static_cast<Standard_Size>(H));
            V3d_ImageDumpOptions o;
            o.Width = W;
            o.Height = H;
            o.BufferType = Graphic3d_BT_RGBA;
            o.ToAdjustAspect = Standard_True;
            if (!view->ToPixMap(pix, o)) return false;
            QImage img(W, H, QImage::Format_RGBX8888);
            const auto rowBytes = static_cast<std::size_t>(W) * 4;
            for (int y = 0; y < H; ++y)  // Row() counts from the top whatever the storage order
                std::memcpy(img.scanLine(y), pix.Row(static_cast<Standard_Size>(y)), rowBytes);
            frame = img;
            dirty = false;
            return true;
        } catch (const Standard_Failure&) {
            return false;
        }
    }

    void applyLook(std::size_t k)
    {
        auto& s = parts[k];
        if (s.IsNull()) return;
        const Look l = looks[k];
        if (l == Look::Hidden) {
            if (ctx->IsDisplayed(s)) ctx->Erase(s, Standard_False);
            return;
        }
        const int mode = l == Look::Wireframe ? AIS_WireFrame : AIS_Shaded;
        const double transparency = l == Look::Shaded ? (k == 0 ? 0.35 : 0.6) : 0.0;
        const bool boundary = l == Look::Solid;
        if (Handle(MeshObject) mo = Handle(MeshObject)::DownCast(s); !mo.IsNull()) {
            mo->setLook(boundary, transparency);  // the presentation is cheap: recompute it
            if (ctx->DisplayStatus(s) == PrsMgr_DisplayStatus_None) {
                s->SetDisplayMode(mode);
                ctx->Display(s, mode, -1, Standard_False);
                return;
            }
            if (s->DisplayMode() != mode) ctx->SetDisplayMode(s, mode, Standard_False);
            ctx->Redisplay(s, Standard_False);
            if (!ctx->IsDisplayed(s)) ctx->Display(s, Standard_False);
            return;
        }
        if (ctx->DisplayStatus(s) == PrsMgr_DisplayStatus_None) {
            // first display: set everything, compute the presentation once
            s->Attributes()->SetFaceBoundaryDraw(boundary);
            s->SetTransparency(transparency);
            s->SetDisplayMode(mode);
            ctx->Display(s, mode, -1, Standard_False);
            return;
        }
        // a look change: only what differs is recomputed
        const bool redo = s->Attributes()->FaceBoundaryDraw() != boundary;
        s->Attributes()->SetFaceBoundaryDraw(boundary);
        if (s->DisplayMode() != mode) ctx->SetDisplayMode(s, mode, Standard_False);
        ctx->SetTransparency(s, transparency, Standard_False);
        if (redo) ctx->Redisplay(s, Standard_False);
        if (!ctx->IsDisplayed(s)) ctx->Display(s, Standard_False);
    }

    Handle(AIS_Shape) makePart(const TopoDS_Shape& shape, std::size_t k)
    {
        Handle(AIS_Shape) s = new AIS_Shape(shape);
        s->SetColor(partColor(k));
        s->SetMaterial(Graphic3d_MaterialAspect(k == 0 ? Graphic3d_NameOfMaterial_Plastified : Graphic3d_NameOfMaterial_Satin));
        const Handle(Prs3d_Drawer)& dr = s->Attributes();
        dr->SetTypeOfDeflection(Aspect_TOD_ABSOLUTE);
        dr->SetMaximalChordialDeviation(k == 0 ? kBladeDeflection : kWallDeflection);
        dr->SetDeviationAngle(k == 0 ? kBladeAngle : kWallAngle);
        dr->SetFaceBoundaryAspect(new Prs3d_LineAspect(Quantity_NOC_GRAY15, Aspect_TOL_SOLID, 1.0));
        const bool ring = scene && scene->copies > 1;
        const Quantity_Color iso = k == 0 ? Quantity_Color(0.45, 0.08, 0.06, Quantity_TOC_sRGB) : Quantity_Color(0.25, 0.27, 0.30, Quantity_TOC_sRGB);
        // u around the profile / around the axis, v along the span / along the wall
        dr->SetUIsoAspect(new Prs3d_IsoAspect(iso, Aspect_TOL_SOLID, 1.0, k == 0 ? 12 : (ring ? 48 : 6)));
        dr->SetVIsoAspect(new Prs3d_IsoAspect(iso, Aspect_TOL_SOLID, 1.0, k == 0 ? 6 : 10));
        return s;
    }

    void setScene(std::shared_ptr<const BladeGen::Design::Scene> s, bool refit)
    {
        if (!ensure()) return;
        for (auto& p : parts)
            if (!p.IsNull()) ctx->Remove(p, Standard_False);
        scene = std::move(s);
        mesh.reset();
        if (!scene) {
            parts = {};
            dirty = true;
            return;
        }
        theta_blade = scene->theta_blade;
        parts[0] = makePart(scene->blades, 0);
        parts[1] = makePart(scene->hub, 1);
        parts[2] = makePart(scene->casing, 2);
        for (std::size_t k = 0; k < 3; ++k) applyLook(k);
        if (refit) direction(Direction::Iso);
        dirty = true;
    }

    void setMesh(std::shared_ptr<const BladeGen::Design::SceneMesh> m, bool refit)
    {
        if (!ensure()) return;
        for (auto& p : parts)
            if (!p.IsNull()) ctx->Remove(p, Standard_False);
        mesh = std::move(m);
        scene.reset();
        parts = {};
        dirty = true;
        if (!mesh) return;
        theta_blade = mesh->theta_blade;
        for (std::size_t k = 0; k < 3; ++k) {
            parts[k] = new MeshObject(mesh, k, partColor(k), k == 0 ? Graphic3d_NameOfMaterial_Plastified : Graphic3d_NameOfMaterial_Satin);
            applyLook(k);
        }
        if (refit) direction(Direction::Iso);
    }

    void direction(Direction d)
    {
        if (!ok()) return;
        const double c = std::cos(theta_blade), sn = std::sin(theta_blade);
        // radial, tangential (+theta) and axial unit vectors at the blade
        const gp_Vec er(0, sn, c), et(0, c, -sn), ex(1, 0, 0);
        gp_Vec eye, up = er;
        switch (d) {
        case Direction::Iso: eye = 0.55 * er - 0.70 * et - 0.45 * ex; break;
        case Direction::Axial: eye = -ex; break;
        case Direction::Meridional: eye = -et; break;
        case Direction::Radial:
            eye = er;
            up = -et;
            break;
        }
        view->SetProj(eye.X(), eye.Y(), eye.Z());
        view->SetUp(up.X(), up.Y(), up.Z());
        view->FitAll(0.04, Standard_False);
        dirty = true;
    }
};

#else  // no OCCT visualization toolkits

struct OcctView::Impl
{
    OcctView* w = nullptr;
    std::array<Look, 3> looks{Look::Solid, Look::Solid, Look::Shaded};
    QString reason = "this build has no OpenCASCADE visualization toolkits (TKV3d, TKOpenGl, TKService)";
    bool ok() const { return false; }
    bool ensure() { return false; }
};

#endif

OcctView::OcctView(QWidget* parent) : QWidget(parent), d_(std::make_unique<Impl>())
{
    d_->w = this;
    setMouseTracking(false);
    setMinimumSize(80, 60);
    setFocusPolicy(Qt::StrongFocus);
}

OcctView::~OcctView() = default;

bool OcctView::available() const { return d_->ensure(); }

QString OcctView::unavailableReason() const { return d_->reason; }

void OcctView::prepare(const BladeGen::Design::Scene& s)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    // the copies share the blade's faces: one triangulation for all of them
    BRepMesh_IncrementalMesh(s.shapes.blade, kBladeDeflection, Standard_False, kBladeAngle, Standard_True);
    BRepMesh_IncrementalMesh(s.hub, kWallDeflection, Standard_False, kWallAngle, Standard_True);
    BRepMesh_IncrementalMesh(s.casing, kWallDeflection, Standard_False, kWallAngle, Standard_True);
#else
    (void)s;
#endif
}

void OcctView::setScene(std::shared_ptr<const BladeGen::Design::Scene> s, bool refit)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    d_->setScene(std::move(s), refit);
#else
    (void)s;
    (void)refit;
#endif
    update();
}

void OcctView::setMeshScene(std::shared_ptr<const BladeGen::Design::SceneMesh> m, bool refit)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    d_->setMesh(std::move(m), refit);
#else
    (void)m;
    (void)refit;
#endif
    update();
}

void OcctView::setMessage(const QString& m)
{
    message_ = m;
    update();
}

void OcctView::setLook(Part p, Look l)
{
    const auto k = static_cast<std::size_t>(p);
    if (p == Part::Blade && l == Look::Hidden) l = Look::Solid;  // the blade is always shown
    d_->looks[k] = l;
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok()) {
        d_->applyLook(k);
        d_->dirty = true;
    }
#endif
    update();
}

OcctView::Look OcctView::look(Part p) const { return d_->looks[static_cast<std::size_t>(p)]; }

void OcctView::setDirection(Direction dir)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    d_->direction(dir);
#else
    (void)dir;
#endif
    update();
}

void OcctView::fitAll()
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok()) {
        d_->view->FitAll(0.04, Standard_False);
        d_->dirty = true;
    }
#endif
    update();
}

QImage OcctView::image()
{
    const double r = devicePixelRatioF();
    const int W = std::max(2, static_cast<int>(std::lround(width() * r))), H = std::max(2, static_cast<int>(std::lround(height() * r)));
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ensure() && d_->render(W, H)) {
        QImage img = d_->frame.convertToFormat(QImage::Format_RGB32);
        if (!message_.isEmpty()) {
            QPainter p(&img);
            p.setPen(QColor(40, 40, 40));
            p.drawText(img.rect().adjusted(8, 6, -8, -6), Qt::AlignTop | Qt::AlignLeft | Qt::TextWordWrap, message_);
        }
        return img;
    }
#else
    (void)W;
    (void)H;
#endif
    return grab().toImage();
}

void OcctView::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    const double r = devicePixelRatioF();
    const int W = static_cast<int>(std::lround(width() * r)), H = static_cast<int>(std::lround(height() * r));
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ensure()) {
        if (d_->dirty || d_->frame.width() != W || d_->frame.height() != H) d_->render(W, H);
        if (!d_->frame.isNull()) {
            QImage f = d_->frame;
            f.setDevicePixelRatio(r);
            p.drawImage(QPointF(0, 0), f);
        }
    }
#else
    (void)W;
    (void)H;
#endif
    if (!d_->ok()) {
        p.fillRect(rect(), QColor(236, 238, 242));
        p.setPen(QColor(120, 30, 30));
        p.drawText(rect().adjusted(12, 12, -12, -12), Qt::AlignCenter | Qt::TextWordWrap,
                   "3-D solid view unavailable:\n" + d_->reason + "\n\nThe 3D View (mesh) and the STEP export still work.");
    }
    if (!message_.isEmpty()) {
        p.setPen(QColor(40, 40, 40));
        p.drawText(rect().adjusted(8, 6, -8, -6), Qt::AlignTop | Qt::AlignLeft | Qt::TextWordWrap, message_);
    }
}

void OcctView::resizeEvent(QResizeEvent*)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    const double r = devicePixelRatioF();
    d_->resize(static_cast<int>(std::lround(width() * r)), static_cast<int>(std::lround(height() * r)));
#endif
}

#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
namespace {
Graphic3d_Vec2i pixel(const QWidget* w, const QPointF& p)
{
    const double r = w->devicePixelRatioF();
    return Graphic3d_Vec2i(static_cast<int>(std::lround(p.x() * r)), static_cast<int>(std::lround(p.y() * r)));
}
Aspect_VKeyMouse buttons(Qt::MouseButtons b)
{
    Aspect_VKeyMouse m = Aspect_VKeyMouse_NONE;
    if (b & Qt::LeftButton) m |= Aspect_VKeyMouse_LeftButton;
    if (b & Qt::MiddleButton) m |= Aspect_VKeyMouse_MiddleButton;
    return m;
}
Aspect_VKeyFlags flags(Qt::KeyboardModifiers k)
{
    Aspect_VKeyFlags f = Aspect_VKeyFlags_NONE;
    if (k & Qt::ShiftModifier) f |= Aspect_VKeyFlags_SHIFT;
    if (k & Qt::ControlModifier) f |= Aspect_VKeyFlags_CTRL;
    if (k & Qt::AltModifier) f |= Aspect_VKeyFlags_ALT;
    return f;
}
}  // namespace
#endif

void OcctView::mousePressEvent(QMouseEvent* e)
{
    press_ = e->pos();
    moved_ = false;
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok() && e->button() != Qt::RightButton) {
        d_->UpdateMouseButtons(pixel(this, e->position()), buttons(e->buttons()), flags(e->modifiers()), false);
        d_->FlushViewEvents(d_->ctx, d_->view, Standard_True);
    }
#endif
}

void OcctView::mouseMoveEvent(QMouseEvent* e)
{
    if ((e->pos() - press_).manhattanLength() > 3) moved_ = true;
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok() && (e->buttons() & (Qt::LeftButton | Qt::MiddleButton))) {
        d_->UpdateMousePosition(pixel(this, e->position()), buttons(e->buttons()), flags(e->modifiers()), false);
        d_->FlushViewEvents(d_->ctx, d_->view, Standard_True);
    }
#endif
}

void OcctView::mouseReleaseEvent(QMouseEvent* e)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok() && e->button() != Qt::RightButton) {
        d_->UpdateMouseButtons(pixel(this, e->position()), buttons(e->buttons()), flags(e->modifiers()), false);
        d_->FlushViewEvents(d_->ctx, d_->view, Standard_True);
    }
#else
    (void)e;
#endif
}

void OcctView::mouseDoubleClickEvent(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton) fitAll();
}

void OcctView::wheelEvent(QWheelEvent* e)
{
#if defined(BLADEGEN_HAVE_OCCT_VIEWER)
    if (d_->ok()) {
        d_->UpdateMouseScroll(Aspect_ScrollDelta(pixel(this, e->position()), e->angleDelta().y() / 8.0, flags(e->modifiers())));
        d_->FlushViewEvents(d_->ctx, d_->view, Standard_True);
    }
#endif
    e->accept();
}

void OcctView::contextMenuEvent(QContextMenuEvent* e) { emit contextMenuWanted(e->globalPos()); }

}  // namespace Designer
