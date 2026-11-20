#include "occt/Fillet.hpp"

#include "blade/Blade.hpp"
#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"
#include "occt/OcctBlade.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_Status.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Solid.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <sstream>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace BladeGen::Occt {

namespace {

using Math::Vec2;
using Math::Vec3;

gp_Pnt pnt(const Vec3& p) { return {p.x(), p.y(), p.z()}; }
Vec3 vec(const gp_Pnt& p) { return {p.X(), p.Y(), p.Z()}; }
Vec3 vec(const gp_Vec& p) { return {p.X(), p.Y(), p.Z()}; }

std::vector<TopoDS_Face> faces_of(const TopoDS_Shape& s)
{
    std::vector<TopoDS_Face> f;
    for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next()) f.push_back(TopoDS::Face(ex.Current()));
    return f;
}

double shape_distance(const TopoDS_Shape& a, const TopoDS_Shape& b)
{
    BRepExtrema_DistShapeShape d(a, b);
    return d.IsDone() ? d.Value() : std::numeric_limits<double>::infinity();
}

Handle(Geom2d_Line) line2d(double u0, double v0, double du, double dv)
{
    return new Geom2d_Line(gp_Pnt2d(u0, v0), gp_Dir2d(du, dv));
}

/// Closed 2-D curve through uv at parameters t (t.size() == uv.size() + 1); the end
/// point is uv[0] shifted by (period_u, 0); periodic when period_u == 0, else C1
/// across the closure (Python closed_pcurve).
Handle(Geom2d_BSplineCurve) closed_pcurve(const std::vector<Vec2>& uv, const std::vector<double>& t, double period_u)
{
    const int n = static_cast<int>(uv.size());
    if (period_u == 0.0) {
        Handle(TColgp_HArray1OfPnt2d) hp = new TColgp_HArray1OfPnt2d(1, n);
        for (int i = 0; i < n; ++i) hp->SetValue(i + 1, gp_Pnt2d(uv[static_cast<std::size_t>(i)].x(), uv[static_cast<std::size_t>(i)].y()));
        Handle(TColStd_HArray1OfReal) hr = new TColStd_HArray1OfReal(1, n + 1);
        for (int i = 0; i <= n; ++i) hr->SetValue(i + 1, t[static_cast<std::size_t>(i)]);
        Geom2dAPI_Interpolate api(hp, hr, Standard_True, 1e-12);
        api.Perform();
        if (!api.IsDone()) throw std::runtime_error("periodic pcurve interpolation failed");
        return api.Curve();
    }
    Handle(TColgp_HArray1OfPnt2d) hp = new TColgp_HArray1OfPnt2d(1, n + 1);
    for (int i = 0; i < n; ++i) hp->SetValue(i + 1, gp_Pnt2d(uv[static_cast<std::size_t>(i)].x(), uv[static_cast<std::size_t>(i)].y()));
    hp->SetValue(n + 1, gp_Pnt2d(uv[0].x() + period_u, uv[0].y()));
    Handle(TColStd_HArray1OfReal) hr = new TColStd_HArray1OfReal(1, n + 1);
    for (int i = 0; i <= n; ++i) hr->SetValue(i + 1, t[static_cast<std::size_t>(i)]);
    const Vec2 shift(period_u, 0.0);
    const Vec2 tan = (uv[1] - uv[static_cast<std::size_t>(n - 1)] + shift) /
                     (t[1] - t[0] + t[static_cast<std::size_t>(n)] - t[static_cast<std::size_t>(n - 1)]);
    Geom2dAPI_Interpolate api(hp, hr, Standard_False, 1e-12);
    api.Load(gp_Vec2d(tan.x(), tan.y()), gp_Vec2d(tan.x(), tan.y()), Standard_False);
    api.Perform();
    if (!api.IsDone()) throw std::runtime_error("pcurve interpolation failed");
    return api.Curve();
}

/// BRepCheck statuses of a face, its wires and edges (for error messages).
std::string describe_invalid(const TopoDS_Shape& f)
{
    BRepCheck_Analyzer an(f);
    std::string out;
    auto add = [&](const TopoDS_Shape& sub, const char* what) {
        const Handle(BRepCheck_Result) r = an.Result(sub);
        if (r.IsNull()) return;
        for (BRepCheck_ListIteratorOfListOfStatus it(r->Status()); it.More(); it.Next())
            if (it.Value() != BRepCheck_NoError) out += std::string(" ") + what + ":" + std::to_string(static_cast<int>(it.Value()));
        r->InitContextIterator();
        for (; r->MoreShapeInContext(); r->NextShapeInContext())
            for (BRepCheck_ListIteratorOfListOfStatus it(r->StatusOnShape()); it.More(); it.Next())
                if (it.Value() != BRepCheck_NoError) out += std::string(" ") + what + "@ctx:" + std::to_string(static_cast<int>(it.Value()));
    };
    add(f, "face");
    for (TopExp_Explorer w(f, TopAbs_WIRE); w.More(); w.Next()) add(w.Current(), "wire");
    for (TopExp_Explorer e(f, TopAbs_EDGE); e.More(); e.Next()) add(e.Current(), "edge");
    for (TopExp_Explorer v(f, TopAbs_VERTEX); v.More(); v.Next()) add(v.Current(), "vertex");
    return out.empty() ? " (no status)" : out;
}

/// UV path of a wire on a face: orientation and pcurve end points of every edge.
std::string wire_uv(const TopoDS_Wire& w, const TopoDS_Face& f)
{
    std::ostringstream o;
    o.precision(5);
    for (TopoDS_Iterator it(w); it.More(); it.Next()) {
        const TopoDS_Edge e = TopoDS::Edge(it.Value());
        double a, b;
        const Handle(Geom2d_Curve) c = BRep_Tool::CurveOnSurface(e, f, a, b);
        if (c.IsNull()) { o << " [no pcurve]"; continue; }
        gp_Pnt2d p0 = c->Value(a), p1 = c->Value(b);
        if (e.Orientation() == TopAbs_REVERSED) std::swap(p0, p1);
        o << " [" << (e.Orientation() == TopAbs_FORWARD ? "F" : "R") << " (" << p0.X() << "," << p0.Y() << ")->(" << p1.X() << "," << p1.Y() << ")]";
    }
    return o.str();
}

/// ShapeFix_Face on a hand-built face that must keep all its (shared) edges.
TopoDS_Face fixed_face(const TopoDS_Face& f)
{
    ShapeFix_Face fix(f);
    fix.Perform();
    TopoDS_Face out = fix.Face();
    for (TopExp_Explorer e(f, TopAbs_EDGE); e.More(); e.Next()) {
        bool kept = false;
        for (TopExp_Explorer x(out, TopAbs_EDGE); x.More() && !kept; x.Next()) kept = e.Current().IsSame(x.Current());
        if (!kept) throw std::runtime_error("ShapeFix replaced a shared edge of a rebuilt face");
    }
    return BRepCheck_Analyzer(out).IsValid() ? out : TopoDS_Face();
}

TopoDS_Shape solid_from(const std::vector<TopoDS_Face>& faces)
{
    BRep_Builder b;
    TopoDS_Shell sh;
    b.MakeShell(sh);
    for (const auto& f : faces) b.Add(sh, f);
    sh.Closed(Standard_True);
    ShapeAnalysis_FreeBounds fb(sh);
    int free_edges = 0;
    for (TopExp_Explorer e(fb.GetClosedWires(), TopAbs_EDGE); e.More(); e.Next()) ++free_edges;
    for (TopExp_Explorer e(fb.GetOpenWires(), TopAbs_EDGE); e.More(); e.Next()) ++free_edges;
    if (free_edges) throw std::runtime_error("fillet shell has " + std::to_string(free_edges) + " free edges");
    TopoDS_Solid so;
    b.MakeSolid(so);
    b.Add(so, sh);
    ShapeFix_Solid fix(so);
    fix.Perform();
    return fix.Solid();
}

}  // namespace

FilletResult rolling_ball_fillet(const TopoDS_Shape& blade, const TopoDS_Face& blade_face, const TopoDS_Shape& hub_solid,
                                 const Blade::MeridionalCurve& hub, const FilletOptions& opt)
{
    FilletResult res;
    const double r = opt.radius;
    if (!(r > 0.0)) throw std::invalid_argument("rolling_ball_fillet: radius must be > 0");
    if (opt.stations < 16 || opt.arc_rows < 3) throw std::invalid_argument("rolling_ball_fillet: too few stations / rows");
    const double tol = opt.tolerance;
    BRep_Builder b;

    // ── identify the faces ───────────────────────────────────────────────────
    const Handle(Geom_Surface) bsurf = BRep_Tool::Surface(blade_face);
    TopoDS_Face L, R, T;
    std::vector<TopoDS_Face> caps;
    for (const auto& f : faces_of(blade)) {
        if (BRep_Tool::Surface(f) == bsurf && L.IsNull()) L = f;
        else caps.push_back(f);
    }
    if (L.IsNull() || caps.size() != 2) throw std::runtime_error("blade solid: expected one lateral face and two caps");
    if (shape_distance(caps[0], hub_solid) <= shape_distance(caps[1], hub_solid)) { R = caps[0]; T = caps[1]; }
    else { R = caps[1]; T = caps[0]; }
    TopoDS_Face H;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& f : faces_of(hub_solid)) {
        const double d = shape_distance(f, R);
        if (d < best) { best = d; H = f; }
    }
    if (best > 1e-5) throw std::runtime_error("the blade root does not lie on the hub");

    // ── fluid-side normal of the blade surface ───────────────────────────────
    double u1, u2, v1, v2;
    bsurf->Bounds(u1, u2, v1, v2);
    auto point_normal = [&](double u, double v) {
        gp_Pnt p;
        gp_Vec du, dv;
        bsurf->D1(u, v, p, du, dv);
        gp_Vec n = du.Crossed(dv);
        n.Normalize();
        return std::pair<Vec3, Vec3>{vec(p), vec(n)};
    };
    double side = 1.0;
    {
        // probe at the middle of the lateral face
        double fu1, fu2, fv1, fv2;
        BRepTools::UVBounds(L, fu1, fu2, fv1, fv2);
        const auto [p, n] = point_normal(0.5 * (fu1 + fu2), 0.5 * (fv1 + fv2));
        BRepClass3d_SolidClassifier cl(blade, pnt(p + 1e-3 * r * n), 1e-9);
        side = cl.State() == TopAbs_OUT ? 1.0 : -1.0;
    }

    // ── hub contour offset by r (meridional plane) ───────────────────────────
    std::vector<double> xo;
    std::vector<Math::Vec<1>> ro;
    const double Lh = hub.length();
    for (int i = 0; i <= 4000; ++i) {
        const double m = -0.5 * Lh + 2.0 * Lh * i / 4000.0;
        const Vec2 q = hub.at_m(m), t = hub.tangent_at_m(std::clamp(m, 0.0, Lh));
        const double x = q.x() - r * t.y();
        if (!xo.empty() && !(x > xo.back())) throw std::runtime_error("hub contour: fillet radius exceeds its curvature radius");
        xo.push_back(x);
        ro.push_back(Math::Vec<1>(q.y() + r * t.x()));
    }
    const auto roff = Math::interpolate_cubic<1>(xo, ro);
    auto gap = [&](const Vec3& c) { return std::hypot(c.y(), c.z()) - roff.eval(c.x())[0]; };
    auto hub_foot = [&](const Vec3& c) {
        auto f = [&](double m) {
            const Vec2 q = hub.at_m(m), t = hub.tangent_at_m(std::clamp(m, 0.0, Lh));
            return q.x() - r * t.y() - c.x();
        };
        const double m = Math::brent_root(f, -0.5 * Lh, 1.5 * Lh, 1e-15);
        const Vec2 q = hub.at_m(m);
        return Blade::cartesian(q.x(), q.y(), std::atan2(c.y(), c.z()));
    };

    // ── spine ─────────────────────────────────────────────────────────────────
    const int ns = opt.stations;
    std::vector<Vec3> spine, pb, ph;
    std::vector<Vec2> uvb;
    for (int i = 0; i < ns; ++i) {
        const double u = u1 + (u2 - u1) * i / ns;
        auto g = [&](double v) {
            const auto [p, n] = point_normal(u, v);
            return gap(p + r * side * n);
        };
        const double vs = Math::first_root(g, v1, v1 + 0.7 * (v2 - v1), 140, 1e-15);
        if (!std::isfinite(vs)) {
            res.message = "no ball position at blade station u = " + std::to_string(u);
            return res;
        }
        const auto [p, n] = point_normal(u, vs);
        const Vec3 c = p + r * side * n;
        spine.push_back(c);
        pb.push_back(p);
        ph.push_back(hub_foot(c));
        uvb.emplace_back(u, vs);
    }
    // ── arcs and the fillet surface ──────────────────────────────────────────
    const int M = opt.arc_rows;
    std::vector<std::vector<Vec3>> rows(static_cast<std::size_t>(M));
    res.min_arc_deg = 1e9;
    for (int i = 0; i < ns; ++i) {
        const auto si = static_cast<std::size_t>(i);
        const Vec3 e1 = (ph[si] - spine[si]) / r;
        const Vec3 w = pb[si] - spine[si];
        Vec3 e2 = w - w.dot(e1) * e1;
        e2.normalize();
        const double phi = std::atan2(w.dot(e2), w.dot(e1));
        res.min_arc_deg = std::min(res.min_arc_deg, Math::rad2deg(phi));
        res.max_arc_deg = std::max(res.max_arc_deg, Math::rad2deg(phi));
        res.max_ball_error = std::max({res.max_ball_error, std::abs((pb[si] - spine[si]).norm() - r),
                                       std::abs((ph[si] - spine[si]).norm() - r)});
        for (int k = 0; k < M; ++k) {
            const double f = phi * k / (M - 1);
            rows[static_cast<std::size_t>(k)].push_back(spine[si] + r * (std::cos(f) * e1 + std::sin(f) * e2));
        }
    }
    for (auto& row : rows) row.push_back(row.front());
    // parametrise every row by the blade station (uniform in blade u): the rows are
    // smooth functions of it, while their chord lengths differ strongly at the LE
    std::vector<double> ustat(static_cast<std::size_t>(ns) + 1);
    for (int i = 0; i <= ns; ++i) ustat[static_cast<std::size_t>(i)] = static_cast<double>(i) / ns;
    const auto fsk = skin_loops(rows, ustat);
    const Handle(Geom_BSplineSurface)& fsurf = fsk.surface;
    const std::vector<double>& t = fsk.u;  // fillet parameter of every station (t.back() = 1)

    // ── contact edges, shared by the three faces ─────────────────────────────
    TopoDS_Edge eh = BRepBuilderAPI_MakeEdge(fsurf->VIso(0.0)).Edge();
    TopoDS_Edge eb = BRepBuilderAPI_MakeEdge(fsurf->VIso(1.0)).Edge();
    const TopoDS_Vertex vh = TopExp::FirstVertex(eh), vb = TopExp::FirstVertex(eb);

    // ── hub face with the contact loop as a hole ─────────────────────────────
    const Handle(Geom_Surface) hsurf = BRep_Tool::Surface(H);
    double hu1, hu2, hv1, hv2;
    hsurf->Bounds(hu1, hu2, hv1, hv2);
    GeomAPI_ProjectPointOnSurf proj;
    proj.Init(hsurf, hu1, hu2, hv1, hv2);
    std::vector<Vec2> uvh;
    double dmax = 0.0;
    for (const auto& q : ph) {
        proj.Perform(pnt(q));
        if (!proj.IsDone() || proj.NbPoints() == 0) throw std::runtime_error("hub contact projection failed");
        double uu, vv;
        proj.LowerDistanceParameters(uu, vv);
        uvh.emplace_back(uu, vv);
        dmax = std::max(dmax, proj.LowerDistance());
    }
    if (dmax > 1e-5) throw std::runtime_error("hub contact curve is off the hub face");
    const auto pch = closed_pcurve(uvh, t, 0.0);
    double sarea = 0.0;
    for (std::size_t i = 0; i < uvh.size(); ++i) sarea += Math::cross2(uvh[i], uvh[(i + 1) % uvh.size()]);
    const double area_h = area(H);
    TopoDS_Face newH;
    TopAbs_Orientation hole_or = TopAbs_FORWARD;
    const TopAbs_Orientation rule = ((sarea > 0.0) == (H.Orientation() == TopAbs_FORWARD)) ? TopAbs_REVERSED : TopAbs_FORWARD;
    for (TopAbs_Orientation o : {rule, rule == TopAbs_REVERSED ? TopAbs_FORWARD : TopAbs_REVERSED}) {
        TopoDS_Face f = TopoDS::Face(H.EmptyCopied());
        for (TopExp_Explorer w(H, TopAbs_WIRE); w.More(); w.Next()) b.Add(f, w.Current());
        b.UpdateEdge(eh, pch, f, tol);
        TopoDS_Wire hw;
        b.MakeWire(hw);
        b.Add(hw, eh.Oriented(o));
        b.Add(f, hw);
        const TopoDS_Face c = fixed_face(f);
        if (!c.IsNull() && area(c) < area_h) { newH = c; hole_or = o; break; }
    }
    if (newH.IsNull()) throw std::runtime_error("could not rebuild the hub face with the contact loop as its hole");

    // ── blade face from the contact loop upward ──────────────────────────────
    TopTools_IndexedMapOfShape root_edges;
    TopExp::MapShapes(R, TopAbs_EDGE, root_edges);
    TopoDS_Edge seam_old;
    for (TopExp_Explorer e(L, TopAbs_EDGE); e.More(); e.Next()) {
        const TopoDS_Edge& ed = TopoDS::Edge(e.Current());
        if (BRep_Tool::IsClosed(ed, L) && !root_edges.Contains(ed)) { seam_old = ed; break; }
    }
    if (seam_old.IsNull()) throw std::runtime_error("blade face: seam not found");
    TopTools_IndexedMapOfShape root_vertices;
    TopExp::MapShapes(R, TopAbs_VERTEX, root_vertices);
    const TopoDS_Vertex s0 = TopExp::FirstVertex(seam_old), s1 = TopExp::LastVertex(seam_old);
    const bool seam_up = root_vertices.Contains(s0);  // seam runs root -> tip
    const TopoDS_Vertex tip_v = seam_up ? s1 : s0;
    GeomAPI_ProjectPointOnSurf pt(BRep_Tool::Pnt(tip_v), bsurf);
    double ut, vt;
    pt.LowerDistanceParameters(ut, vt);
    // stations start on the seam (u = u1); the contact vertex vb sits there
    TopoDS_Edge seam_b = BRepBuilderAPI_MakeEdge(bsurf->UIso(uvb[0].x()), vb, tip_v, uvb[0].y(), vt).Edge();
    auto seam_u = [&](TopAbs_Orientation o) {
        const TopoDS_Edge e = TopoDS::Edge(seam_old.Oriented(o));
        double f0, f1;
        const Handle(Geom2d_Curve) c2 = BRep_Tool::CurveOnSurface(e, L, f0, f1);
        return c2->Value(0.5 * (f0 + f1)).X();
    };
    TopoDS_Face nb0 = TopoDS::Face(L.EmptyCopied());
    // every occurrence of the new seam must lie on the side (u) of the old occurrence it
    // replaces; when the old seam ran tip -> root its orientations are swapped, so are the pcurves
    const double u_fwd = seam_u(seam_up ? TopAbs_FORWARD : TopAbs_REVERSED);
    const double u_rev = seam_u(seam_up ? TopAbs_REVERSED : TopAbs_FORWARD);
    b.UpdateEdge(seam_b, line2d(u_fwd, 0, 0, 1), line2d(u_rev, 0, 0, 1), nb0, tol);
    b.UpdateEdge(eb, closed_pcurve(uvb, t, u2 - u1), nb0, tol);
    TopoDS_Wire wb;
    b.MakeWire(wb);
    TopoDS_Wire wl;
    for (TopExp_Explorer w(L, TopAbs_WIRE); w.More(); w.Next()) wl = TopoDS::Wire(w.Current());
    bool eb_added = false;
    // the stored order and orientations of the old wire (BRepTools_WireExplorer may
    // report a reoriented traversal, which breaks the 2-D path of the seam)
    for (TopoDS_Iterator it(wl); it.More(); it.Next()) {
        const TopoDS_Edge e = TopoDS::Edge(it.Value());
        const TopAbs_Orientation o = e.Orientation();
        if (root_edges.Contains(e)) {
            if (eb_added) continue;
            double f0, f1;
            const Handle(Geom2d_Curve) c2 = BRep_Tool::CurveOnSurface(e, L, f0, f1);
            double du = c2->Value(f0 + 0.01 * (f1 - f0)).X() - c2->Value(f0).X();
            if (o == TopAbs_REVERSED) du = -du;
            b.Add(wb, eb.Oriented(du > 0 ? TopAbs_FORWARD : TopAbs_REVERSED));
            eb_added = true;
        } else if (e.IsSame(seam_old)) {
            b.Add(wb, seam_b.Oriented(seam_up ? o : TopAbs::Reverse(o)));
        } else {
            b.Add(wb, e);
        }
    }
    b.Add(nb0, wb);
    const TopoDS_Face newB = fixed_face(nb0);
    if (newB.IsNull())
        throw std::runtime_error("could not rebuild the blade face from the contact loop:" + describe_invalid(nb0) +
                                 (seam_up ? " (seam root->tip)" : " (seam tip->root)") +
                                 " root edges " + std::to_string(root_edges.Extent()) + wire_uv(wb, nb0) + " | old:" + wire_uv(wl, L));

    // ── fillet face ───────────────────────────────────────────────────────────
    TopoDS_Face ff;
    b.MakeFace(ff, fsurf, tol);
    TopoDS_Edge seam_f = BRepBuilderAPI_MakeEdge(fsurf->UIso(0.0), vh, vb, 0.0, 1.0).Edge();
    b.UpdateEdge(eh, line2d(0, 0, 1, 0), ff, tol);
    b.UpdateEdge(eb, line2d(0, 1, 1, 0), ff, tol);
    b.UpdateEdge(seam_f, line2d(1, 0, 0, 1), line2d(0, 0, 0, 1), ff, tol);
    TopoDS_Wire wf;
    b.MakeWire(wf);
    b.Add(wf, eh.Oriented(TopAbs_FORWARD));
    b.Add(wf, seam_f.Oriented(TopAbs_FORWARD));
    b.Add(wf, eb.Oriented(TopAbs_REVERSED));
    b.Add(wf, seam_f.Oriented(TopAbs_REVERSED));
    b.Add(ff, wf);
    {
        gp_Pnt p;
        gp_Vec a1, a2;
        fsurf->D1(0.37, 0.5, p, a1, a2);
        const Vec3 nrm = vec(a1.Crossed(a2));
        // outward from the material = toward the ball centre
        const auto it = std::min_element(spine.begin(), spine.end(), [&](const Vec3& a, const Vec3& c) {
            return (a - vec(p)).squaredNorm() < (c - vec(p)).squaredNorm();
        });
        if (nrm.dot(*it - vec(p)) < 0.0) ff = TopoDS::Face(ff.Reversed());
    }
    for (const auto* e : {&eh, &eb, &seam_b, &seam_f}) BRepLib::SameParameter(*e, 1e-5);
    if (!BRepCheck_Analyzer(ff).IsValid()) throw std::runtime_error("fillet face is invalid");
    res.fillet_face = ff;

    // ── G1 check at both contact loops ───────────────────────────────────────
    for (int i = 0; i < ns; i += 8) {
        const auto si = static_cast<std::size_t>(i);
        for (const auto& [v, wall] : {std::pair{0.0, ph[si]}, std::pair{1.0, pb[si]}}) {
            gp_Pnt q;
            gp_Vec a1, a2;
            fsurf->D1(t[si], v, q, a1, a2);
            const Vec3 n = vec(a1.Crossed(a2)).normalized();
            const Vec3 radial = (spine[si] - wall) / r;
            res.max_g1_deg = std::max(res.max_g1_deg, Math::rad2deg(std::acos(std::min(1.0, std::abs(n.dot(radial))))));
        }
    }

    // ── solids ────────────────────────────────────────────────────────────────
    std::vector<TopoDS_Face> hb;
    for (const auto& f : faces_of(hub_solid))
        if (!f.IsSame(H)) hb.push_back(f);
    hb.push_back(newH);
    hb.push_back(ff);
    hb.push_back(newB);
    hb.push_back(T);
    res.hub_blade = solid_from(hb);

    TopoDS_Face foot = TopoDS::Face(H.EmptyCopied());
    TopoDS_Wire wfoot;
    b.MakeWire(wfoot);
    b.Add(wfoot, eh.Oriented(hole_or == TopAbs_REVERSED ? TopAbs_FORWARD : TopAbs_REVERSED));
    b.Add(foot, wfoot);
    foot = fixed_face(foot);
    if (foot.IsNull()) throw std::runtime_error("could not build the blade footprint on the hub");
    res.blade_fillet = solid_from({TopoDS::Face(foot.Reversed()), ff, newB, T});

    res.ok = is_valid(res.hub_blade) && is_valid(res.blade_fillet) && volume(res.hub_blade) > 0.0;
    res.message = res.ok ? "rolling-ball fillet r = " + std::to_string(r) + " m, " + std::to_string(ns) + " stations"
                         : "fillet faces built but the solids are not valid";
    return res;
}

}  // namespace BladeGen::Occt
