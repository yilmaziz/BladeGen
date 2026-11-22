#include "occt/OcctBlade.hpp"

#include "occt/Fillet.hpp"

#include "math/Interpolation.hpp"
#include "math/Numerics.hpp"
#include "profile/Profile2D.hpp"

#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <Interface_Static.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <ShapeFix_Shape.hxx>
#include <StlAPI_Writer.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace BladeGen::Occt {

namespace {

gp_Pnt pnt(const Vec3& p) { return {p.x(), p.y(), p.z()}; }
Vec3 vec(const gp_Pnt& p) { return {p.X(), p.Y(), p.Z()}; }

Handle(Geom_BSplineCurve) interpolate_points(const std::vector<Vec3>& pts)
{
    Handle(TColgp_HArray1OfPnt) a = new TColgp_HArray1OfPnt(1, static_cast<int>(pts.size()));
    for (std::size_t i = 0; i < pts.size(); ++i) a->SetValue(static_cast<int>(i) + 1, pnt(pts[i]));
    GeomAPI_Interpolate I(a, Standard_False, 1e-12);
    I.Perform();
    if (!I.IsDone()) throw std::runtime_error("GeomAPI_Interpolate failed");
    return I.Curve();
}

double theta_of(const Vec3& p) { return std::atan2(p.y(), p.z()); }

/// m on a meridional curve where x(m) = x (extrapolating past the ends).
double m_at_x(const Blade::MeridionalCurve& c, double x)
{
    const double L = c.length();
    auto f = [&](double m) { return c.at_m(m).x() - x; };
    const double lo = -2.0 * L, hi = 3.0 * L;
    if ((f(lo) > 0.0) == (f(hi) > 0.0)) throw std::invalid_argument("meridional contour does not cover x = " + std::to_string(x));
    return Math::brent_root(f, lo, hi, 1e-14);
}

/// (x, r) of a contour between x_a and x_b, optionally offset inward by c along its normal.
std::vector<Vec2> contour(const Blade::MeridionalCurve& cv, double xa, double xb, double c, int n = 161)
{
    const double ma = m_at_x(cv, xa) - 2.0 * c, mb = m_at_x(cv, xb) + 2.0 * c;
    std::vector<Vec2> out;
    for (int i = 0; i < n; ++i) {
        const double m = ma + (mb - ma) * i / (n - 1);
        const Vec2 p = cv.at_m(m), t = cv.tangent_at_m(std::clamp(m, 0.0, cv.length()));
        out.emplace_back(p.x() + c * t.y(), p.y() - c * t.x());  // inward = -(-t_r, t_x)
    }
    return out;
}

std::vector<Vec3> in_plane(const std::vector<Vec2>& xr, double theta)
{
    std::vector<Vec3> out;
    for (const auto& q : xr) out.push_back(Blade::cartesian(q.x(), q.y(), theta));
    return out;
}

/// Solid of revolution of the meridional region between `lower` and `upper` (x, r)
/// polylines, from theta0 over dtheta (toward +theta).
TopoDS_Shape revolve_region(const std::vector<Vec2>& lower, const std::vector<Vec2>& upper, double theta0, double dtheta)
{
    const auto lo = in_plane(lower, theta0), up = in_plane(upper, theta0);
    BRepBuilderAPI_MakeWire w;
    w.Add(BRepBuilderAPI_MakeEdge(interpolate_points(lo)).Edge());
    w.Add(BRepBuilderAPI_MakeEdge(pnt(lo.back()), pnt(up.back())).Edge());
    std::vector<Vec3> up_rev(up.rbegin(), up.rend());
    w.Add(BRepBuilderAPI_MakeEdge(interpolate_points(up_rev)).Edge());
    w.Add(BRepBuilderAPI_MakeEdge(pnt(up.front()), pnt(lo.front())).Edge());
    if (!w.IsDone()) throw std::runtime_error("meridional wire failed");
    BRepBuilderAPI_MakeFace f(w.Wire(), Standard_True);
    if (!f.IsDone()) throw std::runtime_error("meridional face failed");
    // a rotation about -x by phi increases theta (y = r sin theta, z = r cos theta)
    BRepPrimAPI_MakeRevol rev(f.Face(), gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(-1, 0, 0)), dtheta);
    rev.Build();
    if (!rev.IsDone()) throw std::runtime_error("revolution failed");
    return rev.Shape();
}

TopoDS_Shape first_solid_containing(const TopoDS_Shape& s, const Vec3& p, bool want_in)
{
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) {
        BRepClass3d_SolidClassifier cl(ex.Current(), pnt(p), 1e-9);
        if ((cl.State() == TopAbs_IN) == want_in) return ex.Current();
    }
    return {};
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────

PeriodicInterpolator::PeriodicInterpolator(const std::vector<double>& u)
    : n_(static_cast<int>(u.size()) - 1), knots_(1, static_cast<int>(u.size())), mults_(1, static_cast<int>(u.size()))
{
    if (n_ < 4) throw std::invalid_argument("PeriodicInterpolator: need at least 4 points");
    for (int j = 0; j <= n_; ++j) {
        if (j > 0 && !(u[static_cast<std::size_t>(j)] > u[static_cast<std::size_t>(j - 1)]))
            throw std::invalid_argument("PeriodicInterpolator: parameters must increase");
        knots_.SetValue(j + 1, u[static_cast<std::size_t>(j)]);
        mults_.SetValue(j + 1, 1);  // all knots simple: C2 everywhere, also at the seam
    }
    TColgp_Array1OfPnt poles(1, n_);
    for (int i = 1; i <= n_; ++i) poles.SetValue(i, gp_Pnt(0, 0, 0));
    Handle(Geom_BSplineCurve) c = new Geom_BSplineCurve(poles, knots_, mults_, 3, Standard_True);
    // the interpolation matrix in OCCT's own periodic basis: A(j, i) = N_i(u_j)
    Eigen::MatrixXd A(n_, n_);
    for (int i = 1; i <= n_; ++i) {
        c->SetPole(i, gp_Pnt(1, 0, 0));
        for (int j = 0; j < n_; ++j) A(j, i - 1) = c->Value(u[static_cast<std::size_t>(j)]).X();
        c->SetPole(i, gp_Pnt(0, 0, 0));
    }
    lu_.compute(A);
    if (!(lu_.rcond() > 1e-13)) throw std::runtime_error("PeriodicInterpolator: singular or ill-conditioned system");
}

Handle(Geom_BSplineCurve) PeriodicInterpolator::curve(const std::vector<Vec3>& pts) const
{
    if (static_cast<int>(pts.size()) < n_) throw std::invalid_argument("PeriodicInterpolator: too few points");
    Eigen::MatrixXd B(n_, 3);
    for (int j = 0; j < n_; ++j) B.row(j) = pts[static_cast<std::size_t>(j)].transpose();
    const Eigen::MatrixXd X = lu_.solve(B);
    TColgp_Array1OfPnt poles(1, n_);
    for (int i = 0; i < n_; ++i) poles.SetValue(i + 1, gp_Pnt(X(i, 0), X(i, 1), X(i, 2)));
    return new Geom_BSplineCurve(poles, knots_, mults_, 3, Standard_True);
}

SkinResult skin_loops(const std::vector<std::vector<Vec3>>& loops, const std::vector<double>& u_given)
{
    const std::size_t K = loops.size();
    if (K < 2) throw std::invalid_argument("skin_loops: need at least 2 loops");
    const std::size_t N = loops[0].size() - 1;
    if (N < 6) throw std::invalid_argument("skin_loops: loops need at least 7 points");
    for (const auto& l : loops) {
        if (l.size() != N + 1) throw std::invalid_argument("skin_loops: every loop needs the same point count");
        if ((l.front() - l.back()).norm() > 1e-9 * (1.0 + l.front().norm()))
            throw std::invalid_argument("skin_loops: loops must be closed (last point = first)");
    }
    SkinResult res;
    // 1. common chordwise parameters: averaged normalised chord length
    res.u.assign(N + 1, 0.0);
    for (const auto& l : loops) {
        std::vector<double> s(N + 1, 0.0);
        for (std::size_t j = 1; j <= N; ++j) s[j] = s[j - 1] + (l[j] - l[j - 1]).norm();
        for (std::size_t j = 0; j <= N; ++j) res.u[j] += s[j] / s[N] / static_cast<double>(K);
    }
    res.u.front() = 0.0;
    res.u.back() = 1.0;
    if (!u_given.empty()) {
        if (u_given.size() != N + 1 || u_given.front() != 0.0 || u_given.back() != 1.0)
            throw std::invalid_argument("skin_loops: given u must have N + 1 values from 0 to 1");
        res.u = u_given;
    }
    const PeriodicInterpolator interp(res.u);
    std::vector<Handle(Geom_BSplineCurve)> curves;
    for (const auto& l : loops) curves.push_back(interp.curve(l));
    const int nu = curves[0]->NbPoles();
    for (const auto& c : curves) {
        if (c->NbPoles() != nu || c->NbKnots() != curves[0]->NbKnots())
            throw std::runtime_error("skin_loops: section curves are not compatible");
        for (int k = 1; k <= c->NbKnots(); ++k)
            if (std::abs(c->Knot(k) - curves[0]->Knot(k)) > 1e-14)
                throw std::runtime_error("skin_loops: section knots differ");
    }
    // 2. span parameters: averaged chord length between neighbouring loops
    res.v.assign(K, 0.0);
    for (std::size_t k = 1; k < K; ++k) {
        double d = 0.0;
        for (std::size_t j = 0; j < N; ++j) d += (loops[k][j] - loops[k - 1][j]).norm();
        res.v[k] = res.v[k - 1] + d / static_cast<double>(N);
    }
    for (auto& v : res.v) v /= res.v.back();
    // 3. spanwise C2 interpolation of every pole column
    std::vector<Math::BSpline3> cols;
    for (int i = 1; i <= nu; ++i) {
        std::vector<Vec3> q;
        for (const auto& c : curves) q.push_back(vec(c->Pole(i)));
        cols.push_back(Math::interpolate_cubic<3>(res.v, q));
    }
    const int nv = cols[0].n_ctrl();
    TColgp_Array2OfPnt poles(1, nu, 1, nv);
    for (int i = 1; i <= nu; ++i)
        for (int k = 1; k <= nv; ++k) poles.SetValue(i, k, pnt(cols[static_cast<std::size_t>(i - 1)].ctrl()[static_cast<std::size_t>(k - 1)]));
    TColStd_Array1OfReal uk(1, curves[0]->NbKnots());
    TColStd_Array1OfInteger um(1, curves[0]->NbKnots());
    curves[0]->Knots(uk);
    curves[0]->Multiplicities(um);
    std::vector<double> vk;
    std::vector<int> vm;
    for (double k : cols[0].knots()) {
        if (!vk.empty() && std::abs(k - vk.back()) < 1e-15) ++vm.back();
        else { vk.push_back(k); vm.push_back(1); }
    }
    TColStd_Array1OfReal vka(1, static_cast<int>(vk.size()));
    TColStd_Array1OfInteger vma(1, static_cast<int>(vk.size()));
    for (std::size_t i = 0; i < vk.size(); ++i) {
        vka.SetValue(static_cast<int>(i) + 1, vk[i]);
        vma.SetValue(static_cast<int>(i) + 1, vm[i]);
    }
    res.surface = new Geom_BSplineSurface(poles, uk, vka, um, vma, curves[0]->Degree(), cols[0].degree(),
                                          Standard_True, Standard_False);
    // 4. interpolation check
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t j = 0; j < N; j += 3)
            res.max_error = std::max(res.max_error, (vec(res.surface->Value(res.u[j], res.v[k])) - loops[k][j]).norm());
    return res;
}

std::vector<std::vector<Vec3>> extend_loops(const std::vector<std::vector<Vec3>>& loops, double before, double after)
{
    const auto sk = skin_loops(loops);
    const double dv = 1.0 / static_cast<double>(loops.size() - 1);
    auto extrap = [&](double v0, double h) {
        std::vector<Vec3> l;
        for (std::size_t j = 0; j + 1 < sk.u.size(); ++j) {
            gp_Pnt P;
            gp_Vec Du, Dv, Duu, Dvv, Duv;
            sk.surface->D2(sk.u[j], v0, P, Du, Dv, Duu, Dvv, Duv);
            l.push_back(vec(P) + h * Vec3(Dv.X(), Dv.Y(), Dv.Z()) + 0.5 * h * h * Vec3(Dvv.X(), Dvv.Y(), Dvv.Z()));
        }
        l.push_back(l.front());
        return l;
    };
    std::vector<std::vector<Vec3>> out;
    if (before > 0.0) out.push_back(extrap(0.0, -before * dv));
    out.insert(out.end(), loops.begin(), loops.end());
    if (after > 0.0) out.push_back(extrap(1.0, after * dv));
    return out;
}

Vec3 inside_point(const Blade::Section3D& s)
{
    const Profile::Profile2D p(s.params);
    const double xm = 0.45;
    const Vec2 a = p.suction_side().eval(p.ss_u_at_mprime(xm));
    const auto& ps = p.pressure_side();
    const double u = Math::first_root([&](double t) { return ps.eval(t).x() - xm; }, 0.0, 1.0, 256);
    const Vec2 b = std::isfinite(u) ? ps.eval(u) : Vec2(xm, 0.0);
    const auto map = Blade::make_section_map(s.streamline, s.placement);
    return map.cartesian(map.conformal(0.5 * (a + b)));
}

BladeShapes build_shapes(const std::vector<std::vector<Vec3>>& loops, const SolidOptions& opt, const Vec3& inside)
{
    BladeShapes out;
    out.skin = skin_loops(loops);
    BRepBuilderAPI_MakeFace mf(out.skin.surface, 1e-9);
    if (!mf.IsDone()) throw std::runtime_error("blade face failed");
    out.blade_face = mf.Face();

    // annulus sector around the blade
    double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin, tmin = xmin, tmax = -xmin;
    const double tref = theta_of(loops.front().front());  // angles relative to it: no wrap at +-pi
    for (const auto& l : loops)
        for (const auto& p : l) {
            xmin = std::min(xmin, p.x());
            xmax = std::max(xmax, p.x());
            const double th = tref + std::remainder(theta_of(p) - tref, 2.0 * Math::kPi);
            tmin = std::min(tmin, th);
            tmax = std::max(tmax, th);
        }
    if (tmax - tmin > 1.5 * Math::kPi) throw std::invalid_argument("blade spans more than 270 deg in theta");
    const double cx = xmax - xmin;
    const double xa = xmin - opt.axial_margin * cx, xb = xmax + opt.axial_margin * cx;
    const double marg = Math::deg2rad(opt.sector_margin_deg);
    const double theta0 = tmin - marg, dtheta = (tmax - tmin) + 2.0 * marg;
    const auto hub = contour(opt.hub, xa, xb, 0.0);
    const auto top = contour(opt.casing, xa, xb, opt.tip_clearance);
    const TopoDS_Shape sector = revolve_region(hub, top, theta0, dtheta);

    BRepAlgoAPI_Splitter sp;
    TopTools_ListOfShape args, tools;
    args.Append(sector);
    tools.Append(out.blade_face);
    sp.SetArguments(args);
    sp.SetTools(tools);
    sp.SetFuzzyValue(1e-9);
    sp.Build();
    if (!sp.IsDone() || sp.HasErrors()) throw std::runtime_error("splitting the annulus by the blade failed");
    out.blade = first_solid_containing(sp.Shape(), inside, true);
    if (out.blade.IsNull())
        throw std::runtime_error("blade solid not found: does the blade surface cross both hub and casing?");
    if (opt.blade_only) return out;
    out.blade_valid = is_valid(out.blade);
    out.blade_volume = volume(out.blade);
    out.blade_area = area(out.blade);

    // hub sector (bore to hub surface)
    double rmin = std::numeric_limits<double>::infinity();
    for (const auto& q : hub) rmin = std::min(rmin, q.y());
    const double rb = opt.hub_bore_radius > 0.0 ? opt.hub_bore_radius : 0.85 * rmin;
    std::vector<Vec2> bore = {{hub.front().x(), rb}, {hub.back().x(), rb}};
    std::vector<Vec2> bore_line;
    for (int i = 0; i <= 8; ++i) bore_line.push_back(bore[0] + (bore[1] - bore[0]) * i / 8.0);
    out.hub = revolve_region(bore_line, hub, theta0, dtheta);

    if (opt.fillet_radius > 0.0) {
        // BRepFilletAPI fails on this junction (walking failure; see Fillet.hpp), so
        // the rolling-ball construction is used directly
        try {
            FilletOptions fo;
            fo.radius = opt.fillet_radius;
            const auto fr = rolling_ball_fillet(out.blade, out.blade_face, out.hub, opt.hub, fo);
            out.fillet_ok = fr.ok;
            out.fillet_message = fr.message;
            out.fillet_g1_deg = fr.max_g1_deg;
            out.fillet_ball_error = fr.max_ball_error;
            out.fillet_face = fr.fillet_face;
            if (fr.ok) {
                out.blade_on_hub = fr.hub_blade;
                out.blade_with_fillet = fr.blade_fillet;
            }
        } catch (const std::exception& e) {
            out.fillet_ok = false;
            out.fillet_message = std::string("fillet failed: ") + e.what();
        } catch (const Standard_Failure& e) {
            out.fillet_ok = false;
            out.fillet_message = std::string("fillet failed (OCCT): ") + e.GetMessageString();
        }
    }
    if (out.blade_on_hub.IsNull()) {
        BRepAlgoAPI_Fuse fu(out.hub, out.blade);
        fu.Build();
        if (fu.IsDone()) out.blade_on_hub = fu.Shape();
    }
    return out;
}

BladeShapes build_blade(const Blade::Blade& blade, double fillet_radius, int n_edge, int n_body)
{
    const auto& def = blade.definition();
    auto secs = blade.sections(n_edge, n_body);
    std::vector<std::vector<Vec3>> loops;
    for (const auto& s : secs) loops.push_back(s.xyz);
    SolidOptions opt;
    const bool explicit_sl = std::all_of(def.stations.begin(), def.stations.end(),
                                         [](const Blade::BladeStation& s) { return s.streamline.has_value(); });
    if (explicit_sl) {
        opt.hub = *def.stations.front().streamline;
        opt.casing = *def.stations.back().streamline;
        loops = extend_loops(loops, 1.0, 1.0);
    } else {
        opt.hub = def.channel.hub();
        opt.casing = def.channel.shroud();
    }
    opt.tip_clearance = def.tip_clearance;
    opt.fillet_radius = fillet_radius;
    return build_shapes(loops, opt, inside_point(secs[secs.size() / 2]));
}

// ─────────────────────────────────────────────────────────────────────────────

TopoDS_Shape fillet_blade_hub(const TopoDS_Shape& blade, const TopoDS_Shape& hub, double radius, bool& ok,
                              std::string& message)
{
    ok = false;
    BRepAlgoAPI_Fuse fu(hub, blade);
    fu.Build();
    if (!fu.IsDone()) {
        message = "fusing blade and hub failed";
        return {};
    }
    fu.SimplifyResult();
    TopoDS_Shape fused = fu.Shape();
    // the junction: edges whose midpoint lies on the blade and on the hub
    TopTools_IndexedDataMapOfShapeListOfShape e2f;
    TopExp::MapShapesAndAncestors(fused, TopAbs_EDGE, TopAbs_FACE, e2f);
    BRepFilletAPI_MakeFillet mk(fused);
    int n = 0;
    for (int i = 1; i <= e2f.Extent(); ++i) {
        const TopoDS_Edge& e = TopoDS::Edge(e2f.FindKey(i));
        if (e2f.FindFromIndex(i).Extent() != 2) continue;
        double f, l;
        const Handle(Geom_Curve) c = BRep_Tool::Curve(e, f, l);
        if (c.IsNull()) continue;
        const gp_Pnt pm = c->Value(0.5 * (f + l));
        const TopoDS_Shape v = BRepBuilderAPI_MakeVertex(pm).Vertex();
        BRepExtrema_DistShapeShape db(v, blade), dh(v, hub);
        if (db.IsDone() && dh.IsDone() && db.Value() < 1e-6 && dh.Value() < 1e-6) {
            // not an edge of the blade cap lying inside the hub: both faces must not be caps
            mk.Add(radius, e);
            ++n;
        }
    }
    if (n == 0) {
        message = "no blade / hub junction edge found";
        return fused;
    }
    try {
        mk.Build();
    } catch (const Standard_Failure& f) {
        message = std::string("BRepFilletAPI failed: ") + f.GetMessageString();
        return fused;
    }
    if (!mk.IsDone()) {
        message = "BRepFilletAPI could not build the fillet";
        return fused;
    }
    TopoDS_Shape res = mk.Shape();
    if (!is_valid(res)) {
        ShapeFix_Shape fix(res);
        fix.Perform();
        res = fix.Shape();
    }
    ok = is_valid(res);
    message = ok ? "fillet r = " + std::to_string(radius) + " on " + std::to_string(n) + " edges"
                 : "fillet built but the result is not a valid solid";
    return ok ? res : fused;
}

// ─────────────────────────────────────────────────────────────────────────────

double volume(const TopoDS_Shape& s)
{
    GProp_GProps g;
    BRepGProp::VolumeProperties(s, g, 1e-7);
    return g.Mass();
}

double area(const TopoDS_Shape& s)
{
    GProp_GProps g;
    BRepGProp::SurfaceProperties(s, g, 1e-7);
    return g.Mass();
}

bool is_valid(const TopoDS_Shape& s) { return !s.IsNull() && BRepCheck_Analyzer(s).IsValid(); }

int count_faces(const TopoDS_Shape& s)
{
    int n = 0;
    for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next()) ++n;
    return n;
}

int count_solids(const TopoDS_Shape& s)
{
    int n = 0;
    for (TopExp_Explorer ex(s, TopAbs_SOLID); ex.More(); ex.Next()) ++n;
    return n;
}

void write_shape(const TopoDS_Shape& s, const std::string& path, Unit unit)
{
    if (s.IsNull()) throw std::invalid_argument("write_shape: null shape");
    std::string ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    const bool mm = unit == Unit::Millimetre;
    TopoDS_Shape out = s;
    if (mm) {
        gp_Trsf t;
        t.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        out = BRepBuilderAPI_Transform(s, t, Standard_True).Shape();
    }
    const char* u = mm ? "MM" : "M";
    if (ext == ".step" || ext == ".stp") {
        STEPControl_Writer w;
        Interface_Static::SetCVal("xstep.cascade.unit", u);
        Interface_Static::SetCVal("write.step.unit", u);
        Interface_Static::SetCVal("write.step.schema", "AP214IS");
        if (w.Transfer(out, STEPControl_AsIs) != IFSelect_RetDone) throw std::runtime_error("STEP transfer failed");
        if (w.Write(path.c_str()) != IFSelect_RetDone) throw std::runtime_error("STEP write failed: " + path);
    } else if (ext == ".igs" || ext == ".iges") {
        IGESControl_Controller::Init();
        Interface_Static::SetCVal("xstep.cascade.unit", u);
        IGESControl_Writer w(u, 1);
        if (!w.AddShape(out)) throw std::runtime_error("IGES transfer failed");
        w.ComputeModel();
        if (!w.Write(path.c_str())) throw std::runtime_error("IGES write failed: " + path);
    } else if (ext == ".stl") {
        BRepMesh_IncrementalMesh(out, mm ? 0.01 : 1e-5, Standard_False, 0.2, Standard_True);
        StlAPI_Writer w;
        w.ASCIIMode() = Standard_False;
        if (!w.Write(out, path.c_str())) throw std::runtime_error("STL write failed: " + path);
    } else if (ext == ".brep") {
        if (!BRepTools::Write(out, path.c_str())) throw std::runtime_error("BREP write failed: " + path);
    } else {
        throw std::invalid_argument("write_shape: unknown extension '" + ext + "' (.step .stp .igs .iges .stl .brep)");
    }
}

TopoDS_Shape read_step(const std::string& path)
{
    STEPControl_Reader r;
    Interface_Static::SetCVal("xstep.cascade.unit", "M");
    if (r.ReadFile(path.c_str()) != IFSelect_RetDone) throw std::runtime_error("cannot read STEP file " + path);
    r.TransferRoots();
    return r.OneShape();
}

Mesh tessellate(const TopoDS_Shape& s, double lin, double ang)
{
    BRepMesh_IncrementalMesh(s, lin, Standard_False, ang, Standard_True);
    Mesh m;
    int fid = 0;
    for (TopExp_Explorer ex(s, TopAbs_FACE); ex.More(); ex.Next(), ++fid) {
        const TopoDS_Face& f = TopoDS::Face(ex.Current());
        TopLoc_Location loc;
        const Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(f, loc);
        if (tri.IsNull()) continue;
        const gp_Trsf tr = loc.Transformation();
        const int base = static_cast<int>(m.vertices.size());
        for (int i = 1; i <= tri->NbNodes(); ++i) m.vertices.push_back(vec(tri->Node(i).Transformed(tr)));
        const bool rev = f.Orientation() == TopAbs_REVERSED;
        for (int i = 1; i <= tri->NbTriangles(); ++i) {
            int a, b, c;
            tri->Triangle(i).Get(a, b, c);
            if (rev) std::swap(b, c);
            m.triangles.push_back({base + a - 1, base + b - 1, base + c - 1});
            m.face_of_triangle.push_back(fid);
        }
    }
    return m;
}

double distance(const TopoDS_Shape& a, const TopoDS_Shape& b)
{
    BRepExtrema_DistShapeShape d(a, b);
    if (!d.IsDone()) throw std::runtime_error("distance computation failed");
    return d.Value();
}

std::vector<double> distances(const std::vector<Vec3>& pts, const TopoDS_Shape& s)
{
    std::vector<double> out;
    out.reserve(pts.size());
    for (const auto& p : pts) out.push_back(distance(BRepBuilderAPI_MakeVertex(pnt(p)).Vertex(), s));
    return out;
}

}  // namespace BladeGen::Occt
