// tests/turbolenza/GeneralFixture.hpp — meshes and manufactured solutions for the general operator tests
// (test_general_*.cpp): single blocks from analytic maps of the unit square, the three-block periodic
// H-domain with turned blocks, the O-block round a body, the manufactured fields and the solve.
#pragma once

#include "turbolenza/GeneralOperator.hpp"
#include "turbolenza/LinearSolver.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace gfix {

namespace T = BladeGen::Turbolenza;

inline constexpr double kPi = 3.141592653589793;

inline double measure(const std::string& name, double value)
{
    std::printf("MEASURE %-60s %.3e\n", name.c_str(), value);
    return value;
}

/// (u, v) in [0, 1]^2 -> (x, y).
using Map = std::function<std::array<double, 2>(double, double)>;

/// A block of map(i / (ni - 1), j / (nj - 1)); m = x, r = r(x, y), b = b(x, y) (default 1).
inline T::MeshBlock block(const Map& map, int ni, int nj, const std::function<double(double, double)>& r = nullptr,
                          const std::function<double(double, double)>& b = nullptr, const std::string& name = "block")
{
    T::MeshBlock B;
    B.name = name;
    B.ni = ni;
    B.nj = nj;
    const std::size_t n = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
    B.x.resize(n);
    B.y.resize(n);
    B.m.resize(n);
    B.r.resize(n);
    B.b.resize(n);
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni; ++i) {
            const auto p = map(static_cast<double>(i) / (ni - 1), static_cast<double>(j) / (nj - 1));
            const std::size_t k = B.at(i, j);
            B.x[k] = p[0];
            B.y[k] = p[1];
            B.m[k] = p[0];
            B.r[k] = r ? r(p[0], p[1]) : 1.0;
            B.b[k] = b ? b(p[0], p[1]) : 1.0;
        }
    return B;
}

inline T::FaceSegment segment(int blk, T::BlockFace face, int begin, int end, T::SegmentKind kind)
{
    T::FaceSegment s;
    s.block = blk;
    s.face = face;
    s.begin = begin;
    s.end = end;
    s.kind = kind;
    return s;
}

/// The whole face as one segment.
inline T::FaceSegment whole(const T::GeneralMesh& mesh, int blk, T::BlockFace face, T::SegmentKind kind)
{
    return segment(blk, face, 0, T::face_length(mesh.blocks[static_cast<std::size_t>(blk)], face) - 1, kind);
}

/// One block with Dirichlet segments on its four faces.
inline T::GeneralMesh dirichlet_box(T::MeshBlock B)
{
    T::GeneralMesh mesh;
    mesh.pitch = 1.0;
    mesh.blocks.push_back(std::move(B));
    for (T::BlockFace f : {T::BlockFace::IMin, T::BlockFace::IMax, T::BlockFace::JMin, T::BlockFace::JMax})
        mesh.segments.push_back(whole(mesh, 0, f, T::SegmentKind::Dirichlet));
    return mesh;
}

/// Connects segments a and b (image as given), finding `reversed` from the coordinates.
inline void connect(T::GeneralMesh& mesh, int a, int b, int image)
{
    const T::FaceSegment& A = mesh.segments[static_cast<std::size_t>(a)];
    const T::FaceSegment& B = mesh.segments[static_cast<std::size_t>(b)];
    const T::MeshBlock& BA = mesh.blocks[static_cast<std::size_t>(A.block)];
    const T::MeshBlock& BB = mesh.blocks[static_cast<std::size_t>(B.block)];
    const auto pa = T::face_node(BA, A.face, A.begin);
    const auto pb = T::face_node(BB, B.face, B.begin);
    const std::size_t qa = BA.at(pa[0], pa[1]), qb = BB.at(pb[0], pb[1]);
    const bool same = std::abs(BA.x[qa] - BB.x[qb]) < 1e-9 && std::abs(BB.y[qb] - BA.y[qa] - image * mesh.pitch) < 1e-9;
    mesh.connections.push_back({a, b, !same, image});
}

/// A manufactured field: psi, its derivatives, the coefficient C = mdot' kappa and its derivatives.
struct Manufactured
{
    std::function<double(double, double)> psi, psi_x, psi_y, psi_xx, psi_yy;
    std::function<double(double, double)> c, c_x, c_y;
    /// div(C grad psi)
    [[nodiscard]] double flux_divergence(double x, double y) const
    {
        return c(x, y) * (psi_xx(x, y) + psi_yy(x, y)) + c_x(x, y) * psi_x(x, y) + c_y(x, y) * psi_y(x, y);
    }
};

/// psi = exp(0.7 x) sin(1.1 y) + 0.3 x^2 y; C = 1 + 0.25 sin(x + 2 y) (or 1 when constant).
inline Manufactured smooth(bool constant_coefficient)
{
    Manufactured M;
    M.psi = [](double x, double y) { return std::exp(0.7 * x) * std::sin(1.1 * y) + 0.3 * x * x * y; };
    M.psi_x = [](double x, double y) { return 0.7 * std::exp(0.7 * x) * std::sin(1.1 * y) + 0.6 * x * y; };
    M.psi_y = [](double x, double y) { return 1.1 * std::exp(0.7 * x) * std::cos(1.1 * y) + 0.3 * x * x; };
    M.psi_xx = [](double x, double y) { return 0.49 * std::exp(0.7 * x) * std::sin(1.1 * y) + 0.6 * y; };
    M.psi_yy = [](double x, double y) { return -1.21 * std::exp(0.7 * x) * std::sin(1.1 * y); };
    if (constant_coefficient) {
        M.c = [](double, double) { return 1.0; };
        M.c_x = M.c_y = [](double, double) { return 0.0; };
    } else {
        M.c = [](double x, double y) { return 1.0 + 0.25 * std::sin(x + 2.0 * y); };
        M.c_x = [](double x, double y) { return 0.25 * std::cos(x + 2.0 * y); };
        M.c_y = [](double x, double y) { return 0.5 * std::cos(x + 2.0 * y); };
    }
    return M;
}

/// Quasi-periodic in y with period P: psi = y / P + 0.05 exp(0.5 x) sin(k y) + 0.2 x^2 + 0.1 x, k = 2 pi / P;
/// C = 1 + 0.25 sin(x + k y).
inline Manufactured quasi_periodic(double P)
{
    const double k = 2.0 * kPi / P;
    Manufactured M;
    M.psi = [=](double x, double y) { return y / P + 0.05 * std::exp(0.5 * x) * std::sin(k * y) + 0.2 * x * x + 0.1 * x; };
    M.psi_x = [=](double x, double y) { return 0.025 * std::exp(0.5 * x) * std::sin(k * y) + 0.4 * x + 0.1; };
    M.psi_y = [=](double x, double y) { return 1.0 / P + 0.05 * k * std::exp(0.5 * x) * std::cos(k * y); };
    M.psi_xx = [=](double x, double y) { return 0.0125 * std::exp(0.5 * x) * std::sin(k * y) + 0.4; };
    M.psi_yy = [=](double x, double y) { return -0.05 * k * k * std::exp(0.5 * x) * std::sin(k * y); };
    M.c = [=](double x, double y) { return 1.0 + 0.25 * std::sin(x + k * y); };
    M.c_x = [=](double x, double y) { return 0.25 * std::cos(x + k * y); };
    M.c_y = [=](double x, double y) { return 0.25 * k * std::cos(x + k * y); };
    return M;
}

struct MmsError
{
    double psi = 0.0, grad = 0.0;  ///< max |psi - exact|, max |(psi_x, psi_y) - exact|
    int unknowns = 0;
    double h = 0.0;                ///< 1 / cells of the first block along i
};

/// The fields of a manufactured solve: coefficient C, source J div(C grad psi) (+ extra), Dirichlet exact,
/// flow-angle tangent -psi_x / psi_y.
inline T::OperatorFields mms_fields(const T::GeneralOperator& op, const Manufactured& M,
                                    const std::vector<double>* extra_source = nullptr)
{
    T::OperatorFields f;
    const int n = op.n_nodes();
    f.coefficient.resize(static_cast<std::size_t>(n));
    f.source.resize(static_cast<std::size_t>(n));
    f.flow_angle_tan.resize(static_cast<std::size_t>(n));
    for (int g = 0; g < n; ++g) {
        const double x = op.x(g), y = op.y(g);
        f.coefficient[static_cast<std::size_t>(g)] = M.c(x, y);
        f.source[static_cast<std::size_t>(g)] =
            op.metric(g).J * M.flux_divergence(x, y) + (extra_source ? (*extra_source)[static_cast<std::size_t>(g)] : 0.0);
        f.flow_angle_tan[static_cast<std::size_t>(g)] = -M.psi_x(x, y) / M.psi_y(x, y);
    }
    f.dirichlet = [&op, M](int g) { return M.psi(op.x(g), op.y(g)); };
    return f;
}

/// Solves the manufactured problem and returns psi at every node.
inline std::vector<double> mms_solve(const T::GeneralOperator& op, const T::OperatorFields& f)
{
    const T::GeneralSystem sys = op.assemble(f);
    return op.expand(T::solve(sys.matrix, sys.rhs), f);
}

inline MmsError mms_error(const T::GeneralOperator& op, const Manufactured& M, const std::vector<double>& psi)
{
    MmsError e;
    const T::StreamGradient G = op.gradient(psi);
    for (int g = 0; g < op.n_nodes(); ++g) {
        const double x = op.x(g), y = op.y(g);
        e.psi = std::max(e.psi, std::abs(psi[static_cast<std::size_t>(g)] - M.psi(x, y)));
        e.grad = std::max(e.grad, std::hypot(G.psi_x[static_cast<std::size_t>(g)] - M.psi_x(x, y),
                                             G.psi_y[static_cast<std::size_t>(g)] - M.psi_y(x, y)));
    }
    e.unknowns = op.n_unknowns();
    e.h = 1.0 / (op.mesh().blocks[0].ni - 1);
    return e;
}

// ── the mesh families of the verification plan (maps of the unit square) ─────────────────────────────
inline double stretch(double u) { return (std::exp(2.0 * u) - 1.0) / (std::exp(2.0) - 1.0); }
inline double stretch_v(double v) { return v + 0.3 * v * (1.0 - v); }

inline Map family(int k)
{
    switch (k) {
    case 1: return [](double u, double v) { return std::array<double, 2>{u, 0.8 * v}; };            // uniform orthogonal
    case 2: return [](double u, double v) { return std::array<double, 2>{u, 0.8 * v + 0.6 * u}; };  // uniform sheared
    case 3: return [](double u, double v) { return std::array<double, 2>{stretch(u), 0.8 * stretch_v(v)}; };
    case 4:
        return [](double u, double v) {  // stretched, non-orthogonal (Aungier-like y0(x) + eta S(x))
            const double X = stretch(u);
            return std::array<double, 2>{X, 0.8 * stretch_v(v) * (1.0 + 0.3 * X) + 0.4 * std::sin(1.5 * X)};
        };
    default:
        return [](double u, double v) {  // curvilinear: curved lines of both families
            const double s = 0.05 * std::sin(2.0 * kPi * u) * std::sin(2.0 * kPi * v);
            return std::array<double, 2>{u + s, 0.8 * (v + s) + 0.2 * u * u};
        };
    }
}

inline const char* family_name(int k)
{
    switch (k) {
    case 1: return "uniform orthogonal";
    case 2: return "uniform sheared";
    case 3: return "stretched orthogonal";
    case 4: return "stretched non-orthogonal";
    default: return "curvilinear";
    }
}

/// The periodic H-domain x in [0, 1.5] (curved stations, curved periodic lines, pitch P) in three blocks
/// of M x M cells: block 0 as is, block 1 with i and j reversed, block 2 turned by 90 degrees (i'' = j,
/// j'' = n_i - 1 - i).  Inlet (block 0 imin) Dirichlet or an Inlet segment; outlet (block 2) Dirichlet
/// or an Outlet segment.
inline T::GeneralMesh three_block_periodic(int M, double P, T::SegmentKind inlet, T::SegmentKind outlet)
{
    auto map = [P](double U, double v) {  // U in [0, 1.5]
        const double s = std::sin(2.0 * kPi * v) * std::sin(2.0 * kPi * U / 1.5);
        return std::array<double, 2>{U + 0.03 * s, P * v + 0.3 * std::sin(2.0 * U) + 0.04 * s};
    };
    T::GeneralMesh mesh;
    mesh.pitch = P;
    const int n = M + 1;
    // block 0: U in [0, 0.5]
    mesh.blocks.push_back(block([&](double u, double v) { return map(0.5 * u, v); }, n, n, nullptr, nullptr, "b0"));
    // block 1: U in [0.5, 1], i' = n - 1 - i, j' = n - 1 - j
    mesh.blocks.push_back(block([&](double u, double v) { return map(0.5 + 0.5 * (1.0 - u), 1.0 - v); }, n, n, nullptr, nullptr, "b1"));
    // block 2: U in [1, 1.5], i'' = j (v), j'' = n - 1 - i (U decreasing)
    mesh.blocks.push_back(block([&](double u, double v) { return map(1.0 + 0.5 * (1.0 - v), u); }, n, n, nullptr, nullptr, "b2"));
    using F = T::BlockFace;
    using K = T::SegmentKind;
    mesh.segments = {
        whole(mesh, 0, F::IMin, inlet),         // 0
        whole(mesh, 0, F::IMax, K::Interface),  // 1  U = 0.5
        whole(mesh, 0, F::JMin, K::Periodic),   // 2  v = 0
        whole(mesh, 0, F::JMax, K::Periodic),   // 3  v = 1
        whole(mesh, 1, F::IMax, K::Interface),  // 4  U = 0.5
        whole(mesh, 1, F::IMin, K::Interface),  // 5  U = 1
        whole(mesh, 1, F::JMin, K::Periodic),   // 6  v = 1
        whole(mesh, 1, F::JMax, K::Periodic),   // 7  v = 0
        whole(mesh, 2, F::JMax, K::Interface),  // 8  U = 1
        whole(mesh, 2, F::JMin, outlet),        // 9  U = 1.5
        whole(mesh, 2, F::IMin, K::Periodic),   // 10 v = 0
        whole(mesh, 2, F::IMax, K::Periodic),   // 11 v = 1
    };
    connect(mesh, 2, 3, 1);
    connect(mesh, 1, 4, 0);
    connect(mesh, 6, 7, -1);
    connect(mesh, 5, 8, 0);
    connect(mesh, 10, 11, 1);
    return mesh;
}

/// The O-block round the ellipse x^2 + (y / 0.6)^2 = 1: i clockwise (theta = -2 pi i / (n_i - 1)), j
/// outward to the ellipse of three times the size, stretched; seam i = 0 / n_i - 1 an interface of the
/// block with itself; Dirichlet on the body and the outer boundary.
inline T::GeneralMesh oblock(int ni, int nj)
{
    T::GeneralMesh mesh;
    mesh.pitch = 1.0;
    mesh.blocks.push_back(block(
        [](double u, double v) {
            const double phi = -2.0 * kPi * u, R = std::exp(std::log(3.0) * v);
            return std::array<double, 2>{R * std::cos(phi) + 0.1 * v * v, 0.6 * R * std::sin(phi)};
        },
        ni, nj, nullptr, nullptr, "o"));
    using F = T::BlockFace;
    using K = T::SegmentKind;
    mesh.segments = {whole(mesh, 0, F::IMin, K::Interface), whole(mesh, 0, F::IMax, K::Interface),
                     whole(mesh, 0, F::JMin, K::Dirichlet), whole(mesh, 0, F::JMax, K::Dirichlet)};
    mesh.connections.push_back({0, 1, false, 0});
    // the seam nodes coincide exactly: copy them (cos(-2 pi) is not exactly 1)
    T::MeshBlock& B = mesh.blocks[0];
    for (int j = 0; j < nj; ++j) {
        B.x[B.at(ni - 1, j)] = B.x[B.at(0, j)];
        B.y[B.at(ni - 1, j)] = B.y[B.at(0, j)];
        B.m[B.at(ni - 1, j)] = B.m[B.at(0, j)];
    }
    return mesh;
}

inline double order(double coarse, double fine) { return std::log2(coarse / fine); }

}  // namespace gfix

namespace gfix {

/// The edge shapes of hgrid_with_blade: 0 zero thickness at the edge stations (h = 0, walls on the
/// periodic lines), 1 sharp (h = t sin(pi (x - 1))), 2 round (h = t sqrt(sin)), 3 a C-infinity bump.
inline const char* edge_name(int kind)
{
    switch (kind) {
    case 0: return "zero thickness";
    case 1: return "sharp edge";
    case 2: return "round edge";
    default: return "C-infinity bump";
    }
}

/// Aungier's H-grid layout on a manufactured problem (the layout of the independent review of 0ede48c):
/// x in [0, 3], stations x_i = (i + 1/2) / n (none on an edge), pitch P, a blade on 1 < x < 2 with the
/// suction side y = 0.2 x + h(x) on j = 0 and the pressure side of the +pitch neighbour P + 0.2 x - h(x) on
/// j = n_j - 1 (h by `kind`, edge_name), periodic side lines outside it; walls, inlet and outlet Dirichlet.
inline T::GeneralMesh hgrid_with_blade(int n, int nj, double P, double t, int kind)
{
    const int ni = 3 * n;
    const double dx = 1.0 / n;
    auto h = [t, kind](double x) {
        if (!(x > 1.0 && x < 2.0) || kind == 0) return 0.0;
        const double s = std::sin(kPi * (x - 1.0));
        if (kind == 1) return t * s;
        if (kind == 2) return t * std::sqrt(s);
        const double q = 2.0 * (x - 1.5);
        return t * std::exp(1.0 - 1.0 / (1.0 - q * q));
    };
    T::GeneralMesh mesh;
    mesh.pitch = P;
    T::MeshBlock B;
    B.name = "h";
    B.ni = ni;
    B.nj = nj;
    const std::size_t nn = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
    B.x.resize(nn);
    B.y.resize(nn);
    B.m.resize(nn);
    B.r.assign(nn, 1.0);
    B.b.assign(nn, 1.0);
    for (int i = 0; i < ni; ++i) {
        const double x = (i + 0.5) * dx;
        for (int j = 0; j < nj; ++j) {
            const double eta = static_cast<double>(j) / (nj - 1);
            const double y0 = 0.2 * x + h(x), y1 = P + 0.2 * x - h(x);
            const std::size_t k = B.at(i, j);
            B.x[k] = x;
            B.y[k] = y0 + eta * (y1 - y0);
            B.m[k] = x;
        }
    }
    mesh.blocks.push_back(std::move(B));
    const int f = n, l = 2 * n - 1;  // the first and last station inside 1 < x < 2
    using F = T::BlockFace;
    using K = T::SegmentKind;
    mesh.segments = {segment(0, F::JMin, 0, f - 1, K::Periodic), segment(0, F::JMin, f, l, K::Dirichlet),
                     segment(0, F::JMin, l + 1, ni - 1, K::Periodic), segment(0, F::JMax, 0, f - 1, K::Periodic),
                     segment(0, F::JMax, f, l, K::Dirichlet), segment(0, F::JMax, l + 1, ni - 1, K::Periodic),
                     segment(0, F::IMin, 0, nj - 1, K::Dirichlet), segment(0, F::IMax, 0, nj - 1, K::Dirichlet)};
    mesh.connections = {{0, 3, false, 1}, {2, 5, false, 1}};
    return mesh;
}

}  // namespace gfix
