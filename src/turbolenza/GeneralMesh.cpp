// The general multi-block mesh (see GeneralMesh.hpp): names, face helpers, the connection index map
// and the validation of every rule of the header.  Not a port of v1.
#include "turbolenza/GeneralMesh.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <map>
#include <tuple>
#include <utility>

namespace BladeGen::Turbolenza {

const char* block_face_name(BlockFace f) noexcept
{
    switch (f) {
    case BlockFace::IMin: return "imin";
    case BlockFace::IMax: return "imax";
    case BlockFace::JMin: return "jmin";
    case BlockFace::JMax: return "jmax";
    }
    return "?";
}

const char* segment_kind_name(SegmentKind k) noexcept
{
    switch (k) {
    case SegmentKind::Blade: return "blade";
    case SegmentKind::Inlet: return "inlet";
    case SegmentKind::Outlet: return "outlet";
    case SegmentKind::Periodic: return "periodic";
    case SegmentKind::Interface: return "interface";
    case SegmentKind::Dirichlet: return "dirichlet";
    }
    return "?";
}

const char* edge_kind_name(EdgeKind e) noexcept { return e == EdgeKind::Leading ? "leading" : "trailing"; }
const char* mesh_surface_name(Surface s) noexcept { return s == Surface::Suction ? "suction" : "pressure"; }

int face_length(const MeshBlock& block, BlockFace face) noexcept
{
    return (face == BlockFace::IMin || face == BlockFace::IMax) ? block.nj : block.ni;
}

std::array<int, 2> face_node(const MeshBlock& block, BlockFace face, int t) noexcept
{
    switch (face) {
    case BlockFace::IMin: return {0, t};
    case BlockFace::IMax: return {block.ni - 1, t};
    case BlockFace::JMin: return {t, 0};
    case BlockFace::JMax: return {t, block.nj - 1};
    }
    return {0, 0};
}

namespace {

std::array<int, 2> outward(BlockFace f) noexcept
{
    switch (f) {
    case BlockFace::IMin: return {-1, 0};
    case BlockFace::IMax: return {1, 0};
    case BlockFace::JMin: return {0, -1};
    case BlockFace::JMax: return {0, 1};
    }
    return {0, 0};
}

std::array<int, 2> tangent(BlockFace f) noexcept
{
    return (f == BlockFace::IMin || f == BlockFace::IMax) ? std::array<int, 2>{0, 1} : std::array<int, 2>{1, 0};
}

double cross(double ax, double ay, double bx, double by) noexcept { return ax * by - ay * bx; }

std::string seg_name(const GeneralMesh& mesh, int s)
{
    const FaceSegment& g = mesh.segments[static_cast<std::size_t>(s)];
    return std::format("segment {} ({} of block {} {} [{}, {}])", s, segment_kind_name(g.kind), g.block,
                       block_face_name(g.face), g.begin, g.end);
}

bool close(double a, double b, double tol) noexcept { return std::abs(a - b) <= tol; }
bool close_rel(double a, double b, double tol) noexcept
{
    return std::abs(a - b) <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

void check_block(const MeshBlock& B, int k)
{
    if (B.ni < 3 || B.nj < 3)
        throw MeshError(std::format("block {} ('{}') has {} x {} nodes: at least 3 x 3 are needed (one-sided second-order "
                                    "differences)",
                                    k, B.name, B.ni, B.nj));
    const std::size_t n = static_cast<std::size_t>(B.ni) * static_cast<std::size_t>(B.nj);
    const std::pair<const std::vector<double>*, const char*> fields[] = {
        {&B.x, "m'"}, {&B.y, "theta"}, {&B.m, "m"}, {&B.r, "r"}, {&B.b, "b"}};
    for (const auto& [f, name] : fields) {
        if (f->size() != n)
            throw MeshError(std::format("block {} ('{}'): {} has {} values, the block {} x {} = {} nodes", k, B.name, name,
                                        f->size(), B.ni, B.nj, n));
        for (std::size_t q = 0; q < n; ++q)
            if (!std::isfinite((*f)[q]))
                throw MeshError(std::format("block {} ('{}'): {} is not finite at node {}", k, B.name, name, q));
    }
    for (std::size_t q = 0; q < n; ++q) {
        if (!(B.r[q] > 0.0)) throw MeshError(std::format("block {} ('{}'): r = {} at node {}: the radius must be positive", k, B.name, B.r[q], q));
        if (!(B.b[q] > 0.0)) throw MeshError(std::format("block {} ('{}'): b = {} at node {}: the thickness must be positive", k, B.name, B.b[q], q));
    }
    // every corner Jacobian of every cell positive (right-handed, ENIGMAgrid ENG-334)
    for (int j = 0; j + 1 < B.nj; ++j)
        for (int i = 0; i + 1 < B.ni; ++i) {
            const std::size_t a = B.at(i, j), bb = B.at(i + 1, j), c = B.at(i, j + 1), d = B.at(i + 1, j + 1);
            const double ix0 = B.x[bb] - B.x[a], iy0 = B.y[bb] - B.y[a];  // i edge at j
            const double ix1 = B.x[d] - B.x[c], iy1 = B.y[d] - B.y[c];    // i edge at j + 1
            const double jx0 = B.x[c] - B.x[a], jy0 = B.y[c] - B.y[a];    // j edge at i
            const double jx1 = B.x[d] - B.x[bb], jy1 = B.y[d] - B.y[bb];  // j edge at i + 1
            const double corners[4] = {cross(ix0, iy0, jx0, jy0), cross(ix0, iy0, jx1, jy1), cross(ix1, iy1, jx0, jy0),
                                       cross(ix1, iy1, jx1, jy1)};
            for (double J : corners)
                if (!(J > 0.0))
                    throw MeshError(std::format("block {} ('{}'): the cell ({}, {}) has a corner Jacobian {:.3e} <= 0: the "
                                                "block must be right-handed in (m', theta) with no folded cell (i along "
                                                "+m', j along +theta; ENIGMAgrid ENG-334)",
                                                k, B.name, i, j, J));
        }
}

}  // namespace

std::array<int, 4> face_transform(BlockFace from, BlockFace to, bool reversed) noexcept
{
    // T [u_from tau_from] = [-u_to  +-tau_to]; [u_from tau_from] is a signed permutation, so its
    // inverse is its transpose: T = N M^T.
    const auto uf = outward(from), tf = tangent(from), ut = outward(to), tt = tangent(to);
    const int s = reversed ? -1 : 1;
    const int N[2][2] = {{-ut[0], s * tt[0]}, {-ut[1], s * tt[1]}};  // columns: images of u_from, tau_from
    const int M[2][2] = {{uf[0], tf[0]}, {uf[1], tf[1]}};
    std::array<int, 4> T{};
    for (int r = 0; r < 2; ++r)
        for (int c = 0; c < 2; ++c) T[static_cast<std::size_t>(2 * r + c)] = N[r][0] * M[c][0] + N[r][1] * M[c][1];
    return T;
}

void validate_mesh(const GeneralMesh& mesh, const MeshCheckOptions& options)
{
    if (!(std::isfinite(mesh.pitch) && mesh.pitch > 0.0))
        throw MeshError(std::format("the pitch is {}: a positive shift in theta [rad] is needed", mesh.pitch));
    if (mesh.blocks.empty()) throw MeshError("the mesh has no block");
    for (std::size_t k = 0; k < mesh.blocks.size(); ++k) check_block(mesh.blocks[k], static_cast<int>(k));

    const int nb = static_cast<int>(mesh.blocks.size());
    // segments: on their face, ordered, no shared cell, every face node covered
    std::map<std::pair<int, int>, std::vector<int>> on_face;
    for (std::size_t s = 0; s < mesh.segments.size(); ++s) {
        const FaceSegment& g = mesh.segments[s];
        if (g.block < 0 || g.block >= nb)
            throw MeshError(std::format("segment {} names block {}, the mesh has {}", s, g.block, nb));
        const int len = face_length(mesh.blocks[static_cast<std::size_t>(g.block)], g.face);
        if (!(0 <= g.begin && g.begin < g.end && g.end < len))
            throw MeshError(std::format("{}: the range must satisfy 0 <= begin < end < {} (the face's nodes)",
                                        seg_name(mesh, static_cast<int>(s)), len));
        on_face[{g.block, static_cast<int>(g.face)}].push_back(static_cast<int>(s));
    }
    for (int k = 0; k < nb; ++k)
        for (int f = 0; f < 4; ++f) {
            const BlockFace face = static_cast<BlockFace>(f);
            const int len = face_length(mesh.blocks[static_cast<std::size_t>(k)], face);
            std::vector<int> list = on_face[{k, f}];
            std::sort(list.begin(), list.end(), [&](int a, int b) {
                return mesh.segments[static_cast<std::size_t>(a)].begin < mesh.segments[static_cast<std::size_t>(b)].begin;
            });
            std::vector<char> covered(static_cast<std::size_t>(len), 0);
            for (std::size_t q = 0; q < list.size(); ++q) {
                const FaceSegment& g = mesh.segments[static_cast<std::size_t>(list[q])];
                if (q > 0) {
                    const FaceSegment& p = mesh.segments[static_cast<std::size_t>(list[q - 1])];
                    if (g.begin < p.end)
                        throw MeshError(std::format("{} and {} share a cell of the face", seg_name(mesh, list[q - 1]),
                                                    seg_name(mesh, list[q])));
                }
                for (int t = g.begin; t <= g.end; ++t) covered[static_cast<std::size_t>(t)] = 1;
            }
            for (int t = 0; t < len; ++t)
                if (!covered[static_cast<std::size_t>(t)])
                    throw MeshError(std::format("node {} of face {} of block {} lies on no segment: every face node needs "
                                                "a boundary kind or a connection",
                                                t, block_face_name(face), k));
        }

    // connections: every periodic and interface segment exactly once, matching, coinciding, orientation kept
    std::vector<int> used(mesh.segments.size(), 0);
    for (std::size_t c = 0; c < mesh.connections.size(); ++c) {
        const Connection& cn = mesh.connections[c];
        const int ns = static_cast<int>(mesh.segments.size());
        if (cn.first < 0 || cn.first >= ns || cn.second < 0 || cn.second >= ns)
            throw MeshError(std::format("connection {} names segments {} and {}, the mesh has {}", c, cn.first, cn.second, ns));
        if (cn.first == cn.second)
            throw MeshError(std::format("connection {} pairs {} with itself", c, seg_name(mesh, cn.first)));
        const FaceSegment& a = mesh.segments[static_cast<std::size_t>(cn.first)];
        const FaceSegment& b = mesh.segments[static_cast<std::size_t>(cn.second)];
        const bool periodic = a.kind == SegmentKind::Periodic;
        if (!((a.kind == SegmentKind::Periodic || a.kind == SegmentKind::Interface) && b.kind == a.kind))
            throw MeshError(std::format("connection {} pairs {} with {}: a connection joins two periodic or two interface "
                                        "segments",
                                        c, seg_name(mesh, cn.first), seg_name(mesh, cn.second)));
        if (periodic && !(cn.image == 1 || cn.image == -1))
            throw MeshError(std::format("periodic connection {} has image {}: +1 or -1 (one pitch)", c, cn.image));
        if (!periodic && cn.image != 0)
            throw MeshError(std::format("interface connection {} has image {}: an interface joins coinciding nodes (0)", c,
                                        cn.image));
        if (a.length() != b.length())
            throw MeshError(std::format("connection {}: {} has {} nodes, {} has {}: nonmatching interfaces are not "
                                        "supported",
                                        c, seg_name(mesh, cn.first), a.length(), seg_name(mesh, cn.second), b.length()));
        const std::array<int, 4> T = face_transform(a.face, b.face, cn.reversed);
        if (T[0] * T[3] - T[1] * T[2] != 1)
            throw MeshError(std::format("connection {} ({} to {}{}) reverses the orientation: one of the two blocks would "
                                        "fold over the other (flip 'reversed')",
                                        c, seg_name(mesh, cn.first), seg_name(mesh, cn.second),
                                        cn.reversed ? ", reversed" : ""));
        ++used[static_cast<std::size_t>(cn.first)];
        ++used[static_cast<std::size_t>(cn.second)];
        const MeshBlock& A = mesh.blocks[static_cast<std::size_t>(a.block)];
        const MeshBlock& B = mesh.blocks[static_cast<std::size_t>(b.block)];
        const double tol = options.coincidence;
        for (int k = 0; k < a.length(); ++k) {
            const auto pa = face_node(A, a.face, a.begin + k);
            const auto pb = face_node(B, b.face, cn.reversed ? b.end - k : b.begin + k);
            const std::size_t qa = A.at(pa[0], pa[1]), qb = B.at(pb[0], pb[1]);
            if (a.block == b.block && qa == qb)
                throw MeshError(std::format("connection {} pairs node ({}, {}) of block {} with itself", c, pa[0], pa[1], a.block));
            const double dy = B.y[qb] - A.y[qa] - static_cast<double>(cn.image) * mesh.pitch;
            if (!close(A.x[qa], B.x[qb], tol) || !(std::abs(dy) <= tol) || !close_rel(A.m[qa], B.m[qb], tol) ||
                !close_rel(A.r[qa], B.r[qb], tol) || !close_rel(A.b[qa], B.b[qb], tol))
                throw MeshError(std::format("connection {}: node ({}, {}) of block {} and node ({}, {}) of block {} do not "
                                            "coincide{} (dm' = {:.3e}, dtheta = {:.3e}; tolerance {:g}): a nonmatching "
                                            "interface or a wrong pairing",
                                            c, pa[0], pa[1], a.block, pb[0], pb[1], b.block,
                                            periodic ? std::format(" after a shift of {} pitch", cn.image) : "",
                                            B.x[qb] - A.x[qa], dy, tol));
        }
    }
    for (std::size_t s = 0; s < mesh.segments.size(); ++s) {
        const SegmentKind k = mesh.segments[s].kind;
        if ((k == SegmentKind::Periodic || k == SegmentKind::Interface) && used[s] != 1)
            throw MeshError(std::format("{} is in {} connections: a {} segment must be paired exactly once{}",
                                        seg_name(mesh, static_cast<int>(s)), used[s], segment_kind_name(k),
                                        used[s] == 0 ? " (unpaired)" : ""));
        if (!(k == SegmentKind::Periodic || k == SegmentKind::Interface) && used[s] != 0)
            throw MeshError(std::format("{} is in a connection: only periodic and interface segments are connected",
                                        seg_name(mesh, static_cast<int>(s))));
    }

    // tags: on their block, on a blade segment of their (blade, surface); one leading and one trailing each
    std::map<std::tuple<int, int, int>, int> count;  // (blade, surface, edge) -> tags
    for (std::size_t t = 0; t < mesh.tags.size(); ++t) {
        const EdgeTag& e = mesh.tags[t];
        if (e.block < 0 || e.block >= nb)
            throw MeshError(std::format("tag {} names block {}, the mesh has {}", t, e.block, nb));
        const MeshBlock& B = mesh.blocks[static_cast<std::size_t>(e.block)];
        if (!(0 <= e.i && e.i < B.ni && 0 <= e.j && e.j < B.nj))
            throw MeshError(std::format("tag {}: node ({}, {}) is not a node of block {} ({} x {})", t, e.i, e.j, e.block,
                                        B.ni, B.nj));
        bool on = false;
        for (const FaceSegment& g : mesh.segments) {
            if (g.block != e.block || g.kind != SegmentKind::Blade || g.blade != e.blade || g.surface != e.surface) continue;
            for (int q = g.begin; q <= g.end && !on; ++q) {
                const auto p = face_node(B, g.face, q);
                on = p[0] == e.i && p[1] == e.j;
            }
            if (on) break;
        }
        if (!on)
            throw MeshError(std::format("tag {} ({} edge of blade {} {} surface) at node ({}, {}) of block {} lies on no "
                                        "blade segment of that blade and surface",
                                        t, edge_kind_name(e.edge), e.blade, mesh_surface_name(e.surface), e.i, e.j, e.block));
        ++count[{e.blade, static_cast<int>(e.surface), static_cast<int>(e.edge)}];
    }
    for (const FaceSegment& g : mesh.segments) {
        if (g.kind != SegmentKind::Blade) continue;
        for (EdgeKind e : {EdgeKind::Leading, EdgeKind::Trailing}) {
            const int n = count[{g.blade, static_cast<int>(g.surface), static_cast<int>(e)}];
            if (n != 1)
                throw MeshError(std::format("the {} surface of blade {} has {} {}-edge tags: exactly one is needed (the "
                                            "Kutta condition and the surface distributions read them)",
                                            mesh_surface_name(g.surface), g.blade, n, edge_kind_name(e)));
        }
    }
}

}  // namespace BladeGen::Turbolenza
