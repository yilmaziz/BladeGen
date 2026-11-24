// The general curvilinear multi-block stream-function operator (see GeneralOperator.hpp).  Not a port
// of v1: the flux form of the user's derivation on any structured mesh.
#include "turbolenza/GeneralOperator.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <string>
#include <utility>

namespace BladeGen::Turbolenza {

namespace {

using T4 = std::array<int, 4>;

std::array<int, 2> times(const T4& A, std::array<int, 2> s) noexcept
{
    return {A[0] * s[0] + A[1] * s[1], A[2] * s[0] + A[3] * s[1]};
}

std::size_t z(int k) noexcept { return static_cast<std::size_t>(k); }

}  // namespace

std::array<double, 3> contravariant_tensor(double x_z, double x_e, double y_z, double y_e, double kappa)
{
    const double J = x_z * y_e - x_e * y_z;
    if (!(std::isfinite(J) && J != 0.0))
        throw GeneralOperatorError(std::format("the Jacobian x_zeta y_eta - x_eta y_zeta is {}: no tensor", J));
    return {kappa * (x_e * x_e + y_e * y_e) / J, -kappa * (x_z * x_e + y_z * y_e) / J, kappa * (x_z * x_z + y_z * y_z) / J};
}

double GeneralOperator::field(int g, std::vector<double> MeshBlock::*f) const
{
    const auto l = location(g);
    const MeshBlock& B = mesh_.blocks[z(l[0])];
    return (B.*f)[B.at(l[1], l[2])];
}

int GeneralOperator::node(int block, int i, int j) const
{
    const MeshBlock& B = mesh_.blocks.at(z(block));
    if (!(0 <= i && i < B.ni && 0 <= j && j < B.nj))
        throw GeneralOperatorError(std::format("node ({}, {}) is not in block {} ({} x {})", i, j, block, B.ni, B.nj));
    return offset_[z(block)] + i + B.ni * j;
}

std::array<int, 3> GeneralOperator::location(int g) const
{
    if (g < 0 || g >= n_nodes_) throw GeneralOperatorError(std::format("node {} is not in the mesh ({} nodes)", g, n_nodes_));
    const auto it = std::upper_bound(offset_.begin(), offset_.end(), g);
    const int b = static_cast<int>(it - offset_.begin()) - 1;
    const int k = g - offset_[z(b)];
    const int ni = mesh_.blocks[z(b)].ni;
    return {b, k % ni, k / ni};
}

std::optional<NodeRef> GeneralOperator::step(int g, int axis, int sign) const
{
    const auto [b, i, j] = location(g);
    const MeshBlock& B = mesh_.blocks[z(b)];
    const int ii = i + (axis == 0 ? sign : 0), jj = j + (axis == 1 ? sign : 0);
    if (0 <= ii && ii < B.ni && 0 <= jj && jj < B.nj) return NodeRef{offset_[z(b)] + ii + B.ni * jj, {1, 0, 0, 1}, 0};
    const BlockFace face = axis == 0 ? (sign < 0 ? BlockFace::IMin : BlockFace::IMax) : (sign < 0 ? BlockFace::JMin : BlockFace::JMax);
    const int t = axis == 0 ? j : i;
    for (int s : face_segs_[z(b)][static_cast<std::size_t>(face)][z(t)]) {
        const int c = conn_of_[z(s)];
        if (c < 0) continue;
        const Connection& cn = mesh_.connections[z(c)];
        const bool first = cn.first == s;
        const FaceSegment& here = mesh_.segments[z(s)];
        const FaceSegment& there = mesh_.segments[z(first ? cn.second : cn.first)];
        const int k = t - here.begin;
        const int pos = cn.reversed ? there.end - k : there.begin + k;
        const MeshBlock& P = mesh_.blocks[z(there.block)];
        const auto p = face_node(P, there.face, pos);
        const T4 T = face_transform(here.face, there.face, cn.reversed);
        const auto d = times(T, axis == 0 ? std::array<int, 2>{sign, 0} : std::array<int, 2>{0, sign});
        const int pi = p[0] + d[0], pj = p[1] + d[1];
        if (!(0 <= pi && pi < P.ni && 0 <= pj && pj < P.nj))
            throw GeneralOperatorError(std::format("the step out of block {} at ({}, {}) lands outside block {}", b, i, j,
                                                   there.block));
        return NodeRef{offset_[z(there.block)] + pi + P.ni * pj, T, first ? cn.image : -cn.image};
    }
    return std::nullopt;
}

NodeEquation GeneralOperator::equation(int g) const { return equation_[z(rep_[z(g)])]; }
int GeneralOperator::row_node(int g) const { return row_node_[z(rep_[z(g)])]; }

const std::vector<FaceGeometry>& GeneralOperator::faces(int g) const { return faces_[z(g)]; }

std::vector<int> GeneralOperator::segment_nodes(int s) const
{
    const FaceSegment& g = mesh_.segments.at(z(s));
    const MeshBlock& B = mesh_.blocks[z(g.block)];
    std::vector<int> out;
    for (int t = g.begin; t <= g.end; ++t) {
        const auto p = face_node(B, g.face, t);
        out.push_back(offset_[z(g.block)] + p[0] + B.ni * p[1]);
    }
    return out;
}

std::vector<int> GeneralOperator::segment_nodes(SegmentKind kind) const
{
    std::vector<int> out;
    std::vector<char> seen(z(n_nodes_), 0);
    for (std::size_t s = 0; s < mesh_.segments.size(); ++s) {
        if (mesh_.segments[s].kind != kind) continue;
        for (int g : segment_nodes(static_cast<int>(s)))
            if (!seen[z(g)]) {
                seen[z(g)] = 1;
                out.push_back(g);
            }
    }
    return out;
}

double GeneralOperator::apply(const Stencil& s, const std::vector<double>& values, double jump) const
{
    double sum = 0.0;
    for (const StencilEntry& e : s) sum += e.w * (values[z(e.node)] - static_cast<double>(e.images) * jump);
    return sum;
}

GeneralOperator::GeneralOperator(GeneralMesh mesh, const MeshCheckOptions& options) : mesh_(std::move(mesh))
{
    validate_mesh(mesh_, options);
    const int nb = static_cast<int>(mesh_.blocks.size());
    offset_.resize(z(nb));
    for (int b = 0; b < nb; ++b) {
        offset_[z(b)] = n_nodes_;
        n_nodes_ += mesh_.blocks[z(b)].nodes();
    }
    const std::size_t N = z(n_nodes_);
    face_segs_.resize(z(nb));
    for (int b = 0; b < nb; ++b)
        for (int f = 0; f < 4; ++f)
            face_segs_[z(b)][z(f)].resize(z(face_length(mesh_.blocks[z(b)], static_cast<BlockFace>(f))));
    conn_of_.assign(mesh_.segments.size(), -1);
    for (std::size_t c = 0; c < mesh_.connections.size(); ++c) {
        conn_of_[z(mesh_.connections[c].first)] = static_cast<int>(c);
        conn_of_[z(mesh_.connections[c].second)] = static_cast<int>(c);
    }
    wall_.assign(N, 0);
    wall_face_.assign(N, 0);
    inlet_.assign(N, 0);
    outlet_.assign(N, 0);
    for (std::size_t s = 0; s < mesh_.segments.size(); ++s) {
        const FaceSegment& g = mesh_.segments[s];
        for (int t = g.begin; t <= g.end; ++t) face_segs_[z(g.block)][static_cast<std::size_t>(g.face)][z(t)].push_back(static_cast<int>(s));
        for (int n : segment_nodes(static_cast<int>(s))) {
            if (g.kind == SegmentKind::Blade || g.kind == SegmentKind::Dirichlet) {
                wall_[z(n)] = 1;
                wall_face_[z(n)] = static_cast<char>(wall_face_[z(n)] | (1 << static_cast<int>(g.face)));
            }
            if (g.kind == SegmentKind::Inlet) inlet_[z(n)] = 1;
            if (g.kind == SegmentKind::Outlet) outlet_[z(n)] = 1;
        }
    }

    // node classes: union-find with offsets, psi_g = psi_parent + off_g
    std::vector<int> parent(N), off(N, 0);
    for (std::size_t g = 0; g < N; ++g) parent[g] = static_cast<int>(g);
    auto find = [&](int g) {
        int root = g, sum = 0;
        while (parent[z(root)] != root) {
            sum += off[z(root)];
            root = parent[z(root)];
        }
        // compress
        int cur = g, acc = sum;
        while (parent[z(cur)] != cur) {
            const int next = parent[z(cur)], o = off[z(cur)];
            parent[z(cur)] = root;
            off[z(cur)] = acc;
            acc -= o;
            cur = next;
        }
        return std::pair<int, int>{root, sum};
    };
    for (std::size_t c = 0; c < mesh_.connections.size(); ++c) {
        const Connection& cn = mesh_.connections[c];
        const std::vector<int> a = segment_nodes(cn.first), bnodes = segment_nodes(cn.second);
        for (std::size_t k = 0; k < a.size(); ++k) {
            const int na = a[k], nbn = bnodes[cn.reversed ? bnodes.size() - 1 - k : k];
            const auto [ra, oa] = find(na);
            const auto [rb, ob] = find(nbn);
            if (ra == rb) {
                if (ob != oa + cn.image)
                    throw GeneralOperatorError(std::format("connection {} closes a loop of periodic images with {} pitch(es) "
                                                           "instead of {}: inconsistent periodicity",
                                                           c, ob - oa, cn.image));
                continue;
            }
            parent[z(rb)] = ra;
            off[z(rb)] = oa + cn.image - ob;
        }
    }
    rep_.assign(N, -1);
    images_.assign(N, 0);
    std::vector<int> root_min(N, std::numeric_limits<int>::max()), root_off(N, 0);
    std::vector<std::pair<int, int>> rooted(N);
    for (std::size_t g = 0; g < N; ++g) {
        rooted[g] = find(static_cast<int>(g));
        if (static_cast<int>(g) < root_min[z(rooted[g].first)]) {
            root_min[z(rooted[g].first)] = static_cast<int>(g);
            root_off[z(rooted[g].first)] = rooted[g].second;
        }
    }
    for (std::size_t g = 0; g < N; ++g) {
        rep_[g] = root_min[z(rooted[g].first)];
        images_[g] = rooted[g].second - root_off[z(rooted[g].first)];
    }

    // fixed values and the relation of each class
    equation_.assign(N, NodeEquation::Pde);
    row_node_.assign(N, -1);
    fixed_.assign(N, std::numeric_limits<double>::quiet_NaN());
    dirichlet_.assign(N, 0);
    std::vector<char> fixed_set(N, 0);
    for (std::size_t s = 0; s < mesh_.segments.size(); ++s) {
        const FaceSegment& g = mesh_.segments[s];
        if (g.kind != SegmentKind::Blade && g.kind != SegmentKind::Dirichlet) continue;
        for (int n : segment_nodes(static_cast<int>(s))) {
            const std::size_t r = z(rep_[z(n)]);
            if (g.kind == SegmentKind::Blade) {
                const double v = static_cast<double>(g.blade - images_[z(n)]);
                if (fixed_set[r] && !dirichlet_[r] && fixed_[r] != v)
                    throw GeneralOperatorError(std::format("node {} lies on blade {} but its class is fixed to psi = {} "
                                                           "already: two walls of different psi share a node",
                                                           n, g.blade, fixed_[r] + images_[z(n)]));
                if (!fixed_set[r] || dirichlet_[r]) {
                    fixed_[r] = v;
                    dirichlet_[r] = 0;
                    row_node_[r] = n;
                }
            } else if (!fixed_set[r]) {
                dirichlet_[r] = 1;
                row_node_[r] = n;
            }
            fixed_set[r] = 1;
            equation_[r] = NodeEquation::Fixed;
        }
    }
    for (std::size_t g = 0; g < N; ++g) {
        const std::size_t r = z(rep_[g]);
        if (equation_[r] == NodeEquation::Fixed) continue;
        if (inlet_[g] || outlet_[g]) {
            if (equation_[r] == NodeEquation::FlowAngle && (inlet_[z(row_node_[r])] != inlet_[g]))
                throw GeneralOperatorError(std::format("node {} joins an inlet and an outlet node in one class", g));
            if (equation_[r] != NodeEquation::FlowAngle) {
                equation_[r] = NodeEquation::FlowAngle;
                row_node_[r] = static_cast<int>(g);
            }
        }
    }
    unknown_.assign(N, -1);
    for (std::size_t g = 0; g < N; ++g) {
        if (rep_[g] != static_cast<int>(g)) continue;
        if (equation_[g] == NodeEquation::Pde) row_node_[g] = static_cast<int>(g);
        if (equation_[g] != NodeEquation::Fixed) unknown_[g] = n_unknowns_++;
    }

    // nodal differences
    auto describe = [&](int g) {
        const auto l = location(g);
        return std::format("node ({}, {}) of block {}", l[1], l[2], l[0]);
    };
    // Along a boundary face F of X (tangential to the difference), the neighbour Q must be a wall node of
    // F's continuation exactly when X is a wall node of F (Aungier's runs: no difference across a
    // periodic / wall transition).  F's continuation in Q's block is the face whose outward normal is
    // the image of F's under the step's index map.
    auto on_face = [&](int g, BlockFace f) {
        const auto [b, i, j] = location(g);
        const MeshBlock& B = mesh_.blocks[z(b)];
        switch (f) {
        case BlockFace::IMin: return i == 0;
        case BlockFace::IMax: return i == B.ni - 1;
        case BlockFace::JMin: return j == 0;
        case BlockFace::JMax: return j == B.nj - 1;
        }
        return false;
    };
    auto face_of_normal = [](std::array<int, 2> n) {
        if (n[0] != 0) return n[0] < 0 ? BlockFace::IMin : BlockFace::IMax;
        return n[1] < 0 ? BlockFace::JMin : BlockFace::JMax;
    };
    auto wall_on = [&](int g, BlockFace f) { return (wall_face_[z(g)] >> static_cast<int>(f)) & 1; };
    auto compatible = [&](int X, const NodeRef& Q, int axis) {
        for (int side : {-1, 1}) {
            const std::array<int, 2> n = axis == 0 ? std::array<int, 2>{0, side} : std::array<int, 2>{side, 0};
            const BlockFace F = face_of_normal(n);
            if (!on_face(X, F)) continue;
            const BlockFace FQ = face_of_normal(times(Q.T, n));
            const int wq = on_face(Q.node, FQ) ? wall_on(Q.node, FQ) : 0;
            if (wall_on(X, F) != wq) return false;
        }
        return true;
    };
    auto unit_axis = [](std::array<int, 2> d) { return d[0] != 0 ? std::pair<int, int>{0, d[0]} : std::pair<int, int>{1, d[1]}; };
    deriv_.assign(2 * N, {});
    for (int g = 0; g < n_nodes_; ++g)
        for (int a = 0; a < 2; ++a) {
            const auto plus = step(g, a, 1), minus = step(g, a, -1);
            const bool cp = plus && compatible(g, *plus, a), cm = minus && compatible(g, *minus, a);
            Stencil& S = deriv_[z(2 * g + a)];
            if (cp && cm) {
                S = {{plus->node, 0.5, plus->images}, {minus->node, -0.5, minus->images}};
                continue;
            }
            bool done = false;
            for (int sign : {1, -1}) {
                const std::optional<NodeRef>& q1 = sign > 0 ? plus : minus;
                if (done || !(sign > 0 ? cp : cm)) continue;
                const auto [a1, s1] = unit_axis(times(q1->T, a == 0 ? std::array<int, 2>{sign, 0} : std::array<int, 2>{0, sign}));
                const auto q2 = step(q1->node, a1, s1);
                if (!q2 || !compatible(q1->node, *q2, a1)) continue;
                const double w = static_cast<double>(sign);
                S = {{g, -1.5 * w, 0}, {q1->node, 2.0 * w, q1->images}, {q2->node, -0.5 * w, q1->images + q2->images}};
                done = true;
            }
            if (!done)
                throw GeneralOperatorError(std::format("{} has no second-order difference along {} (no central and no "
                                                       "one-sided stencil of three nodes in one run)",
                                                       describe(g), a == 0 ? "i" : "j"));
        }
    std::vector<double> X(N), Y(N);
    for (int g = 0; g < n_nodes_; ++g) {
        X[z(g)] = x(g);
        Y[z(g)] = y(g);
    }
    metric_.resize(N);
    for (int g = 0; g < n_nodes_; ++g) {
        NodalMetric& mt = metric_[z(g)];
        mt.x_z = apply(derivative(g, 0), X, 0.0);
        mt.x_e = apply(derivative(g, 1), X, 0.0);
        mt.y_z = apply(derivative(g, 0), Y, mesh_.pitch);
        mt.y_e = apply(derivative(g, 1), Y, mesh_.pitch);
        mt.J = mt.x_z * mt.y_e - mt.x_e * mt.y_z;
        if (!(mt.J > 0.0))
            throw GeneralOperatorError(std::format("{}: the nodal Jacobian is {:.3e} <= 0", describe(g), mt.J));
    }

    // The control volume of each PDE class (ours, the cell-vertex dual of the user's scheme): around the
    // class's point, the quarter of every cell incident to a member (in that member's block, shifted into
    // the representative's frame by the member's images).  Its boundary, counterclockwise, runs from cell
    // centroid to edge midpoint to the next centroid.  The two half-edges on the edge to one neighbour make
    // one face from centroid C1 to centroid C2: t = C2 - C1, d = Q - P, and the tangential difference
    // psi(C2) - psi(C1) with psi(C) the mean of the cell's four corners.  For a node with four neighbours
    // in its own block this is the nine-point scheme of the file comment exactly.  Where the half-edges
    // of two members end at different nodes (an H-grid transition station: the two images of the last
    // inlet / first exit station meet the suction wall node (f, 0) and the pressure wall node (f, J) one
    // pitch down), the two half-edges and the gap between their midpoints are closed by one face from
    // C1 to C2 to the mean of the two nodes, Q* = (Q1 + Q2) / 2, psi(Q*) = (psi(Q1) + psi(Q2)) / 2: the
    // control volume is a closed polygon, so the discrete metric identities hold to round-off there too.
    members_.assign(N, {});
    for (int g = 0; g < n_nodes_; ++g) members_[z(rep_[z(g)])].push_back(g);
    faces_.assign(N, {});
    source_weight_.assign(N, 0.0);
    struct HalfEdge
    {
        int q = 0, q_images = 0;    // the neighbour and the images of its block relative to the representative
        int key_rep = 0, key_img = 0;
        bool ends_at_mid = false;   // centroid -> midpoint (else midpoint -> centroid)
        double cx = 0.0, cy = 0.0;  // the cell's centroid (representative's frame)
        double mx = 0.0, my = 0.0;  // the edge midpoint
        Stencil cpsi;               // psi at the centroid
        int member = 0, axis = 0, sign = 0;
    };
    for (int R = 0; R < n_nodes_; ++R) {
        if (rep_[z(R)] != R || equation_[z(R)] != NodeEquation::Pde) continue;
        std::vector<HalfEdge> half;
        for (int M : members_[z(R)]) {
            const auto [bM, iM, jM] = location(M);
            const MeshBlock& BM = mesh_.blocks[z(bM)];
            source_weight_[z(M)] = ((iM == 0 || iM == BM.ni - 1) ? 0.5 : 1.0) * ((jM == 0 || jM == BM.nj - 1) ? 0.5 : 1.0);
            const int m = images_[z(M)];
            const double shift = static_cast<double>(m) * mesh_.pitch;
            const int base = offset_[z(bM)], ni = BM.ni;
            auto at = [base, ni](int i, int j) { return base + i + ni * j; };
            for (int di : {-1, 1})
                for (int dj : {-1, 1}) {
                    if (!(0 <= iM + di && iM + di < BM.ni && 0 <= jM + dj && jM + dj < BM.nj)) continue;
                    const int A = at(iM + di, jM), Bn = at(iM, jM + dj), D = at(iM + di, jM + dj);
                    const double cx = 0.25 * (X[z(M)] + X[z(A)] + X[z(Bn)] + X[z(D)]);
                    const double cy = 0.25 * (Y[z(M)] + Y[z(A)] + Y[z(Bn)] + Y[z(D)]) - shift;
                    const Stencil cpsi = {{M, 0.25, m}, {A, 0.25, m}, {Bn, 0.25, m}, {D, 0.25, m}};
                    // counterclockwise round M in a right-handed cell: first the neighbour along i when
                    // di dj > 0, else the one along j; the first half-edge leaves its midpoint, the second
                    // ends at its midpoint
                    const bool i_first = di * dj > 0;
                    for (int k = 0; k < 2; ++k) {
                        const bool along_i = (k == 0) == i_first;
                        HalfEdge h;
                        h.q = along_i ? A : Bn;
                        h.q_images = m;
                        h.key_rep = rep_[z(h.q)];
                        h.key_img = images_[z(h.q)] - m;
                        h.ends_at_mid = k == 1;
                        h.cx = cx;
                        h.cy = cy;
                        h.mx = 0.5 * (X[z(M)] + X[z(h.q)]);
                        h.my = 0.5 * (Y[z(M)] + Y[z(h.q)]) - shift;
                        h.cpsi = cpsi;
                        h.member = M;
                        h.axis = along_i ? 0 : 1;
                        h.sign = along_i ? di : dj;
                        half.push_back(std::move(h));
                    }
                }
        }
        const double px = X[z(R)], py = Y[z(R)];
        auto make_face = [&](const HalfEdge& in, const HalfEdge& out, bool closure) {
            FaceGeometry f;
            f.closure = closure;
            if (!closure && in.member == R) {
                f.axis = in.axis;
                f.sign = in.sign;
            } else if (!closure && out.member == R) {
                f.axis = out.axis;
                f.sign = out.sign;
            }
            const double w = closure ? 0.5 : 1.0;
            const double yin = Y[z(in.q)] - static_cast<double>(in.q_images) * mesh_.pitch;
            const double yout = Y[z(out.q)] - static_cast<double>(out.q_images) * mesh_.pitch;
            const double qx = closure ? 0.5 * (X[z(in.q)] + X[z(out.q)]) : X[z(in.q)];
            const double qy = closure ? 0.5 * (yin + yout) : yin;
            f.neighbours = closure ? std::vector<int>{in.q, out.q} : std::vector<int>{in.q};
            f.x_n = qx - px;
            f.y_n = qy - py;
            f.x_t = out.cx - in.cx;
            f.y_t = out.cy - in.cy;
            f.J = f.x_n * f.y_t - f.y_n * f.x_t;
            if (!(f.J > 0.0))
                throw GeneralOperatorError(std::format("{}: a face of its control volume has the Jacobian {:.3e} <= 0",
                                                       describe(R), f.J));
            f.G_nn = (f.x_t * f.x_t + f.y_t * f.y_t) / f.J;
            f.G_nt = -(f.x_n * f.x_t + f.y_n * f.y_t) / f.J;
            f.dpsi = {{R, -1.0, 0}, {in.q, w, in.q_images}};
            if (closure) f.dpsi.push_back({out.q, w, out.q_images});
            f.kappa = {{R, 0.5}, {in.q, 0.5 * w}};
            if (closure) f.kappa.emplace_back(out.q, 0.5 * w);
            for (const StencilEntry& e : out.cpsi) f.psi_t.push_back(e);
            for (const StencilEntry& e : in.cpsi) f.psi_t.push_back({e.node, -e.w, e.images});
            faces_[z(R)].push_back(std::move(f));
        };
        // pair the two half-edges of each neighbour; what remains is closed across the gaps
        std::vector<char> used(half.size(), 0);
        for (std::size_t a = 0; a < half.size(); ++a) {
            if (used[a] || !half[a].ends_at_mid) continue;
            int match = -1, count = 0;
            for (std::size_t c = 0; c < half.size(); ++c)
                if (!used[c] && !half[c].ends_at_mid && half[c].key_rep == half[a].key_rep &&
                    half[c].key_img == half[a].key_img) {
                    match = static_cast<int>(c);
                    ++count;
                }
            if (count > 1)
                throw GeneralOperatorError(std::format("{}: the control volume meets one neighbour on more than two "
                                                       "half-edges (a singular vertex is not supported)",
                                                       describe(R)));
            if (match < 0) continue;
            used[a] = used[z(match)] = 1;
            make_face(half[a], half[z(match)], false);
        }
        std::vector<std::size_t> ends, starts;
        for (std::size_t a = 0; a < half.size(); ++a)
            if (!used[a]) (half[a].ends_at_mid ? ends : starts).push_back(a);
        if (ends.size() != starts.size())
            throw GeneralOperatorError(std::format("{}: its control volume does not close ({} open half-edges in, {} out)",
                                                   describe(R), ends.size(), starts.size()));
        for (std::size_t a : ends) {  // the next half-edge counterclockwise starts at the nearest midpoint
            std::size_t best = starts.size();
            double dist = 0.0;
            for (std::size_t k = 0; k < starts.size(); ++k) {
                if (used[starts[k]]) continue;
                const double d = std::hypot(half[starts[k]].mx - half[a].mx, half[starts[k]].my - half[a].my);
                if (best == starts.size() || d < dist) {
                    best = k;
                    dist = d;
                }
            }
            used[a] = used[starts[best]] = 1;
            make_face(half[a], half[starts[best]], true);
        }
    }
}

namespace {

// psi of node n of a fixed class (rep r): the wall's value or the caller's Dirichlet value.
double fixed_value(int r, int row, const std::vector<double>& fixed, const std::vector<char>& dirichlet,
                   const std::vector<int>& images, const OperatorFields& fields)
{
    if (!dirichlet[static_cast<std::size_t>(r)]) return fixed[static_cast<std::size_t>(r)];
    if (!fields.dirichlet) throw GeneralOperatorError("the mesh has Dirichlet segments but no Dirichlet values were given");
    const double v = fields.dirichlet(row);
    if (!std::isfinite(v)) throw GeneralOperatorError(std::format("the Dirichlet value at node {} is not finite", row));
    return v - static_cast<double>(images[static_cast<std::size_t>(row)]);
}

void check_fields(const OperatorFields& f, int n, bool need_tan)
{
    auto check = [&](const std::vector<double>& v, const char* name, bool required) {
        if (v.empty() && !required) return;
        if (static_cast<int>(v.size()) != n)
            throw GeneralOperatorError(std::format("{} has {} values, the mesh {} nodes", name, v.size(), n));
        for (double x : v)
            if (!std::isfinite(x)) throw GeneralOperatorError(std::format("{} is not finite", name));
    };
    check(f.coefficient, "the coefficient mdot' / (rho b)", true);
    for (double c : f.coefficient)
        if (!(c > 0.0)) throw GeneralOperatorError("the coefficient mdot' / (rho b) must be positive at every node");
    check(f.source, "the source", false);
    check(f.flow_angle_tan, "the flow-angle tangent", need_tan);
}

}  // namespace

// The relation of the class of representative r, as sum w psi(node) - images, handed to emit(node, w,
// images); returns the right-hand side before the constants (the source of a PDE row, else 0).
template <class Emit>
static double relation(const GeneralOperator& op, int r, const OperatorFields& fields, Emit&& emit)
{
    switch (op.equation(r)) {
    case NodeEquation::Pde: {
        for (const FaceGeometry& f : op.faces(r)) {
            double kf = 0.0;
            for (const auto& [node, w] : f.kappa) kf += w * fields.coefficient[static_cast<std::size_t>(node)];
            for (const StencilEntry& e : f.dpsi) emit(e.node, kf * f.G_nn * e.w, e.images);
            for (const StencilEntry& e : f.psi_t) emit(e.node, kf * f.G_nt * e.w, e.images);
        }
        double source = 0.0;
        if (!fields.source.empty())
            for (int M : op.members(r)) source += op.source_weight(M) * fields.source[static_cast<std::size_t>(M)];
        return source;
    }
    case NodeEquation::FlowAngle: {
        const int R = op.row_node(r);
        const NodalMetric& mt = op.metric(R);
        const double T = fields.flow_angle_tan[static_cast<std::size_t>(R)];
        const double az = mt.y_e - T * mt.x_e, ae = T * mt.x_z - mt.y_z;
        for (const StencilEntry& e : op.derivative(R, 0)) emit(e.node, az * e.w, e.images);
        for (const StencilEntry& e : op.derivative(R, 1)) emit(e.node, ae * e.w, e.images);
        return 0.0;
    }
    case NodeEquation::Fixed: break;
    }
    return 0.0;
}

GeneralSystem GeneralOperator::assemble(const OperatorFields& fields) const
{
    check_fields(fields, n_nodes_, !segment_nodes(SegmentKind::Inlet).empty() || !segment_nodes(SegmentKind::Outlet).empty());
    GeneralSystem sys;
    sys.rhs.assign(z(n_unknowns_), 0.0);
    std::vector<Eigen::Triplet<double, int>> triplets;
    std::vector<std::pair<int, double>> entries;
    for (int r = 0; r < n_nodes_; ++r) {
        if (rep_[z(r)] != r || unknown_[z(r)] < 0) continue;
        const int row = unknown_[z(r)];
        double constant = 0.0;
        entries.clear();
        const double source = relation(*this, r, fields, [&](int node, double w, int img) {
            const int rn = rep_[z(node)];
            const double shift = static_cast<double>(images_[z(node)] - img);
            const int u = unknown_[z(rn)];
            if (u >= 0) {
                entries.emplace_back(u, w);
                constant += w * shift;
            } else {
                constant += w * (fixed_value(rn, row_node_[z(rn)], fixed_, dirichlet_, images_, fields) + shift);
            }
        });
        // each row divided by its own unknown's coefficient (the analogue of Aungier's division of Eq.
        // 5-44 by -H: centre 1; it changes no solution, only the conditioning of the factorization)
        double diagonal = 0.0;
        for (const auto& [u, w] : entries)
            if (u == row) diagonal += w;
        if (!(std::isfinite(diagonal) && diagonal != 0.0))
            throw GeneralOperatorError(std::format("the relation of node {} has the coefficient {} on its own unknown", r,
                                                   diagonal));
        for (const auto& [u, w] : entries) triplets.emplace_back(row, u, w / diagonal);
        sys.rhs[z(row)] = (source - constant) / diagonal;
    }
    sys.matrix.resize(n_unknowns_, n_unknowns_);
    sys.matrix.setFromTriplets(triplets.begin(), triplets.end());
    sys.matrix.prune([](int, int, double v) { return v != 0.0; });
    sys.matrix.makeCompressed();
    return sys;
}

std::vector<double> GeneralOperator::expand(const std::vector<double>& unknowns, const OperatorFields& fields) const
{
    if (static_cast<int>(unknowns.size()) != n_unknowns_)
        throw GeneralOperatorError(std::format("{} unknowns given, the system has {}", unknowns.size(), n_unknowns_));
    std::vector<double> psi(z(n_nodes_));
    for (int g = 0; g < n_nodes_; ++g) {
        const int r = rep_[z(g)];
        const int u = unknown_[z(r)];
        const double base = u >= 0 ? unknowns[z(u)] : fixed_value(r, row_node_[z(r)], fixed_, dirichlet_, images_, fields);
        psi[z(g)] = base + static_cast<double>(images_[z(g)]);
    }
    return psi;
}

std::vector<double> GeneralOperator::residual(const std::vector<double>& psi, const OperatorFields& fields) const
{
    check_fields(fields, n_nodes_, !segment_nodes(SegmentKind::Inlet).empty() || !segment_nodes(SegmentKind::Outlet).empty());
    if (static_cast<int>(psi.size()) != n_nodes_) throw GeneralOperatorError("psi must have one value per node");
    std::vector<double> per_rep(z(n_nodes_), 0.0), out(z(n_nodes_), 0.0);
    for (int r = 0; r < n_nodes_; ++r) {
        if (rep_[z(r)] != r) continue;
        if (equation_[z(r)] == NodeEquation::Fixed) {
            const int row = row_node_[z(r)];
            per_rep[z(r)] = psi[z(row)] - (fixed_value(r, row, fixed_, dirichlet_, images_, fields) +
                                           static_cast<double>(images_[z(row)]));
            continue;
        }
        // values in the frame of the relation's node: psi(node) - images
        double sum = 0.0;
        const double source = relation(*this, r, fields, [&](int node, double w, int img) {
            sum += w * (psi[z(node)] - static_cast<double>(img));
        });
        // the PDE row's own node (r) and the flow-angle node R are in their own frames already, except
        // that psi(R) carries images(R) relative to r: a constant shift that the differences cancel.
        per_rep[z(r)] = sum - source;
    }
    for (int g = 0; g < n_nodes_; ++g) out[z(g)] = per_rep[z(rep_[z(g)])];
    return out;
}

StreamGradient GeneralOperator::gradient(const std::vector<double>& psi) const
{
    if (static_cast<int>(psi.size()) != n_nodes_) throw GeneralOperatorError("psi must have one value per node");
    StreamGradient G;
    G.psi_x.assign(z(n_nodes_), 0.0);
    G.psi_y.assign(z(n_nodes_), 0.0);
    for (int g = 0; g < n_nodes_; ++g) {
        if (rep_[z(g)] != g) continue;
        const NodalMetric& mt = metric_[z(g)];
        const double pz = apply(derivative(g, 0), psi, 1.0), pe = apply(derivative(g, 1), psi, 1.0);
        G.psi_x[z(g)] = (mt.y_e * pz - mt.y_z * pe) / mt.J;
        G.psi_y[z(g)] = (mt.x_z * pe - mt.x_e * pz) / mt.J;
    }
    for (int g = 0; g < n_nodes_; ++g) {
        const int r = rep_[z(g)];
        G.psi_x[z(g)] = G.psi_x[z(r)];
        G.psi_y[z(g)] = G.psi_y[z(r)];
    }
    return G;
}

std::vector<double> GeneralOperator::rotation_source(double omega) const
{
    std::vector<double> r2(z(n_nodes_)), out(z(n_nodes_), 0.0);
    for (int g = 0; g < n_nodes_; ++g) r2[z(g)] = r(g) * r(g);
    if (omega == 0.0) return out;
    for (int g = 0; g < n_nodes_; ++g) {
        const NodalMetric& mt = metric_[z(g)];
        out[z(g)] = omega * (mt.y_e * apply(derivative(g, 0), r2, 0.0) - mt.y_z * apply(derivative(g, 1), r2, 0.0));
    }
    return out;
}

std::vector<std::array<double, 2>> GeneralOperator::gcl_residual() const
{
    std::vector<std::array<double, 2>> out(z(n_nodes_), {0.0, 0.0});
    for (int g = 0; g < n_nodes_; ++g)
        for (const FaceGeometry& f : faces_[z(g)]) {  // the outward normals (t_y, -t_x) of the boundary
            out[z(g)][0] += f.y_t;
            out[z(g)][1] -= f.x_t;
        }
    for (int g = 0; g < n_nodes_; ++g) out[z(g)] = out[z(rep_[z(g)])];
    return out;
}

}  // namespace BladeGen::Turbolenza
