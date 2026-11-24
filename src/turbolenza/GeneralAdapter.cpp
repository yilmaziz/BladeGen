// TURBOlenza's Aungier grid as a GeneralMesh (see GeneralAdapter.hpp).
#include "turbolenza/GeneralAdapter.hpp"

#include <cmath>
#include <format>

namespace BladeGen::Turbolenza {

GeneralMesh mesh_from_grid(const Grid& grid, double c_ref)
{
    if (!(std::isfinite(c_ref) && c_ref > 0.0)) throw MeshError(std::format("c_ref = {}: a positive chord [m] is needed", c_ref));
    if (grid.touching()) throw MeshError("the touching edge-node layout is not adapted (zero-thickness wall stations)");
    const int nm = grid.n_m(), ne = grid.n_eta(), f = grid.first_blade_i, l = grid.last_blade_i;
    if (!(f >= 2 && l <= nm - 3 && f < l))
        throw MeshError(std::format("the grid's blade run [{}, {}] of {} stations leaves fewer than two periodic stations "
                                    "on a side",
                                    f, l, nm));
    const WallIdentity& id = grid.wall_identity;
    auto one_surface = [&](const std::vector<Surface>& s, const char* wall) {
        if (s.empty()) throw MeshError(std::format("the {} wall has no identity", wall));
        for (Surface x : s)
            if (x != s.front()) throw MeshError(std::format("the {} wall's nodes lie on both surfaces", wall));
        return s.front();
    };
    const Surface lower = one_surface(id.lower.surface, "lower"), upper = one_surface(id.upper.surface, "upper");
    GeneralMesh mesh;
    mesh.description = std::format("TURBOlenza Aungier grid {} x {} (blade stations {} .. {}), c_ref = {} m", nm, ne, f, l, c_ref);
    mesh.pitch = grid.pitch();
    MeshBlock B;
    B.name = "aungier";
    B.ni = nm;
    B.nj = ne;
    const std::size_t n = static_cast<std::size_t>(nm) * static_cast<std::size_t>(ne);
    B.x.resize(n);
    B.y.resize(n);
    B.m.resize(n);
    B.r.assign(n, c_ref);
    B.b.assign(n, 1.0);
    for (int j = 0; j < ne; ++j)
        for (int i = 0; i < nm; ++i) {
            const std::size_t k = B.at(i, j);
            B.x[k] = grid.node_m(i, j);
            B.y[k] = grid.node_y(i, j);
            B.m[k] = c_ref * grid.node_m(i, j);
        }
    mesh.blocks.push_back(std::move(B));
    auto seg = [&](BlockFace face, int a, int b, SegmentKind kind, int blade = 0, Surface s = Surface::Suction) {
        FaceSegment g;
        g.block = 0;
        g.face = face;
        g.begin = a;
        g.end = b;
        g.kind = kind;
        g.blade = blade;
        g.surface = s;
        mesh.segments.push_back(g);
    };
    seg(BlockFace::JMin, 0, f - 1, SegmentKind::Periodic);
    seg(BlockFace::JMin, f, l, SegmentKind::Blade, 0, lower);
    seg(BlockFace::JMin, l + 1, nm - 1, SegmentKind::Periodic);
    seg(BlockFace::JMax, 0, f - 1, SegmentKind::Periodic);
    seg(BlockFace::JMax, f, l, SegmentKind::Blade, 1, upper);
    seg(BlockFace::JMax, l + 1, nm - 1, SegmentKind::Periodic);
    seg(BlockFace::IMin, 0, ne - 1, SegmentKind::Inlet);
    seg(BlockFace::IMax, 0, ne - 1, SegmentKind::Outlet);
    mesh.connections.push_back({0, 3, false, 1});
    mesh.connections.push_back({2, 5, false, 1});
    for (EdgeKind e : {EdgeKind::Leading, EdgeKind::Trailing}) {
        const int i = e == EdgeKind::Leading ? f : l;
        mesh.tags.push_back({e, 0, lower, 0, i, 0});
        mesh.tags.push_back({e, 1, upper, 0, i, ne - 1});
    }
    validate_mesh(mesh);
    return mesh;
}

}  // namespace BladeGen::Turbolenza
