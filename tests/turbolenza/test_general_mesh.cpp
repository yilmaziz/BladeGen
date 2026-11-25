// The general mesh (GeneralMesh.hpp), its file format "turbolenza.mesh" version 1 (MeshFormat.hpp) and
// the adapter of TURBOlenza's Aungier grid (GeneralAdapter.hpp): the index map of a connection, the
// refusals of validate_mesh and of the reader, the bit-exact round trip, the adapted grid.
#include "CascadeFixture.hpp"
#include "GeneralFixture.hpp"

#include "turbolenza/GeneralAdapter.hpp"
#include "turbolenza/GeneralMesh.hpp"
#include "turbolenza/Grid.hpp"
#include "turbolenza/MeshFormat.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace BladeGen::Turbolenza;

namespace {

/// A valid one-passage mesh with blades: the Aungier grid of the synthetic blade (no reference data).
GeneralMesh synthetic_mesh()
{
    const Cascade c = make_cascade(fixture::synthetic_section(), 1.25, 30.0);
    const StandardGrid g = build_grid(c, "very_coarse", 45.0 * fixture::kPi / 180.0, std::nullopt);
    return mesh_from_grid(g.grid, 0.5);
}

void same(const GeneralMesh& a, const GeneralMesh& b)
{
    REQUIRE(a.blocks.size() == b.blocks.size());
    CHECK(a.pitch == b.pitch);
    CHECK(a.description == b.description);
    for (std::size_t k = 0; k < a.blocks.size(); ++k) {
        const MeshBlock &A = a.blocks[k], &B = b.blocks[k];
        CHECK(A.name == B.name);
        CHECK(A.ni == B.ni);
        CHECK(A.nj == B.nj);
        CHECK(A.x == B.x);  // bit for bit (17 significant digits)
        CHECK(A.y == B.y);
        CHECK(A.m == B.m);
        CHECK(A.r == B.r);
        CHECK(A.b == B.b);
    }
    REQUIRE(a.segments.size() == b.segments.size());
    for (std::size_t k = 0; k < a.segments.size(); ++k) {
        const FaceSegment &A = a.segments[k], &B = b.segments[k];
        CHECK((A.block == B.block && A.face == B.face && A.begin == B.begin && A.end == B.end && A.kind == B.kind));
        if (A.kind == SegmentKind::Blade) CHECK((A.blade == B.blade && A.surface == B.surface));
    }
    REQUIRE(a.connections.size() == b.connections.size());
    for (std::size_t k = 0; k < a.connections.size(); ++k) {
        const Connection &A = a.connections[k], &B = b.connections[k];
        CHECK((A.first == B.first && A.second == B.second && A.reversed == B.reversed && A.image == B.image));
    }
    REQUIRE(a.tags.size() == b.tags.size());
    for (std::size_t k = 0; k < a.tags.size(); ++k) {
        const EdgeTag &A = a.tags[k], &B = b.tags[k];
        CHECK((A.edge == B.edge && A.blade == B.blade && A.surface == B.surface && A.block == B.block && A.i == B.i &&
               A.j == B.j));
    }
}

std::string replaced(std::string s, const std::string& from, const std::string& to)
{
    const auto at = s.find(from);
    REQUIRE(at != std::string::npos);
    return s.replace(at, from.size(), to);
}

}  // namespace

TEST_CASE("general mesh: the index map of a connection", "[turbolenza][general]")
{
    using F = BlockFace;
    // jmin -> jmax (an H block's periodic faces): the identity
    CHECK(face_transform(F::JMin, F::JMax, false) == std::array<int, 4>{1, 0, 0, 1});
    // imax -> imin (an O seam): the identity
    CHECK(face_transform(F::IMax, F::IMin, false) == std::array<int, 4>{1, 0, 0, 1});
    // imax -> imax reversed (a block turned by 180 degrees): -I
    CHECK(face_transform(F::IMax, F::IMax, true) == std::array<int, 4>{-1, 0, 0, -1});
    // imin -> jmax reversed (a block turned by 90 degrees): a rotation
    CHECK(face_transform(F::IMin, F::JMax, true) == std::array<int, 4>{0, -1, 1, 0});
    // jmin -> jmax reversed folds (determinant -1)
    const auto T = face_transform(F::JMin, F::JMax, true);
    CHECK(T[0] * T[3] - T[1] * T[2] == -1);
}

TEST_CASE("general mesh: validate_mesh refuses what breaks the rules", "[turbolenza][general]")
{
    const GeneralMesh good = synthetic_mesh();
    REQUIRE_NOTHROW(validate_mesh(good));
    {
        GeneralMesh m = good;  // an unpaired periodic segment
        m.connections.pop_back();
        CHECK_THROWS_WITH(validate_mesh(m), Catch::Matchers::ContainsSubstring("unpaired"));
    }
    {
        // interface node-count mismatch (nonmatching): block 1's U = 0.5 face split in two
        GeneralMesh m = gfix::three_block_periodic(6, 0.8, SegmentKind::Inlet, SegmentKind::Outlet);
        REQUIRE_NOTHROW(validate_mesh(m));
        m.segments[4].end -= 2;
        m.segments.push_back(gfix::segment(1, BlockFace::IMax, m.segments[4].end, m.segments[4].end + 2, SegmentKind::Interface));
        CHECK_THROWS_WITH(validate_mesh(m), Catch::Matchers::ContainsSubstring("nonmatching"));
    }
    {
        GeneralMesh m = good;  // a folded cell (negative Jacobian)
        MeshBlock& B = m.blocks[0];
        std::swap(B.y[B.at(5, 2)], B.y[B.at(5, 3)]);
        CHECK_THROWS_WITH(validate_mesh(m), Catch::Matchers::ContainsSubstring("Jacobian"));
    }
    {
        GeneralMesh m = good;  // a missing blade tag
        m.tags.pop_back();
        CHECK_THROWS_WITH(validate_mesh(m), Catch::Matchers::ContainsSubstring("edge tags"));
    }
    {
        GeneralMesh m = good;  // a tag off its blade segment
        m.tags[0].i -= 1;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
    {
        GeneralMesh m = good;  // periodic images that do not coincide
        m.pitch *= 1.001;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
    {
        GeneralMesh m = good;  // a connection that folds the orientation
        m.connections[0].reversed = true;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
    {
        GeneralMesh m = good;  // a face node on no segment
        m.segments[6].end -= 1;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
    {
        GeneralMesh m = good;  // two segments sharing a cell
        m.segments[0].end += 1;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
    {
        GeneralMesh m = good;  // a non-positive radius, a too small block
        m.blocks[0].r[3] = 0.0;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
        GeneralMesh t = gfix::dirichlet_box(gfix::block(gfix::family(1), 2, 5));
        CHECK_THROWS_AS(validate_mesh(t), MeshError);
    }
    {
        GeneralMesh m = good;  // a periodic connection without an image, an interface with one
        m.connections[0].image = 0;
        CHECK_THROWS_AS(validate_mesh(m), MeshError);
    }
}

TEST_CASE("general mesh: the file format round trip is bit for bit", "[turbolenza][general]")
{
    for (int k = 0; k < 2; ++k) {
        const GeneralMesh mesh =
            k == 0 ? synthetic_mesh() : gfix::three_block_periodic(6, 0.8, SegmentKind::Inlet, SegmentKind::Outlet);
        const std::string text = write_mesh_json(mesh);
        const GeneralMesh back = read_mesh_json(text);
        same(mesh, back);
        CHECK(write_mesh_json(back) == text);
    }
    // through a file
    const GeneralMesh mesh = synthetic_mesh();
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "turbolenza_general_roundtrip.mesh.json";
    write_mesh_file(mesh, path.string());
    same(mesh, read_mesh_file(path.string()));
    std::filesystem::remove(path);
    CHECK_THROWS_AS(read_mesh_file((std::filesystem::temp_directory_path() / "no_such_mesh.json").string()), MeshFormatError);
}

TEST_CASE("general mesh: the reader refuses an inconsistent or foreign file", "[turbolenza][general]")
{
    const std::string text = write_mesh_json(synthetic_mesh());
    CHECK_THROWS_AS(read_mesh_json("{"), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"turbolenza.mesh\"", "\"other.mesh\"")), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"version\": 1", "\"version\": 2")), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"theta\": \"rad\"", "\"theta\": \"deg\"")), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"pitch\":", "\"extra\": 1, \"pitch\":")), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"kind\": \"periodic\", \"segments\": [0, 3], \"reversed\": false, \"image\": 1},", "")),
                    MeshError);  // unpaired periodic segments
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"kind\": \"inlet\"", "\"kind\": \"dirichlet\"")), MeshFormatError);
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"surface\": \"suction\", \"block\"", "\"surface\": \"pressure\", \"block\"")),
                    MeshError);  // a tag of a surface with no segment
    CHECK_THROWS_AS(read_mesh_json(replaced(text, "\"reversed\": false, \"image\": 1", "\"reversed\": false")), MeshFormatError);
    // a key twice in one object; integers that are floats or out of the range of int
    CHECK_THROWS_WITH(read_mesh_json(replaced(text, "\"pitch\":", "\"pitch\": 1.0, \"pitch\":")),
                      Catch::Matchers::ContainsSubstring("twice"));
    CHECK_THROWS_WITH(read_mesh_json(replaced(text, "\"block\": 0, \"face\": \"imin\"", "\"block\": 0, \"block\": 0, \"face\": \"imin\"")),
                      Catch::Matchers::ContainsSubstring("twice"));
    CHECK_THROWS_WITH(read_mesh_json(replaced(text, "\"version\": 1", "\"version\": 1.0")), Catch::Matchers::ContainsSubstring("integer"));
    CHECK_THROWS_WITH(read_mesh_json(replaced(text, "\"block\": 0, \"face\": \"imin\"", "\"block\": 4294967296, \"face\": \"imin\"")),
                      Catch::Matchers::ContainsSubstring("range"));
    CHECK_THROWS_WITH(read_mesh_json(replaced(text, "\"block\": 0, \"face\": \"imin\"", "\"block\": -4294967296, \"face\": \"imin\"")),
                      Catch::Matchers::ContainsSubstring("range"));
    // a Dirichlet segment is verification only: never written
    CHECK_THROWS_AS(write_mesh_json(gfix::dirichlet_box(gfix::block(gfix::family(1), 4, 4))), MeshFormatError);
}

TEST_CASE("general mesh: the Aungier grid adapted node for node", "[turbolenza][general]")
{
    const Cascade c = make_cascade(fixture::synthetic_section(), 1.25, 30.0);
    const StandardGrid sg = build_grid(c, "very_coarse", 45.0 * fixture::kPi / 180.0, std::nullopt);
    const Grid& grid = sg.grid;
    const double c_ref = 0.5;
    const GeneralMesh mesh = mesh_from_grid(grid, c_ref);
    REQUIRE(mesh.blocks.size() == 1);
    const MeshBlock& B = mesh.blocks[0];
    CHECK(B.ni == grid.n_m());
    CHECK(B.nj == grid.n_eta());
    bool nodes = true;
    for (int i = 0; i < grid.n_m(); ++i)
        for (int j = 0; j < grid.n_eta(); ++j) {
            const std::size_t k = static_cast<std::size_t>(general_node_of(grid, i, j));
            nodes = nodes && B.x[k] == grid.node_m(i, j) && B.y[k] == grid.node_y(i, j) && B.m[k] == c_ref * grid.node_m(i, j) &&
                    B.r[k] == c_ref && B.b[k] == 1.0;
        }
    CHECK(nodes);
    CHECK(mesh.pitch == grid.pitch());
    const GeneralOperator op(mesh);
    // the walls are fixed (blade 0: psi 0, blade 1: psi 1), the periodic nodes one class with offset 1
    const int f = grid.first_blade_i, l = grid.last_blade_i, J = grid.n_eta() - 1;
    CHECK(op.equation(op.node(0, f, 0)) == NodeEquation::Fixed);
    CHECK(op.equation(op.node(0, l, J)) == NodeEquation::Fixed);
    CHECK(op.representative(op.node(0, 1, J)) == op.node(0, 1, 0));
    CHECK(op.images(op.node(0, 1, J)) == 1);
    CHECK(op.equation(op.node(0, 0, 3)) == NodeEquation::FlowAngle);
    CHECK(op.equation(op.node(0, f - 1, 0)) == NodeEquation::Pde);
    // the transition: the periodic node next to the wall has a one-sided difference along the side line
    CHECK(op.derivative(op.node(0, f - 1, 0), 0).size() == 3);
    CHECK(op.derivative(op.node(0, f, 0), 0).size() == 3);
    CHECK(op.derivative(op.node(0, f + 1, 0), 0).size() == 2);
    // the run rule in action: a field 0 on the walls and 1 off them has no difference along the wall at the
    // first blade station, nor along the side line at the last inlet station (central: -1/2 and +1/2)
    std::vector<double> off_wall(static_cast<std::size_t>(op.n_nodes()), 1.0);
    for (int g = 0; g < op.n_nodes(); ++g)
        if (op.wall(g)) off_wall[static_cast<std::size_t>(g)] = 0.0;
    CHECK(op.apply(op.derivative(op.node(0, f, 0), 0), off_wall, 0.0) == 0.0);
    CHECK(op.apply(op.derivative(op.node(0, l, J), 0), off_wall, 0.0) == 0.0);
    CHECK(op.apply(op.derivative(op.node(0, f - 1, 0), 0), off_wall, 0.0) == 0.0);
    // the transition classes get a closure face across the blade thickness, and their control volume closes
    const auto gcl = op.gcl_residual();
    for (int i : {f - 1, l + 1}) {
        const int r = op.representative(op.node(0, i, 0));
        int closures = 0;
        for (const FaceGeometry& fc : op.faces(r)) closures += fc.closure ? 1 : 0;
        CHECK(closures == 1);
        CHECK(std::abs(gcl[static_cast<std::size_t>(r)][0]) + std::abs(gcl[static_cast<std::size_t>(r)][1]) <= 1e-15);
    }
    // unknowns: every node but the walls, one per periodic pair
    const int periodic_stations = grid.n_m() - (l - f + 1);
    CHECK(op.n_unknowns() == grid.n_m() * grid.n_eta() - 2 * (l - f + 1) - periodic_stations);
    Grid touching = grid;
    touching.edge_nodes = EdgeNodes::Touching;
    CHECK_THROWS_AS(mesh_from_grid(touching, c_ref), MeshError);
    CHECK_THROWS_AS(mesh_from_grid(grid, -1.0), MeshError);
}

TEST_CASE("general mesh: the example of docs/turbolenza/mesh_format.md reads", "[turbolenza][general]")
{
    std::ifstream in(std::string(BLADEGEN_SOURCE_DIR) + "/docs/turbolenza/mesh_format.md", std::ios::binary);
    REQUIRE(in);
    std::ostringstream s;
    s << in.rdbuf();
    const std::string doc = s.str();
    const auto section = doc.find("## 9.");
    REQUIRE(section != std::string::npos);
    const auto begin = doc.find("```json\n", section);
    const auto end = doc.find("```", begin + 8);
    REQUIRE((begin != std::string::npos && end != std::string::npos));
    const GeneralMesh mesh = read_mesh_json(doc.substr(begin + 8, end - begin - 8));
    CHECK(mesh.blocks.size() == 1);
    CHECK(mesh.connections.size() == 1);
    const GeneralOperator op(mesh);
    CHECK(op.n_unknowns() == 6);  // 9 nodes, the 3 upper ones the images of the lower ones
}
