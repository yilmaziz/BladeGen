/// @file turbolenza/GeneralMesh.hpp
/// @brief A structured multi-block mesh of one blade-to-blade stream sheet in the conformal plane, the
///        input of TURBOlenza's general curvilinear potential-flow operator (GeneralOperator.hpp).
///
/// Not a port of v1 (v1 knows only Aungier's sheared H-grid, GridData.hpp).  The file format of this
/// data is docs/turbolenza/mesh_format.md (MeshFormat.hpp); the vocabulary follows ENIGMAgrid's
/// core.block / core.grid (faces IMin, IMax, JMin, JMax; boundary segments on node ranges of a face;
/// node-to-node connections, periodic by a translation in theta), read from ENIGMAgrid's headers.
///
/// Coordinates (the user's derivation; Aungier 2003 Sec. 5.2, Eqs. 5-1, 5-2, 5-9 for the stream sheet):
///  * the stream sheet of revolution with meridional arc length m [metres] and radius r(m) [metres];
///    the conformal coordinates x = m' = int dm / r (dimensionless) and y = theta [rad], in which
///    ds^2 = r^2 (dx^2 + dy^2), so angles in (x, y) are the physical flow angles;
///  * a planar cascade is the sheet r = const = c_ref: x = m / c_ref, y = (tangential distance) / c_ref
///    (TURBOlenza's chord-normalized (m, y) of GridData.hpp, with r = c_ref the true chord in metres);
///  * b(m): the stream-sheet thickness relative to a reference thickness (dimensionless).
///
/// Blocks: n_i x n_j nodes, node (i, j) stored at i + n_i j (i fastest: ENIGMAgrid's and CGNS's
/// order).  Every block is right-handed in the (x, y) plane (ENIGMAgrid ENG-334: i along +m', j along
/// +theta for H blocks; i clockwise round a body, j outward, for O blocks): every cell corner Jacobian
/// (x_i y_j - x_j y_i of the cell's edges) is positive, or validate_mesh refuses the mesh.
///
/// Faces and segments: the face IMin (i = 0) and IMax (i = n_i - 1) are indexed along j, JMin and JMax
/// along i.  A segment is the node range [begin, end] (begin < end, both included) of one face with one
/// kind.  Two segments of a face may share their end node, never a cell; every node of every face lies
/// on at least one segment.  A cell between two segments that share no node (Aungier's H-grid: the
/// last inlet station and the first blade station) is a transition cell and belongs to neither.
///  * Blade: a blade wall, psi = blade (the blade index: 0 the reference blade, 1 its +pitch image, ...;
///    psi in units of the passage mass flow), with its surface (suction or pressure);
///  * Inlet, Outlet: the far boundaries, where the flow angle is prescribed;
///  * Periodic: paired by a Connection with another Periodic segment one pitch away (image +-1);
///  * Interface: paired by a Connection with another Interface segment, node to node, coinciding;
///  * Dirichlet: psi prescribed node by node by the caller (verification only; not in the file format).
///
/// Connections (ENIGMAgrid core.grid Patch, 2-D): node k of the first segment (begin + k) is paired with
/// node k of the second (begin + k, or end - k when reversed).  A step out of the first segment's block
/// through its face is a step into the second's block (the index map is the 2-D signed permutation that
/// maps the first face's outward normal to the second face's inward normal and the first face's
/// tangent to the second's, reversed or not); it must keep the orientation (determinant +1).  Periodic:
/// theta(second) = theta(first) + image * pitch and psi(second) = psi(first) + image (Aungier 2003
/// Eq. 5-29), image = +1 or -1; Interface: image = 0, the nodes coincide.  Every Periodic and Interface
/// segment is in exactly one connection.  Nonmatching interfaces are refused.
///
/// Tags: the leading- and trailing-edge node of each blade surface (blade, surface): the nodes the Kutta
/// condition compares (Aungier 2003 Eq. 5-34) and post-processing starts and ends at.  Each Blade
/// segment's (blade, surface) has exactly one leading and one trailing tag on one of its nodes.
#pragma once

#include "turbolenza/GridData.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace BladeGen::Turbolenza {

/// A mesh that breaks the rules above.
class MeshError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

enum class BlockFace : int { IMin = 0, IMax = 1, JMin = 2, JMax = 3 };
enum class SegmentKind : int { Blade = 0, Inlet = 1, Outlet = 2, Periodic = 3, Interface = 4, Dirichlet = 5 };
enum class EdgeKind : int { Leading = 0, Trailing = 1 };

[[nodiscard]] const char* block_face_name(BlockFace f) noexcept;     ///< "imin", "imax", "jmin", "jmax"
[[nodiscard]] const char* segment_kind_name(SegmentKind k) noexcept;  ///< "blade", "inlet", ...
[[nodiscard]] const char* edge_kind_name(EdgeKind e) noexcept;        ///< "leading", "trailing"
[[nodiscard]] const char* mesh_surface_name(Surface s) noexcept;      ///< "suction", "pressure"

/// One structured block: node fields of n_i n_j values, i fastest.
struct MeshBlock
{
    std::string name;
    int ni = 0, nj = 0;
    std::vector<double> x;  ///< m' (dimensionless)
    std::vector<double> y;  ///< theta [rad]
    std::vector<double> m;  ///< meridional arc length [m] (reports; r and b are its functions)
    std::vector<double> r;  ///< radius of the stream sheet [m], > 0
    std::vector<double> b;  ///< relative stream-sheet thickness, > 0

    [[nodiscard]] std::size_t at(int i, int j) const noexcept
    {
        return static_cast<std::size_t>(i) + static_cast<std::size_t>(ni) * static_cast<std::size_t>(j);
    }
    [[nodiscard]] int nodes() const noexcept { return ni * nj; }
};

/// A node range of one block face with one kind (see the file comment).
struct FaceSegment
{
    int block = 0;
    BlockFace face = BlockFace::JMin;
    int begin = 0, end = 0;          ///< along the face: j on IMin / IMax, i on JMin / JMax
    SegmentKind kind = SegmentKind::Blade;
    int blade = 0;                   ///< Blade: the blade index, psi = blade
    Surface surface = Surface::Suction;  ///< Blade

    [[nodiscard]] int length() const noexcept { return end - begin + 1; }
};

/// A node-to-node pairing of two segments (see the file comment).
struct Connection
{
    int first = 0, second = 0;  ///< segment indices
    bool reversed = false;
    int image = 0;              ///< Periodic: +1 or -1; Interface: 0
};

/// The leading- or trailing-edge node of one blade surface.
struct EdgeTag
{
    EdgeKind edge = EdgeKind::Trailing;
    int blade = 0;
    Surface surface = Surface::Suction;
    int block = 0, i = 0, j = 0;
};

/// The mesh of one blade passage (or any multi-block domain of one stream sheet).
struct GeneralMesh
{
    std::string description;
    double pitch = 0.0;  ///< the periodic shift in theta [rad] (2 pi / N, or s / c_ref for a planar cascade)
    std::vector<MeshBlock> blocks;
    std::vector<FaceSegment> segments;
    std::vector<Connection> connections;
    std::vector<EdgeTag> tags;
};

/// Tolerances of validate_mesh.
struct MeshCheckOptions
{
    double coincidence = 1e-9;  ///< paired nodes: |dx|, |dy - image pitch| and the relative |dm|, |dr|, |db|
};

/// The number of nodes along a face (n_j for IMin / IMax, n_i for JMin / JMax).
[[nodiscard]] int face_length(const MeshBlock& block, BlockFace face) noexcept;
/// The (i, j) of the node at position t along a face.
[[nodiscard]] std::array<int, 2> face_node(const MeshBlock& block, BlockFace face, int t) noexcept;

/// The 2-D index map of a connection from `from` to `to` (row major, s_to = T s_from for an index step
/// s): the outward normal of `from` goes to the inward normal of `to`, the tangent of `from` to the
/// tangent of `to` (negated when reversed).  Its determinant is +1 or -1 (validate_mesh refuses -1).
[[nodiscard]] std::array<int, 4> face_transform(BlockFace from, BlockFace to, bool reversed) noexcept;

/// Checks every rule of the file comment; MeshError naming the first one broken.
void validate_mesh(const GeneralMesh& mesh, const MeshCheckOptions& options = {});

}  // namespace BladeGen::Turbolenza
