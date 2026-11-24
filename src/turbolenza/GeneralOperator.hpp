/// @file turbolenza/GeneralOperator.hpp
/// @brief The stream-function potential-flow operator of one stream sheet on any structured multi-block
///        mesh (GeneralMesh.hpp): the conservative (flux) form, its boundary relations, the assembled
///        sparse system and the velocity recovery.  Not a port of v1: v1 (and Aungier) discretize the
///        expanded form on the sheared H-grid only (Discretization.hpp).
///
/// The equation (the user's derivation, checked step by step; Aungier 2003 Sec. 5.3, Eqs. 5-10 to 5-21):
///  * continuity on the sheet, d/dm (r rho b W_m) + d/dtheta (rho b W_theta) = 0, is satisfied by the
///    stream function Psi with Psi_theta = r rho b W_m, Psi_m = -rho b W_theta; Psi = mdot' psi, mdot' the
///    passage mass flow per unit b, psi dimensionless (0 on blade 0, +1 per pitch: Eqs. 5-16, 5-17, 5-29);
///  * zero absolute vorticity normal to the sheet (Eq. 5-11, C = W + omega r e_theta), on the conformal
///    plane x = m' (dx = dm / r), y = theta:
///        d/dx (kappa Psi_x) + d/dy (kappa Psi_y) = omega d(r^2)/dx,    kappa = 1 / (rho b),
///        W_m = (kappa / r) Psi_y,   W_theta = -(kappa / r) Psi_x;
///    omega is the angular speed in the +theta sense (U = omega r e_theta, Aungier's Eq. 5-11 frame);
///  * on a block x(zeta, eta), y(zeta, eta) (computational coordinates = the node indices, h = k = 1):
///        d_zeta (K11 Psi_zeta + K12 Psi_eta) + d_eta (K12 Psi_zeta + K22 Psi_eta) = J f,
///        J = x_zeta y_eta - x_eta y_zeta,  K = kappa J g^ab = (kappa / J) [[x_eta^2 + y_eta^2, -(x_zeta x_eta
///        + y_zeta y_eta)], [., x_zeta^2 + y_zeta^2]],  f = omega d(r^2)/dx,  J f = omega (y_eta (r^2)_zeta -
///        y_zeta (r^2)_eta);
///  * on Aungier's mapping x = m'(m), y = theta0(m) + eta Theta(m) with (m, eta) computational, J = S / r^2,
///    J g = [[S, -T], [-T, (1 + T^2) / S]], T = tan(beta) (Eq. 5-3), S = r Theta (Eq. 5-9); mdot' kappa J g =
///    [[C, -B], [-B, A]] (Eqs. 5-22 to 5-24), the flux form is Eq. 5-20 and J f = 2 S omega sin(phi).
///
/// The discretization (ours; the user's scheme, the nine-point table re-derived in the tests):
///  * the residual of a class is the outward flux through the boundary of its control volume, the dual
///    cell round the node: from the centroid of each incident cell to the midpoint of the next edge to the
///    next centroid, counterclockwise.  Sum of the face fluxes = (J f)_P (the source shared by the members);
///  * the face toward a neighbour Q runs between the centroids C1, C2 of the two cells on the edge PQ:
///        F = kappa_f (|t|^2 dpsi - (d . t) dpsi_t) / (d x t),  d = Q - P, t = C2 - C1,
///        dpsi = psi_Q - psi_P, dpsi_t = psi(C2) - psi(C1), psi(C) the mean of the cell's four corners,
///    kappa_f = (kappa_P + kappa_Q) / 2: the gradient at the face from its two directional differences.
///    For the face E this is F_E = K11_E (psi_E - psi_P) + K12_E (psi_NE - psi_SE + psi_N - psi_S) / 4
///    with K11 = kappa (x_eta^2 + y_eta^2) / J, K12 = -kappa (x_zeta x_eta + y_zeta y_eta) / J of the face
///    metrics x_zeta = x_E - x_P, x_eta = (x_N + x_NE - x_S - x_SE) / 4 (the central nodal differences of P
///    and E averaged), and for a node with four neighbours in its block it is the nine-point stencil
///    P: -(pE + pW + pN + pS), E: pE + cN - cS, W: pW - cN + cS, N: pN + cE - cW, S: pS - cE + cW,
///    NE: cE + cN, NW: -cW - cN, SE: -cE - cS, SW: cW + cS, p = K_nn of the face, c = K_nt / 4;
///  * a periodic or interface class collects the cells of all its members (each in its own block, shifted
///    into the representative's frame by its images): one control volume across blocks and the periodic
///    cut.  Where the half-edges of two members end at different nodes (the last inlet and the first exit
///    station of Aungier's H-grid: the two images of the periodic node meet the suction wall node (f, 0)
///    and the pressure wall node (f, J) one pitch down, a blade thickness h_f apart), one closure face runs
///    across the gap, from the last centroid of one member to the first of the other, toward the mean of
///    the two wall nodes (psi their mean).  Without it the control volume would leave a gap of h_f / 2 in
///    its boundary and the discrete metric identities would fail there by h_f / 2 (measured, see the
///    tests); with it every control volume is a closed polygon, the metric identities
///    d_zeta (J zeta_x) + d_eta (J eta_x) = 0 (likewise for y) hold to round-off at every PDE class, and
///    a uniform flow (psi linear in x, y, kappa constant) is an exact discrete solution of every PDE row;
///  * every face is computed once per class from the same nodes from both sides, so the fluxes telescope:
///    conservation across interfaces and the periodic cut;
///  * nodal differences (flow-angle rows, velocity recovery, sources): central (Q+ - Q-) / 2 when both
///    neighbours exist; otherwise second-order one-sided (-3 P + 4 Q1 - Q2) / 2 into the block (Aungier
///    Eqs. 5-40 to 5-43).  Along a boundary face a nodal difference does not cross between a wall node and
///    a non-wall node (Aungier's runs, Metrics.hpp): it is one-sided on both sides of such a transition.
///    Neighbours across periodic and interface connections are the nodes of the other side (their theta
///    and psi shifted by the image).
/// Accuracy at H-grid transitions (measured on a manufactured problem with Aungier's layout,
/// test_general_operator.cpp): psi second order with zero thickness at the edge stations, with a sharp edge
/// (thickness ~ h at the first station) and with a C-infinity bump; first order with a round edge
/// (thickness ~ sqrt(h)), whose velocity next to the edge does not converge in the max norm.  Without the
/// closure face (0ede48c: the two halves only) the GCL failed there by h_f / 2 and the orders were 1 and
/// 0.57.  Aungier's own transition rows (v1), re-implemented by the independent review of 0ede48c, are
/// first order with a sharp edge and about half order with a round one; on the manufactured problem both
/// schemes converge to the exact solution (two slow schemes with different error constants, not
/// different limits).
///
/// Boundary relations (one equation per unknown):
///  * blade walls psi = blade (the index, Eqs. 5-27, 5-28); Dirichlet segments psi given (verification);
///  * periodic and interface nodes are one unknown per node class: psi(node) = u(class) + images (Eq. 5-29,
///    images the pitches between the node and the class representative); the class's PDE is the flux
///    balance of its control volume (above), its source the members' (J f) weighted 1/2 per block face
///    each lies on;
///  * inlet and outlet: the flow angle T_f = tan(beta_f) prescribed, Psi_x + T_f Psi_y = 0 (from
///    tan(beta) = W_theta / W_m = -Psi_x / Psi_y), multiplied by J:
///        (y_eta - T_f x_eta) psi_zeta + (T_f x_zeta - y_zeta) psi_eta = 0
///    with the nodal (one-sided normal) differences; on Aungier's aligned grid (x_eta = 0, y_zeta / x_zeta =
///    T_f at the far station) it is psi_m = 0, Eq. 5-33.
///
/// Velocity recovery: psi_x = (y_eta psi_zeta - y_zeta psi_eta) / J, psi_y = (x_zeta psi_eta - x_eta psi_zeta) / J
/// from the nodal differences, at every node, each class taking its representative's values (Eqs. 5-30 to
/// 5-32 exactly).  Refusals (std::invalid_argument types): an invalid mesh (MeshError); a node without a
/// nodal difference in a direction; a nodal or face Jacobian <= 0; a control volume that does not close or
/// meets one neighbour on more than two half-edges (a singular vertex); inputs of the wrong size or not
/// finite (GeneralOperatorError).
#pragma once

#include "turbolenza/GeneralMesh.hpp"
#include "turbolenza/MatrixAssembly.hpp"

#include <array>
#include <functional>
#include <optional>
#include <stdexcept>
#include <vector>

namespace BladeGen::Turbolenza {

class GeneralOperatorError : public std::invalid_argument
{
public:
    using std::invalid_argument::invalid_argument;
};

/// A node seen from an origin node: s_node = T s_origin for an index step s (T row major), and the
/// pitches between them: theta_origin = theta_node - images pitch, psi_origin = psi_node - images.
struct NodeRef
{
    int node = -1;
    std::array<int, 4> T{1, 0, 0, 1};
    int images = 0;
};

/// sum w (value(node) - images jump): jump = pitch for theta, 1 for psi, 0 for a field without one.
struct StencilEntry
{
    int node = 0;
    double w = 0.0;
    int images = 0;
};
using Stencil = std::vector<StencilEntry>;

enum class NodeEquation { Fixed, FlowAngle, Pde };

/// Nodal metric (own block frame): x_zeta, x_eta, y_zeta, y_eta and J = x_zeta y_eta - x_eta y_zeta.
struct NodalMetric
{
    double x_z = 0.0, x_e = 0.0, y_z = 0.0, y_e = 0.0, J = 0.0;
};

/// One face of the control volume of a PDE class (see the file comment), in the representative's frame:
/// the outward flux is kappa_f (G_nn dpsi + G_nt psi_t), kappa_f = sum w kappa(node).
struct FaceGeometry
{
    int axis = -1, sign = 0;          ///< the representative's index direction toward the neighbour (-1: none)
    bool closure = false;             ///< a face closing the gap between two half-edges (H-grid transition)
    std::vector<int> neighbours;      ///< the node across the face (two for a closure face)
    double x_n = 0.0, y_n = 0.0;      ///< d = Q - P (Q the mean of the two nodes of a closure face)
    double x_t = 0.0, y_t = 0.0;      ///< t = C2 - C1 between boundary corners (cell centroids), counterclockwise
    double J = 0.0;                   ///< d x t > 0
    double G_nn = 0.0, G_nt = 0.0;    ///< |t|^2 / J, -(d . t) / J
    Stencil dpsi;                     ///< psi(Q) - psi(P)
    Stencil psi_t;                    ///< psi(C2) - psi(C1), psi(C) the mean of the cell's corners
    std::vector<std::pair<int, double>> kappa;  ///< kappa_f = (kappa_P + kappa_Q) / 2
};

/// The per-pass inputs of the operator, one value per node (n_nodes()).
struct OperatorFields
{
    std::vector<double> coefficient;     ///< mdot' kappa = mdot' / (rho b): K = coefficient J g^ab
    std::vector<double> source;          ///< (J f) at the node; empty: 0
    std::vector<double> flow_angle_tan;  ///< T_f, read at FlowAngle nodes; empty: none needed
    std::function<double(int)> dirichlet;  ///< psi at the nodes of Dirichlet segments
};

struct GeneralSystem
{
    SparseMatrix matrix;
    std::vector<double> rhs;
};

/// psi_x, psi_y at every node.
struct StreamGradient
{
    std::vector<double> psi_x, psi_y;
};

class GeneralOperator
{
public:
    /// Validates the mesh (validate_mesh) and builds the node classes, differences and face geometry.
    explicit GeneralOperator(GeneralMesh mesh, const MeshCheckOptions& options = {});

    [[nodiscard]] const GeneralMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] int n_nodes() const noexcept { return n_nodes_; }
    [[nodiscard]] int node(int block, int i, int j) const;
    [[nodiscard]] std::array<int, 3> location(int g) const;  ///< (block, i, j)
    [[nodiscard]] double x(int g) const { return field(g, &MeshBlock::x); }
    [[nodiscard]] double y(int g) const { return field(g, &MeshBlock::y); }
    [[nodiscard]] double m(int g) const { return field(g, &MeshBlock::m); }
    [[nodiscard]] double r(int g) const { return field(g, &MeshBlock::r); }
    [[nodiscard]] double b(int g) const { return field(g, &MeshBlock::b); }

    /// The neighbour one index step away (axis 0: zeta = i, 1: eta = j; sign +-1), across a connection
    /// when the step leaves the block; empty at a boundary without one.
    [[nodiscard]] std::optional<NodeRef> step(int g, int axis, int sign) const;
    /// The difference d/d(axis) at node g in its own block frame (see the file comment).
    [[nodiscard]] const Stencil& derivative(int g, int axis) const { return deriv_[static_cast<std::size_t>(2 * g + axis)]; }
    [[nodiscard]] const NodalMetric& metric(int g) const { return metric_[static_cast<std::size_t>(g)]; }
    [[nodiscard]] bool wall(int g) const { return wall_[static_cast<std::size_t>(g)] != 0; }

    // node classes
    [[nodiscard]] int representative(int g) const { return rep_[static_cast<std::size_t>(g)]; }
    [[nodiscard]] int images(int g) const { return images_[static_cast<std::size_t>(g)]; }
    [[nodiscard]] NodeEquation equation(int g) const;
    /// The node whose relation the class of g gets (its representative for the PDE, its inlet / outlet
    /// member for the flow angle, its wall member when fixed).
    [[nodiscard]] int row_node(int g) const;
    /// The unknown of g's class, -1 when fixed.
    [[nodiscard]] int unknown(int g) const { return unknown_[static_cast<std::size_t>(rep_[static_cast<std::size_t>(g)])]; }
    [[nodiscard]] int n_unknowns() const noexcept { return n_unknowns_; }
    /// The nodes of the segments of `kind` (each once, in segment order).
    [[nodiscard]] std::vector<int> segment_nodes(SegmentKind kind) const;
    /// The nodes of segment s, in its order.
    [[nodiscard]] std::vector<int> segment_nodes(int s) const;

    /// The faces of the control volume of a PDE class, at its representative (empty for other nodes).
    [[nodiscard]] const std::vector<FaceGeometry>& faces(int g) const;
    /// The members of the class whose representative is r (r itself first is not guaranteed).
    [[nodiscard]] const std::vector<int>& members(int r) const { return members_[static_cast<std::size_t>(r)]; }
    /// The share of a member's (J f) in its class's equation (1, 1/2 on a block face, 1/4 at a corner).
    [[nodiscard]] double source_weight(int g) const { return source_weight_[static_cast<std::size_t>(g)]; }

    /// The assembled system in the unknowns, each row divided by the coefficient of its own unknown (the
    /// analogue of Aungier's Eq. 5-44 normalization; GeneralOperatorError for broken inputs or a zero one).
    [[nodiscard]] GeneralSystem assemble(const OperatorFields& fields) const;
    /// psi at every node from the unknowns (fixed classes from the walls and fields.dirichlet).
    [[nodiscard]] std::vector<double> expand(const std::vector<double>& unknowns, const OperatorFields& fields) const;
    /// The residual of every class's relation at psi (one per node: the relation of its class; the
    /// fixed ones psi - value).  The discrete equations are residual = 0.
    [[nodiscard]] std::vector<double> residual(const std::vector<double>& psi, const OperatorFields& fields) const;
    /// psi_x, psi_y at every node (each class its representative's).
    [[nodiscard]] StreamGradient gradient(const std::vector<double>& psi) const;
    /// (J f) of the rotation term at every node: omega (y_eta (r^2)_zeta - y_zeta (r^2)_eta).
    [[nodiscard]] std::vector<double> rotation_source(double omega) const;
    /// The discrete metric identities of every PDE class (each member gets its class's; 0 elsewhere): the
    /// sum of the outward normals (t_y, -t_x) over the faces of the control volume, = d_zeta (J zeta_x) +
    /// d_eta (J eta_x) and likewise for y; zero to round-off because the control volume is closed.
    [[nodiscard]] std::vector<std::array<double, 2>> gcl_residual() const;
    /// Applies a stencil to a node field (jump: the value of one image: pitch for theta, 1 for psi).
    [[nodiscard]] double apply(const Stencil& s, const std::vector<double>& values, double jump) const;

private:
    double field(int g, std::vector<double> MeshBlock::*f) const;

    GeneralMesh mesh_;
    int n_nodes_ = 0, n_unknowns_ = 0;
    std::vector<int> offset_;  // first global node of each block
    std::vector<std::array<std::vector<std::vector<int>>, 4>> face_segs_;  // [block][face][t]: segments on node t
    std::vector<int> conn_of_;            // per segment: its connection, -1
    std::vector<char> inlet_, outlet_;
    std::vector<char> wall_, wall_face_;  // wall node; bit f: on a wall segment of face f
    std::vector<int> rep_, images_, unknown_, row_node_;
    std::vector<NodeEquation> equation_;  // per representative
    std::vector<double> fixed_;            // per representative: psi of the representative (blade walls)
    std::vector<char> dirichlet_;          // per representative: fixed by a Dirichlet segment (row_node_ the member)
    std::vector<Stencil> deriv_;
    std::vector<NodalMetric> metric_;
    std::vector<std::vector<FaceGeometry>> faces_;
    std::vector<std::vector<int>> members_;
    std::vector<double> source_weight_;
};

/// The contravariant tensor K = kappa J g^ab of the metric (x_zeta, x_eta, y_zeta, y_eta) and kappa:
/// {K11, K12, K22}.
[[nodiscard]] std::array<double, 3> contravariant_tensor(double x_z, double x_e, double y_z, double y_e, double kappa);

}  // namespace BladeGen::Turbolenza
