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
///  * the residual at a node P is the balance of the four face fluxes, (F_E - F_W) + (F_N - F_S) = (J f)_P;
///  * F through the face between P and its neighbour Q in index direction n (t the other one) is
///        F = K_nn (psi_Q - psi_P) + K_nt (psi_t(P) + psi_t(Q)) / 2,
///    with K at the face from face metrics: x_n = x_Q - x_P (the compact difference across the face) and
///    x_t = (x_t(P) + x_t(Q)) / 2 (the nodal tangential differences averaged), kappa_f = (kappa_P +
///    kappa_Q) / 2; at nodes whose eight neighbours all have central differences this is the nine-point
///    stencil P: -(pE + pW + pN + pS), E: pE + cN - cS, W: pW - cN + cS, N: pN + cE - cW, S: pS - cE + cW,
///    NE: cE + cN, NW: -cW - cN, SE: -cE - cS, SW: cW + cS, p = K_nn of the face, c = K_nt / 4;
///  * nodal differences: central (Q+ - Q-) / 2 when both neighbours exist; otherwise second-order
///    one-sided (-3 P + 4 Q1 - Q2) / 2 into the block (Aungier Eqs. 5-40 to 5-43).  Along a boundary face a
///    difference does not cross between a wall node and a non-wall node (Aungier's runs, Metrics.hpp):
///    it is one-sided on both sides of such a transition.  Neighbours across periodic and interface
///    connections are the nodes of the other side (their theta and psi shifted by the image);
///  * the face fluxes telescope: a face shared by two blocks carries one flux, which leaves one control
///    volume and enters the other (conservation across interfaces and the periodic cut); the discrete
///    metric identities d_zeta (y_eta) - d_eta (y_zeta) = 0 hold to round-off wherever all eight
///    neighbours have central differences, so a uniform flow (psi linear in x, y, kappa constant) is an
///    exact discrete solution there.
///
/// Boundary relations (one equation per unknown):
///  * blade walls psi = blade (the index, Eqs. 5-27, 5-28); Dirichlet segments psi given (verification);
///  * periodic and interface nodes are one unknown per node class: psi(node) = u(class) + images (Eq. 5-29,
///    images the pitches between the node and the class representative).  The class's PDE is the flux
///    balance of its whole control volume, each member contributing the part in its own block (ours):
///    its faces into the block, those along a block face it lies on with half their flux, and its share of
///    the source (1/2 per block face it lies on).  For a regular periodic or interface node the parts add
///    up to the nine-point equation; where the two images see different neighbours (the transition
///    stations of Aungier's H-grid, next to the leading and trailing edge) the equation averages the
///    suction- and the pressure-side halves;
///  * inlet and outlet: the flow angle T_f = tan(beta_f) prescribed, Psi_x + T_f Psi_y = 0 (from
///    tan(beta) = W_theta / W_m = -Psi_x / Psi_y), multiplied by J:
///        (y_eta - T_f x_eta) psi_zeta + (T_f x_zeta - y_zeta) psi_eta = 0
///    with the nodal (one-sided normal) differences; on Aungier's aligned grid (x_eta = 0, y_zeta / x_zeta =
///    T_f at the far station) it is psi_m = 0, Eq. 5-33.
///
/// Velocity recovery: psi_x = (y_eta psi_zeta - y_zeta psi_eta) / J, psi_y = (x_zeta psi_eta - x_eta psi_zeta) / J
/// from the nodal differences, at every node, each class taking its representative's values (Eqs. 5-30 to
/// 5-32 exactly).  Refusals (std::invalid_argument types): an invalid mesh (MeshError), a PDE node
/// without four neighbours, a node without a difference in a direction, a face or nodal Jacobian <= 0,
/// inputs of the wrong size or not finite (GeneralOperatorError).
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

/// The geometry of one face of a PDE node: F = kappa_f (G_nn dpsi + G_nt psi_t,f).
struct FaceGeometry
{
    int axis = 0, sign = 1;  ///< n = axis, the face toward sign
    double weight = 1.0;     ///< 1/2 along a block face the node lies on (its block's half of the face)
    NodeRef q;               ///< the neighbour across the face (always in the node's block)
    double x_n = 0.0, y_n = 0.0, x_t = 0.0, y_t = 0.0, J = 0.0;
    double G_nn = 0.0, G_nt = 0.0;  ///< K / kappa_f: K_nn = kappa_f G_nn, K_nt = kappa_f G_nt
    Stencil psi_t;           ///< (psi_t(P) + psi_t(Q)) / 2 in P's frame
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

    /// The in-block faces of a node of a PDE class (empty for other nodes), in the order W, E, S, N.
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
    /// The discrete metric identities at every PDE node (0 elsewhere): {sum of s sigma y_t,f, sum of
    /// -s sigma x_t,f} over the four faces = d_zeta (J zeta_x) + d_eta (J eta_x), likewise for y.
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
