/// @file turbolenza/GeneralAdapter.hpp
/// @brief TURBOlenza's own Aungier grid (GridData.hpp) as a one-block GeneralMesh on the same nodes, so
///        that the Aungier solver and the general operator can be compared node by node.
///
/// The planar cascade is the stream sheet r = c_ref (GeneralMesh.hpp): m' = m* (the chord-normalized m of
/// the grid), theta = y* (its chord-normalized tangential coordinate), m = c_ref m*, r = c_ref, b = 1,
/// pitch = grid.pitch() (in chords).  Block 0 has n_i = n_m (i along m) and n_j = n_eta (j along eta), so
/// general node i + n_m j is Aungier node (i, j).  Segments (Aungier 2003 Fig. 5-3, Eqs. 5-27 to 5-33):
///   0 jmin periodic [0, f - 1]    1 jmin blade 0, the lower wall's surface [f, l]    2 jmin periodic [l + 1, n_m - 1]
///   3 jmax periodic [0, f - 1]    4 jmax blade 1, the upper wall's surface [f, l]    5 jmax periodic [l + 1, n_m - 1]
///   6 imin inlet [0, n_eta - 1]   7 imax outlet [0, n_eta - 1]
/// connections 0-3 and 2-5 with image +1, the leading-edge tags at station f and the trailing-edge tags at
/// station l, on both walls.  The transition cells f - 1 .. f and l .. l + 1 of the two side faces lie
/// between a periodic and a wall segment.  Refused (MeshError): the EdgeNodes::Touching layout, a wall
/// whose nodes do not all lie on one surface, fewer than two periodic stations on a side.
#pragma once

#include "turbolenza/GeneralMesh.hpp"
#include "turbolenza/GridData.hpp"

namespace BladeGen::Turbolenza {

[[nodiscard]] GeneralMesh mesh_from_grid(const Grid& grid, double c_ref);

/// The general node of Aungier node (i, j) in mesh_from_grid's mesh.
[[nodiscard]] inline int general_node_of(const Grid& grid, int i, int j) noexcept { return i + grid.n_m() * j; }

}  // namespace BladeGen::Turbolenza
