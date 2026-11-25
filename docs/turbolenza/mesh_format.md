# The `turbolenza.mesh` file format, version 1

This is the contract between a mesh producer (ENIGMAgrid's writer, TURBOlenza's own adapter of its
Aungier grid) and TURBOlenza's general curvilinear potential-flow solver (`src/turbolenza/General*`).
Reader and writer: `src/turbolenza/MeshFormat.hpp`. In-memory form: `GeneralMesh` (`GeneralMesh.hpp`).
The reader is strict: it refuses anything it does not understand, and anything inconsistent, with a
message that names the offending item. It never repairs a mesh.

## 1. What the file holds

One blade-to-blade stream sheet, meshed by one or more structured 2-D blocks in the **conformal plane**:

| quantity | key | unit | meaning |
|---|---|---|---|
| x | `m_prime` | 1 (dimensionless) | m' = ∫ dm / r along the stream sheet |
| y | `theta` | rad | the circumferential angle θ |
| m | `m` | m | meridional arc length (for reports; r and b are functions of it) |
| r | `r` | m | radius of the stream sheet, > 0 |
| b | `b` | 1 (relative) | stream-sheet thickness relative to a reference thickness, > 0 |

In (m', θ) the metric is ds² = r²(dm'² + dθ²), so angles measured in this plane are the physical flow
angles. A **planar cascade** is the sheet r = const = c_ref (the chord in metres): m' = m/c_ref and θ =
(tangential distance)/c_ref, the pitch in θ is s/c_ref. TURBOlenza's own Aungier grid is written this
way (`mesh_from_grid`).

Numbers are written with 17 significant digits (`{:.17g}`), so a double reads back bit for bit.

## 2. Top level

```json
{
  "format": "turbolenza.mesh",
  "version": 1,
  "units": {"m_prime": "1", "theta": "rad", "m": "m", "r": "m", "b": "1"},
  "description": "free text",
  "pitch": 0.6283185307179586,
  "blocks": [ ... ],
  "segments": [ ... ],
  "connections": [ ... ],
  "tags": [ ... ]
}
```

- `format` must be `"turbolenza.mesh"`, `version` must be `1`, `units` must be exactly the object above.
- `pitch` is the periodic shift in θ [rad] between a passage and its neighbour (2π/N for N blades).
- `description` is optional (default empty). No other key is allowed (at any level), and no key may
  appear twice in one object.
- Every integer (`version`, `ni`, `nj`, `block`, `range`, `segments`, `image`, `blade`, `node`) is a JSON
  integer within the range of a 32-bit signed int; `3.0` or `3e9` is refused.

## 3. Blocks

```json
{"name": "passage", "ni": 141, "nj": 29,
 "m_prime": [...], "theta": [...], "m": [...], "r": [...], "b": [...]}
```

- `name` is optional (default empty); every other key of a block is required.
- `ni × nj` nodes, `ni, nj ≥ 3`. Each array has `ni·nj` finite values, **i fastest**: node (i, j) is
  element `i + ni·j` (ENIGMAgrid's and CGNS's order, without ghost nodes).
- **Orientation** (ENIGMAgrid ENG-334): every block is right-handed in (m', θ). Every corner Jacobian of
  every cell, (x_i y_j − x_j y_i) of the cell's two edges at that corner, must be positive. For H blocks
  this is i along +m', j along +θ; for O blocks i runs clockwise round the body and j outward.
- Blocks are numbered by their position in the array, from 0.

## 4. Segments (boundary kinds of face parts)

```json
{"block": 0, "face": "jmin", "range": [50, 91], "kind": "blade", "blade": 0, "surface": "suction"}
```

- `face`: `imin` (i = 0), `imax` (i = ni − 1), `jmin` (j = 0), `jmax` (j = nj − 1).
- `range: [begin, end]`: node positions along the face, both included, `0 ≤ begin < end < length`, where
  the position along `imin`/`imax` is j and along `jmin`/`jmax` is i.
- Two segments of one face may share their end node, never a cell. Every node of every face must lie on
  at least one segment. A cell between two segments that share no node is a *transition cell* (for
  example Aungier's grid between the last inlet station and the first blade station) and belongs to
  neither.
- `kind`:
  - `blade`: a blade wall, ψ = `blade` (an integer: 0 the reference blade, 1 its neighbour one pitch
    toward +θ, −1 one pitch toward −θ). `surface` is `suction` or `pressure`. Only blade segments carry
    `blade` and `surface`.
  - `inlet`, `outlet`: the far boundaries. The solver prescribes the flow angle there
    (Ψ_x + tan β Ψ_y = 0; the inlet angle from the operating point, the outlet angle from the Kutta loop).
  - `periodic`: paired with another periodic segment one pitch away (section 5).
  - `interface`: paired node to node with another interface segment of the same or another block.

ψ is the stream function in units of the passage mass flow: it rises by 1 across one pitch.

## 5. Connections

```json
{"kind": "periodic",  "segments": [0, 3], "reversed": false, "image": 1}
{"kind": "interface", "segments": [8, 9], "reversed": true}
```

- `segments: [first, second]`: two segments of the connection's kind, with the same number of nodes.
  Node k of `first` (position `begin + k`) is paired with node k of `second`: position `begin + k`, or
  `end − k` when `reversed` is true.
- Every periodic and every interface segment is in exactly one connection. A connection is listed once;
  the reverse direction is implied.
- `periodic`: `image` is +1 or −1, and θ(second) = θ(first) + image·pitch, ψ(second) = ψ(first) + image.
  m', m, r and b of paired nodes are equal.
- `interface`: no `image`; the paired nodes coincide (m', θ, m, r, b).
- Coincidence is checked to 1e-9 (absolute in m' and θ, relative in m, r, b). **Nonmatching interfaces
  (different node counts, or nodes that do not coincide) are refused.**
- Orientation: a step out of `first`'s block through its face is a step into `second`'s block. The index
  map T (a 2-D signed permutation) sends the outward normal of `first`'s face to the inward normal of
  `second`'s face and the tangent of `first`'s face to the tangent of `second`'s (negated when
  `reversed`). T must have determinant +1; a connection that would fold one block over the other is
  refused.

## 6. Tags

```json
{"edge": "trailing", "blade": 0, "surface": "suction", "block": 0, "node": [91, 0]}
```

- The leading- and trailing-edge node of each blade surface. Each (blade, surface) that has a blade
  segment needs exactly one `leading` and one `trailing` tag, on a node of a blade segment of that blade
  and surface.
- The Kutta condition (Aungier 2003 Eq. 5-34) compares the relative speeds at the `suction` and the
  `pressure` trailing-edge tags; the solver needs exactly one of each. They are found by the tags, never
  by array position.
- On Aungier's own grid, where no node lies on the leading or trailing edge itself, the tags are the wall
  nodes of the first and last blade stations (Aungier 2003 pp. 89-90).

## 7. What the reader refuses

Schema (`MeshFormatError`): not JSON; a key twice in one object; another format name or version; other
units; a missing or unknown key; a wrong type; an integer out of range; an unknown face, kind, surface or edge; a connection whose kind differs
from its segments' kind; an `image` on an interface; `blade`/`surface` on a non-blade segment.

Consistency (`MeshError`, `validate_mesh`): a block smaller than 3 × 3; array sizes that differ from
`ni·nj`; a non-finite value; r ≤ 0 or b ≤ 0; a corner Jacobian ≤ 0 (folded or left-handed cell); a
segment outside its face or with `begin ≥ end`; two segments sharing a cell; a face node on no segment;
an unpaired (or twice paired) periodic or interface segment; a periodic image other than ±1; node counts
that differ; paired nodes that do not coincide; an orientation-reversing connection; a missing,
duplicated or misplaced leading- or trailing-edge tag.

## 8. Writing from ENIGMAgrid

ENIGMAgrid's `core::Block` / `core::Grid` map onto this format as follows (read from ENIGMAgrid's
`include/enigma/core/block.hpp` and `grid.hpp`):

| ENIGMAgrid | `turbolenza.mesh` |
|---|---|
| a 2-D block (nk = 1) in the section plane {m', θ, 0}, core nodes only | a block; `m_prime` = x, `theta` = y; add `m`, `r`, `b` |
| `Face::IMin` .. `JMax` | `imin` .. `jmax` |
| `BoundarySegment` range (two corners on the face) | `range` = the tangential index of the two corners |
| `BoundaryKind::Wall` on a blade | `blade` with the blade index (ψ) and its surface |
| `BoundaryKind::Inflow`, `Outflow` | `inlet`, `outlet` |
| `Patch` `Periodic`, `RigidTransform::translation({0, Δθ, 0})` | `periodic` segments on both faces, a connection receiver → `first`, donor → `second`, `image` = Δθ / pitch |
| `Patch` `OneToOne` | `interface` segments, a connection |
| `Patch::orientation` (CGNS Transform) | `reversed` = the sign with which the receiver's tangential direction maps to the donor's |

Notes for the writer:
- The 2-D orientation rule of ENG-334 (J > 0) is exactly this format's; ENIGMAgrid's patch rule "a step
  out of the receiver is a step into the donor" is the index map of section 5.
- ENIGMAgrid's H-grid splits the walls at the least and greatest m' points and lets the periodic and wall
  segments share that node; this is allowed (a shared end node; the node is a wall node).
- Give each blade surface its leading- and trailing-edge tag (the split points of the wall chains).
- For a planar cascade write r = c_ref and b = 1 at every node.

## 9. Example (a 3 × 3 block, both faces periodic, Dirichlet-free)

```json
{
  "format": "turbolenza.mesh",
  "version": 1,
  "units": {"m_prime": "1", "theta": "rad", "m": "m", "r": "m", "b": "1"},
  "description": "minimal",
  "pitch": 1,
  "blocks": [
    {"name": "b", "ni": 3, "nj": 3,
      "m_prime": [0, 0.5, 1, 0, 0.5, 1, 0, 0.5, 1],
      "theta": [0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1],
      "m": [0, 0.5, 1, 0, 0.5, 1, 0, 0.5, 1],
      "r": [1, 1, 1, 1, 1, 1, 1, 1, 1],
      "b": [1, 1, 1, 1, 1, 1, 1, 1, 1]}
  ],
  "segments": [
    {"block": 0, "face": "imin", "range": [0, 2], "kind": "inlet"},
    {"block": 0, "face": "imax", "range": [0, 2], "kind": "outlet"},
    {"block": 0, "face": "jmin", "range": [0, 2], "kind": "periodic"},
    {"block": 0, "face": "jmax", "range": [0, 2], "kind": "periodic"}
  ],
  "connections": [
    {"kind": "periodic", "segments": [2, 3], "reversed": false, "image": 1}
  ],
  "tags": []
}
```

(Without a blade the stream function is fixed only up to a constant: such a mesh is valid as a file,
but the solver needs a wall or a prescribed value to pin ψ.)
