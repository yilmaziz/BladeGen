# Changelog

## v1.0.4 — Stacking on a reference point (DLR variants 1 and 2, centre of gravity)

### New
- **Stacking of the sections** (Initial Meridional Configuration; Vos et al. 2025, Sec. 4). Two ways to fix
  where each section sits on its stream surface:
  - **prescribed LE / TE lines** (Sec. 4.2): what BladeGen did up to v1.0.3, and still the default. Each
    section's LE sits on the meridional LE line and its TE on the TE line, so the chord follows.
  - **reference point** (Sec. 4.1): each section's **centre of gravity**, **LE** or **TE** sits on a stacking
    line, and the chord is free. Two variants place the reference points:
    - **variant 1, same Z and θ**: every reference point at (Z_ref, θ_ref), so the stacking line is the
      radial line through that point. With the centre of gravity this is the radial CG stack.
    - **variant 2, closest point**: each layer's reference point is the point of its streamline closest to
      the previous layer's, at the same θ. It solves (C(m′) − Q)·dC/dm′ = 0, so the line follows the normals
      of the stream surfaces. In the default flaring channel the tip centre of gravity is then 4.7 mm
      downstream of variant 1's.

    On cylindrical walls both variants give the same blade.
  - The centre of gravity is the physical one, the area centroid on the stream surface. The conformal map
    does not carry it over from the profile's own frame except on cylinders, so the section is shifted until
    its centre of gravity is on the target. On the default hub the correction is 0.011 in m′ (about 0.4 mm).
  - *Take Z_ref, θ_ref from the current hub section* (also done automatically when you switch from the LE /
    TE lines) keeps the hub section where it is. It takes the hub's reference point before the re-stagger, so a
    re-staggered hub stays put too, and pressing it again changes nothing.
  - A reference point may lie anywhere on the stream surface, extensions included. A point outside the
    channel is reported as a design error.
- **Chord scale** per layer (Layers and Profiles → Placement): with reference-point stacking,
  c′ = scale × the chord of the prescribed edges. θ LE (lean) is not used there; θ_ref places the sections.
- The **Meridional Contour View** draws the stacking line (dash-dot) with the layers' reference points. It
  fills the blade between its actual LE and TE, which leave the LE / TE lines when the blade is stacked on a
  point or re-staggered.
- The information tables give each layer's stacking point (Z, θ) and chord scale.
- The re-stagger axis through the hub, mid-span or tip section follows the stacked section.
- Design files: a `"stacking"` block (`mode`, `variant`, `reference`, `z_ref`, `theta_ref_deg`) and each
  layer's `chord_scale`. Files without them load with the prescribed edges, as before.
- Python:
  - `stack_z`, `stack_r`, `stack_theta_deg` and `stack_mprime` for each layer;
  - `blade_le`, `blade_te`, `stacking_xr` and `stacking_xyz` for the mesh rows;
  - `DesignEvaluator.hub_reference_point(design_json, reference)`: the (Z_ref, θ_ref) of the button above.

## v1.0.3 — Blade Sections view (M-R·Theta, Z-R·Theta, axial X-Y) and centres of gravity

### New
- **Blade Sections (all layers)**: a new view content (*New Auxiliary View* or *Auxiliary View Content*). It
  draws every layer's section (Hub, 25 %, Mid, 75 %, Tip) in one plot: suction side blue, pressure side red,
  LE / TE closures green and the camber line thin dark. Green **LE and TE lines** join the leading and the
  trailing edges of all sections. Right-click → *Blade Sections Coordinates*:
  - **M, RTheta Plot**: m from each layer's LE line against r θ (mm);
  - **Z, RTheta Plot**: axial position Z against r θ (mm);
  - **X, Y Plot (axial view)**: the sections seen along the machine axis from downstream,
    X = r sin θ and Y = r cos θ (mm), with the hub and casing circles at the LE and TE radii dashed.
- The window title shows the chosen coordinates.
- **Centres of gravity**:
  - **Sections**: each layer's area centroid in its stream surface's (m, r θ) plane, mapped back to the
    surface. It is exact on cylindrical layers (the unrolled cylinder).
  - **Blade**: the volume centroid of the closed blade (the mesh rows with flat hub and tip caps, by the
    divergence theorem; uniform density).
  - Shown in the Blade Information tables (centre of gravity Z, R, θ for every layer), as a CG marker in the
    Blade-to-Blade View, and in the Blade Sections view as markers joined by a dashed line, with the blade's
    centre of gravity and volume in the status box.
  - In the 3D View: a dashed line through the section centres of gravity and a cross at the blade's.
  - Python: `section_area`, `cg_xyz`, `cg_z`, `cg_r`, `cg_theta_deg` and `cg_mprime` for each layer;
    `blade_volume` and `blade_cg` for the design.
- The information tables' area is the section's area on its stream surface; it was the chord-frame area
  scaled by the chord squared.
- The axis triads of the 3D View and the 3D Solid view are labelled in the designer's machine convention:
  - **Z axial**, as in the meridional input;
  - **X = r sin θ** and **Y = r cos θ**, as in the axial view;
  - X red, Y green, Z blue.

  They were labelled with the internal frame, where x is axial. Only the labels change (a cyclic relabelling,
  still right-handed); the geometry, the STEP / OBJ exports and the Python arrays keep the internal
  (axial, r sin θ, r cos θ) order.

### Fixed
- Re-staggered blades: the blade-to-blade view drew the camber line and the TE marker from the LE line rather
  than from the turned LE point. A blade without a re-stagger was not affected.
- The blade-to-blade status box gives the stagger actually used, including the re-stagger, and notes the
  re-stagger.

### Tests
- Centres of gravity:
  - polygon area and centroid against closed forms (rectangle, triangle, either orientation, far from the
    origin);
  - the closed-surface volume and centroid of an N-gon prism and a sheared prism in closed form;
  - on cylindrical layers, the section centre of gravity equals the chord-frame centroid mapped by the
    placement (10⁻¹²), and its area is the chord-frame area × (c′ r)²;
  - a re-stagger about the radial axis through the mid layer's centre of gravity leaves that centre of gravity
    in place;
  - the blade's centre of gravity is within 0.05 mm of the exact OCCT solid's centre of mass (measured
    0.016 mm), and its volume within 1 % (+0.47 %, from the flat caps);
  - pytest: the new fields.
- The designer self-test renders the Blade Sections view in all three coordinate systems (the content loop
  also renders it for every layer) and the tip section with a 15° re-stagger. That is 113 view renders.

## v1.0.2 — BladeGen Designer: meridional channel, five layers, conformal 3-D blade

BladeGen gets a native **C++ / Qt designer** with the layout of a professional blade-design program: MDI views
whose content is picked from a right-click menu, live plots in thin lines, and a meshed 3-D blade. The 3-D
part follows the DLR BladeGen parameterization (Voß et al., Computers and Fluids 292 (2025) 106588):
profiles in the (m′, θ) plane, stacked on the stream surfaces of a meridional channel.

### New

**Initial meridional configuration** (`src/blade/FlowPath`)
- Eight cylindrical coordinates of an axial compressor: hub and casing (Z, R) at the LE and the TE.
- Per wall, **Bell & Mehta** or **Linear**:
  - Bell & Mehta (NASA CR-182747, 1988): the fifth-order contour R = R_LE + ΔR (6ξ⁵ − 15ξ⁴ + 10ξ³) between
    the LE and TE, with zero slope and curvature at both ends, joining cylinders upstream and downstream (G2);
  - Linear: the straight line through the LE and TE points, continued straight.
- Inlet and outlet extensions as fractions of each wall's axial length.
- The walls are **exact** degree-5 B-splines; layer streamlines are their exact blend, so every layer's LE and
  TE lie on the straight LE / TE lines.

**Five layers**: Hub, 25 %, Mid, 75 % and Tip (the 25 % and 75 % layers are new), each with its own
construction method (NACA 65, Custom, NACA 65 + Custom) and its placement (stagger, lean θ_LE, 2-D solidity).

**Conformal mapping to 3-D** (`StreamSurface`): m′ = ∫ dm / r measured from the layer LE.
- **Cylindrical wrapping** (closed form x = x_LE + r m′) on constant-radius layers (straight walls).
- **Mapping onto the surface of revolution** r(z) everywhere else.
- (m′, θ) is the main design coordinate system; the blade-to-blade view also shows (m, rθ).

**Metal angle and thickness** (`src/profile/Camber`): the camber line is the locus of the centres of the
inscribed circles, smoothed by a least-squares cubic B-spline pinned at the LE / TE points. It gives:
- the local metal angle β along the chord of every layer;
- θ along the camber;
- the thickness t/c along the chord;
- the **metal-angle contour in the meridional cut** (and a thickness contour) on the (Z, R) plane.

**Meshed 3-D geometry** (`src/design/BladeDesign`) instead of a rendered OCCT solid:
- structured blade rows between the layers, clustered toward hub and tip (two-sided tanh, Vinokur 1983);
- hub and casing grids over one pitch;
- exports: OBJ mesh, layer CSV, and the OCCT blade solid as STEP / IGES / BREP / STL.

**BladeGen Designer** (`apps/designer`, `bladegen_designer`, Qt 6 Widgets)
- MDI sub-windows. **Right-click a view** for:
  - its content: 3D View, Blade-to-Blade View, Profile Section, Meridional Contour, Meridional Metal Angle /
    Thickness Contour, Meridional Curvature, Blade-to-Blade Curvature (SS / PS), LE/TE Theta, LE/TE Beta,
    LE/TE Parameters, Airfoil Area, Blade Information Table (current / all layers), Blade Angles, Beta / Theta /
    Beta vs. Theta / Thickness of all layers;
  - its layer;
  - the horizontal axis of the angle / thickness plots;
  - the blade-to-blade coordinates;
  - the plot options: legend, points, status strip, swap / proportional / flipped axes, symmetric-log scale;
  - Number of Points, Save Graph Data, Copy View Image, Zoom Fit.
- The plots are drawn in C++ (QPainter), so there is no matplotlib:
  - thin antialiased lines: **suction side blue, pressure side red, LE / TE closures green**;
  - wheel zoom, drag pan, double-click fit;
  - the cursor position is shown bottom left.
- The 3-D view (no OpenGL needed):
  - shaded or wireframe blade (red), hub (green) and casing (black) grids, tip cap (magenta), rotation axis;
  - drag to rotate.
- Parameter panel:
  - *Initial Meridional Configuration*, with the Bell & Mehta / Linear buttons per wall;
  - *Layers and Profiles*: the layer buttons, the method, knobs by family and the NACA 65 + Custom *Fit*;
  - *Blade and Mesh*.
- The *Geometry checks* dock lists errors and warnings per layer.
- Geometry runs on a worker thread; the newest request replaces a pending one. The five layers are built in
  parallel, and unchanged layer sections come from a cache.
- Files: `bladegen-design-2` JSON. `bladegen-design-1` files (the former section window) open, with the 25 % and 75 % layers
  blended.
- Launchers: `tools\designer.cmd` (Windows) and `tools/designer.sh` (Linux / macOS) build the designer on
  first use, then open it.
- `--screenshot out.png` and `--self-test dir`: CI renders every view content for every layer offscreen.

**Python**:
- `MeridionalConfig`, `WallPoints`, `WallLaw`, `FlowPath`, `StreamSurface`, `bell_mehta`;
- `DesignEvaluator.evaluate(design_json)`, which returns layers, mesh and fields as numpy arrays, and
  `DesignEvaluator.write_step(...)`;
- `design_defaults()`, `design_normalize()`.

**Build / CI**:
- Qt 6 in `environment.yml` (`qt6-main`) and in all CI jobs (apt `qt6-base-dev` on Linux);
- `BLADEGEN_BUILD_DESIGNER=AUTO|ON|OFF`: without Qt, everything else still builds;
- version 1.0.2.

### Designer refinements: draggable control points, thickness law, Rotor 67 default layers

**New**
- **Profile Section view: B-spline control points.** For Custom layers the control polygons of the suction
  side, the pressure side and both edge connectors are drawn (dotted), and the control points can be
  **dragged with the mouse**. Each drag changes the corresponding parameter, and the knobs follow on release:
  - SS: the LE and TE handles move along their fixed directions (`ss.handle_le_frac`, `ss.handle_te_frac`);
    CP2 moves freely (`ss.ctrl2_mprime`, `ss.ctrl2_theta`);
  - PS: the LE and TE handles (`ps.handle_*_frac`); a free PS point moves its m′ station (`ps.mprime_ctrl`),
    because its height is the thickness law;
  - edge points and connector control points follow from the edge parameters and are shown fixed.
- **Thickness plots show the thickness law** of a Custom layer with the layer's numbers:
  - t(ξ) = t_LE + (t_max − t_LE)(ξ/ξ_max)^a forward of ξ_max and t_TE + (t_max − t_TE)((1 − ξ)/(1 − ξ_max))^b
    aft of it;
  - the law is drawn as a dashed curve next to the measured inscribed-circle thickness.
- **Default design: NASA Rotor 67 custom layers.** The five layers start from the 14 TP-2879 sections fitted
  with the custom profile (`data/rotor67/rotor67_bladegen.json`), blended to 0 / 25 / 50 / 75 / 100 % span,
  with their stagger and LE angle. The profiles are chord-normalised, so they scale with the chosen
  meridional path (c′ cos γ = Δm′ of each layer; angles are kept). *File → New design (NACA 65 layers)*
  gives the v1.0.2 NACA 65 default.
- Python: `rotor67_layer(span)`, `design_naca65_defaults()`.

**Changed**
- Geometry-check markers are small x marks, black for warnings and red for errors (no more large orange
  triangles).

**Tests**
- Rotor 67 layers exact at the 14 sections, blended between them, with stagger rising toward the tip and
  invalid spans rejected.
- The default design builds without errors on every layer. Its control polygons match the parameters, and
  the law maximum equals t_max.
- Doubling the axial chord doubles the tip chord and keeps the angles. The NACA 65 default still evaluates.
- The designer self-test drags CP2 of the mid layer and checks that the design and the geometry follow.

### Designer update: OpenCASCADE solid view, re-stagger, equal-size layouts

**New**
- **3D Solid (OpenCASCADE)**, a view content and a toolbar button. It shows the blade with the hub and casing,
  drawn by OCCT's own viewer (V3d / AIS, OpenGL), inside the MDI window like the other views.
  - **Geometry** (right-click):
    - *Design mesh* (the default) updates with every parameter change, in 3 ms for one blade and 11 ms for
      24, before a 5–20 ms redraw. It is the design's own mesh rows as triangles with smooth normals and flat
      hub and tip caps, handed to OpenGL as primitive arrays. There is no OCCT modelling, Boolean or BRepMesh.
      It passes through every mesh point and is within 0.054 mm of the exact C2 skin through the same rows
      (1/800 of the chord), and its enclosed volume is within 0.5 % of the exact solid's.
    - *Exact OCCT solid* is the trimmed B-spline solid of the STEP export, built on a worker thread in about
      1.5–2 s per change.
  - **Blades** (right-click): a *single blade* on one pitch of hub and casing, centred on the blade, or the
    *whole annulus*, all n_blades copies at 2π / n_blades with the hub and casing over 360°. The copies share
    the blade's geometry and triangulation.
  - **Blade, Hub and Casing Display** (right-click), each set on its own:
    - *wireframe*: edges and iso lines;
    - *solid*: opaque, lit, with the face boundaries;
    - *shaded*: lit, smooth, see-through;
    - *hidden* (hub and casing only).
    The defaults are a solid blade and hub and a see-through casing.
  - **Standard View**: isometric, axial, meridional (from the side), radial (onto the blade-to-blade plane).
    Left drag rotates, middle drag (or shift + left) pans, the wheel zooms and a double click fits.
  - The status line in the view gives the update time.
  - The view needs the OCCT visualization toolkits (TKService, TKV3d, TKOpenGl, detected by CMake) and OpenGL
    2.1 with frame buffers. Without them it says why, and the mesh 3D View and the STEP export still work.
    For example, the software driver of a remote desktop or a CI machine (Windows "GDI Generic", OpenGL 1.1)
    gets the message instead of blank images. `BLADEGEN_NO_OPENGL=1` switches the view off.
- **Re-stagger the whole blade** (Initial Meridional Configuration → *Re-stagger the whole blade*): an angle
  Δγ turns every section about one radial axis.
  - + closes the blade (higher stagger: more stall margin), − opens it (more flow).
  - The axis is the radial line through the **hub** section (the rotor default), the **mid-span** or the
    **tip** section at a chord fraction x/c (0 = LE, 1 = TE), or through any **point** given by its axial
    position Z and angle θ.
  - Each section turns rigidly in its conformal (m′, θ) plane, by Δγ, about the point where the axis crosses
    its stream surface. So the sections stay on their stream surfaces and keep their chord, and the LE and TE
    leave the meridional LE / TE lines as a re-staggered blade does. The map is conformal, so **every metal
    angle on every layer changes by exactly Δγ**. On a cylindrical layer this is the rigid rotation of the
    unrolled section about the axis.
  - Saved in the design file as `"restagger": {"angle_deg", "axis": "hub" | "mid" | "tip" | "point",
    "chord_fraction", "z", "theta_deg"}`. Older files read as no re-stagger.
- **Window layouts** (toolbar and View menu): **6 views** (3 × 2), **4 views** (2 × 2) and **2 views** (side
  by side), all in equal cells that follow the window size.
  - Views beyond the count close; the active view is kept. Missing views open with the default contents.
  - The default layout is the 6-view grid.
- Python: each layer's dict has `mprime_le` and `theta_le_deg` (the placed LE point).
- **New default channel and blade count**:
  - casing LE (Z −11, R 106) and TE (10, 102);
  - hub LE (−15, 40) and TE (15, 43);
  - both walls Bell & Mehta, 18 blades;
  - the mesh is unchanged: 25 rows, 48 body and 16 edge points per side, tanh d = 1.5.

  The Rotor 67 layers scale with the new channel. Tests that need cylindrical walls set them explicitly.

**Changed**
- `LayerResult::mprime`, `m` and `mfrac` are measured from the placed LE point, and `mfrac` runs 0 → 1 between
  the placed LE and TE. Without a re-stagger nothing changes.
- `Occt::SolidOptions::blade_only` builds only the skin and the trimmed blade (no checks, properties, hub
  sector or fuse). The exact view uses it with the design's own mesh resolution, which cuts the build from
  about 9 s to about 1.5–2 s. The STEP export is unchanged.
- The exact view meshes its solid at 0.05 mm on the blade and 0.2 mm on the walls (it was 0.01 / 0.05 mm),
  and the display geometry is built once instead of twice. That cuts the time the window freezes from about
  2 s to about 0.5 s.
- `Design::build_scene_mesh`, `SceneMesh`, `TriangleMesh` and `blade_theta`: the fast display, in the core
  library (tested without OpenGL).

**Tests**
- Placement: the similarity map against its closed form, and the re-stagger as a rotation about the pivot
  (also mirrored).
- Re-stagger about the hub LE: the axis passes through the hub LE, and the hub LE stays. On every layer the
  re-staggered loop is the nominal loop rotated by Δγ about the axis crossing point (10⁻¹² in m′, θ), and every
  metal angle is +Δγ. The mesh rows at the layer spans equal the placed layers.
- The same checks for the mid-span, tip, hub-TE and arbitrary-point axes (the tip LE and the hub TE stay
  where they are).
- Invalid input: |Δγ| > 30°, a non-finite value, a pivot outside the chord, a total stagger above 80°, and an
  axis outside the channel (reported as a design error). Also JSON round trip, old files and unknown axis
  names.
- Scene:
  - a single blade on one pitch centred on it; the annulus has n_blades copies, each the blade located by
    k pitches about the axis (exact transform), and integrated volume and centroid check three of them;
  - neighbours do not touch;
  - hub and casing lie on their contours to 10⁻⁷ m;
  - 360° walls have n_blades times the area of one pitch (the area of a surface of revolution is the sweep
    × ∫ r ds);
  - invalid input.
- Scene mesh:
  - every design mesh point is a vertex, and the normals are unit length and point out of the blade;
  - the mesh is closed, and its enclosed volume (divergence theorem) is within 1 % of the exact solid's
    (measured +0.47 %);
  - the largest gap to the exact C2 skin at the cell centres is below 0.1 mm (measured 0.054 mm);
  - the walls lie on their contours, span exactly one pitch centred on the blade, and their normals face the
    passage;
  - the annulus and invalid input.
- The designer self-test:
  - re-staggers by 3° and checks every layer's β_LE;
  - builds the design mesh and the exact solid; where OpenGL is available, it renders the single blade, the
    annulus, a wireframe look and the exact solid;
  - checks that the 6, 4 and 2 view layouts give that many equal cells.
- pytest: re-stagger through the Python bindings.

### Tests
92 ctest entries in the `python` preset.

New C++ tests:
- the Bell & Mehta polynomial and its derivatives;
- the wall B-splines against the closed form, with G2 at their joins;
- linear walls;
- layer LE / TE on the LE / TE lines;
- the cylindrical wrap in closed form;
- m′ of a cone in closed form;
- the Hermite map against the adaptive integral (< 10⁻⁶ mm);
- conformality in 3-D (10⁻⁵ °);
- camber of concentric arcs, of a symmetric section and of a NACA 65 section;
- natural-spline weights and tanh rows;
- five layers placed and meshed onto hub and casing;
- the metal angle re-measured on the 3-D surface;
- cylindrical layers;
- JSON round trip and the bladegen-design-1 import;
- invalid input;
- the section cache;
- NACA 65 + Custom inside a design.

Also new:
- pytest for the new bindings;
- the designer self-test and a bad-file test.

### Known limitations
- The metal angle is the slope of the inscribed-circle camber line. On NACA 65 (A10) layers it is steep at the
  LE and TE, because the A10 slope is logarithmic there. The information table also lists the equivalent
  circular-arc angles.
- Stacking is at the LE point on the LE line (prescribed radial edges, Voß et al. Sec. 4.2, straight lines).
  Stacking at the centre of gravity, bow and sweep are not yet inputs. Lean is θ_LE per layer.
- The mesh rows interpolate the layer profiles pointwise in span. The STEP export skins those rows, so it
  matches the mesh but not a parametric spanwise law.
- The STEP export takes a few seconds and blocks the window.
- The OCCT solid view needs a real OpenGL 2.1 context. CI machines (offscreen, software GL) build the solid
  but cannot render it, so the view's rendering is tested only where a GPU driver is present.
- The re-stagger turns each section in its conformal plane, about the point where the radial axis crosses its
  stream surface. On a cylinder this is the rigid rotation of the unrolled section about the axis. On a
  conical stream surface it differs slightly from the rigid 3-D rotation of the blade, which would move
  sections off their stream surfaces.
- DCA and C4 are still listed as planned.

### Next steps
- DCA and C4 profiles; stacking options (centre of gravity, bow, sweep, curved LE / TE lines).
- Tip clearance and the hub fillet in the designer's solid export.
- The ENIGMAgrid structured grid from the designer's stream surfaces; TURBOlenza export of all five layers.

### NACA 65 series (from TURBOlenza)

The first integration step of the BladeGen → ENIGMAgrid → TURBOlenza chain: TURBOlenza's blade geometry, the C2
NACA 65 cascade section, now lives in BladeGen (it was prepared as v1.0.1, which was never released and is
withdrawn: its Tk / matplotlib section window is replaced by BladeGen Designer).

#### New

**Parametric NACA 65 series** (`src/naca/Naca65`, migrated from TURBOlenza `b2b/geometry.py`)
- NACA 65-(Cl0 A10) t/c for Cl0 0–2.5 and t/c 0.02–0.20, both variants: **Kovach–Sandercock** (taper from
  station 0.60 to an R = 0.008 c TE circle) and **Emery** (sharp TE).
- Reviewed before porting. **Kept:**
  - Table 4-1 with the 2.5 % misprint corrected;
  - the closed-form A10 mean line and the arc-length thickness station;
  - the C4 quintic thickness in √station;
  - the nose circle;
  - the C2 **monotone-curvature leading edge** (Whewell-form fair, the curvature falls from 1/R_LE without
    overshoot or S-bend);
  - the equivalent-arc angles.

  **Dropped:**
  - the PCHIP thickness (curvature sawtooth);
  - the legacy LE circle with its polar fair;
  - kernel smoothing and the least-squares C2 refit.

  **Replaced by BladeGen's own constructions:**
  - a C2 blend at the Kovach start;
  - a G2 minimum-jerk TE connector through the TE circle's point with its tangent and curvature;
  - a C2 closing piece to the sharp Emery TE.
- Agreement with TURBOlenza's frozen C2 geometry on five reference sections:

  | variant | body | LE | TE | difference at the TE |
  |---|---|---|---|---|
  | Emery | ≤ 5·10⁻⁸ c | ≤ 4·10⁻⁸ c | ≤ 1.1·10⁻⁵ c | |
  | Kovach | ≤ 3·10⁻⁶ c | ≤ 5·10⁻⁸ c | ≤ 1.9·10⁻⁴ c | the different, G2 TE connector |
- About 0.25 s per section including the quality check (the same checks as custom profiles; the sharp Emery
  TE is reported as information).
- `math/Constrained.hpp`: Lawson–Hanson NNLS (with the anti-cycling safeguard; without it the LE fair took 7 s),
  least-distance programming, LSEI, Gauss–Legendre of any order.

**NACA 65 + Custom** (`fit/NacaFit`, `fit_naca65`)
- The custom BladeGen parameters are inverse-fitted to a NACA 65 section, starting from the 65-(12)10 fit scaled
  to Cl0 and t. The result is an ordinary `ProfileParams`, so every custom knob then edits the NACA geometry.
- Accuracy and time: RMS ≈ 0.011 % c for Kovach sections in 3–5 s.

**Other**
- CLI: `bladegen naca65 <cl0> <t> [--variant …] [--csv sides.csv] [--fit custom.json]`; `bladegen version` now
  reports the CMake project version.
- Python: `Naca65Params`, `Naca65Section` (metadata, samples, camber line), `check_quality(section)`,
  `naca65_start`, `fit_naca65`.
- Tests:
  - Catch2 (NACA 65 and the solvers): Table 4-1, A10, the thickness law, the LE fair (monotone, G2), both TEs,
    the quality check over a Cl0 × t × variant grid, the TURBOlenza reference, NACA 65 + Custom, the solvers,
    invalid input.
  - pytest: the bindings.
  - CLI tests.

#### Known limitations (NACA 65)
- **Kovach station 0.60:** the table/taper corner (0.08°) is blended over ±0.005 in station. This keeps the
  surface within 3·10⁻⁶ c of TURBOlenza's, but the curvature has a short excursion there (TURBOlenza's frozen
  C2 surface has the same, κ −0.8…1.9 on the suction side), reported as suction-side inflection warnings. A
  wider blend removes it but moves away from TURBOlenza's geometry.
- **Pressure-side waviness:** the curvature of the thin NACA 65 pressure side changes sign a few times. This is
  inherited from the tabulated thickness and reported as a warning, not an error.
- **Emery in NACA 65 + Custom:** a sharp Emery TE cannot be matched by the round custom TE (RMS ≈ 0.07 % c,
  locally up to 0.36 % c at the TE).
- **DCA and C4** profiles are not implemented yet (listed as planned in the designer).
- All v1.0.0 limitations still apply.

#### Next steps (NACA 65)
- DCA and C4 sections.
- An option for a curvature-fair Kovach start (wider blend) and a sharp / cusped TE in the custom profile, so
  NACA 65 + Custom can match the Emery TE.
- Export of stations in ENIGMAgrid's geometry format, and a TURBOlenza run from the designer (solver option).

## v1.0.0 — first release

The C++23 / OpenCASCADE successor of the `blade2d` Python blade generator: a parametric C2 blade engine
with curvature-quality checks, inverse design from published coordinates, C2 NURBS solids, a rolling-ball
hub fillet, Python bindings and interactive notebooks. Tested on Windows, Linux and macOS.

### What is in it

**2-D profile engine** (`src/profile`)
- Port of the blade2d construction: DLR-style characteristic points at LE and TE, cubic B-spline suction
  side, pressure side coupled to it through a thickness law. The Python-compatible mode (`G2Legacy`)
  reproduces blade2d to 2·10⁻⁶ chord (regression test against the Python engine).
- New default LE / TE connector `Fair`: curvature-continuous (G2) at both junctions and the fairest such
  curve (minimum jerk energy, one KKT solve). Also `C2Parametric` (one parametrically C2 loop).
- 26 scalar parameters in six families (leading edge, trailing edge, suction side, thickness distribution,
  pressure side, connector) with generic access, spanwise blending and JSON I/O; blade2d
  `*_blade_params.json` files load directly.

**Geometry quality check** (`check_quality`) — the curvature plots the Python engine lacked
- Curvature of both sides from the LE point to the TE point; S-shaped / looping edges, curvature dips and
  double peaks, inflections, curvature oscillation, self-intersection, negative thickness, orientation,
  junction continuity. Problems are reported as errors / warnings with their location. The geometry is
  always built.

**3-D blades** (`src/blade`, `src/occt`)
- Meridional channel, exact arc length m and conformal m′ = ∫ dm / r (closed-form straight extrapolation
  past the ends); conformal placement with stagger, lean and mirror; natural-spline interpolation of every
  parameter in span.
- C2 NURBS skin (Piegl–Tiller skinning) that is **C2 across the periodic seam**: OCCT's periodic
  interpolation has a double knot there (C1). BladeGen solves the periodic system itself.
- Blade solid trimmed exactly on the hub and on the casing offset by the tip clearance; STEP AP214, IGES,
  STL, BREP export; tessellation for viewing.
- **Rolling-ball hub fillet** with an explicit B-rep (`BRepFilletAPI` fails on this junction). Ball error
  ~10⁻¹⁶ m, tangency 0.2–0.4°.

**Inverse design** (`src/fit`)
- Side-aware residuals, bounded Levenberg–Marquardt (Eigen), optional fairness penalty on edge defects.
- NACA 65-(12)10: RMS 0.0095 % chord in about 5 s.
- **NASA Rotor 67**: all 14 published sections (TP-2879) rebuilt on their stream surfaces with 3-D RMS
  1.6–16 µm and no geometry errors; C2 blade, 1.016 mm tip clearance, 1.78 mm hub fillet, STEP. The fitted
  definition is in `data/rotor67/rotor67_bladegen.json`.

**Tools**
- `bladegen` CLI (profile check, build, fit, demo) and `rotor67` (the whole Rotor 67 reconstruction).
- Python module `bladegen` (pybind11) and four notebooks: profile designer with knobs grouped by family and
  live curvature plots, 3-D blade with an **Add hub fillet** button, inverse fit, Rotor 67.
- `tools/jupyter.cmd` / `tools/jupyter.sh`: one command to build the module and open the notebooks.

**Quality assurance**
- 56 tests (Catch2, pytest, notebooks executed headless, CLI), each against a closed form, a published value
  or its defining equation, plus invalid input.
- CI on Windows (MSVC, OCCT 7.9), Linux (GCC, OCCT 7.6 from apt and an ASan/UBSan job) and macOS (AppleClang,
  OCCT 7.9). Zero warnings with warnings-as-errors everywhere.

### Findings on the Python engine and on OCCT
- The blade2d NACA 65-(12)10 fit has an S-shaped TE connector a few µm long. It is caused by fitting a round
  TE to the cusped NACA 65 trailing edge and is visible only in the curvature plot. The default
  `naca65_1210()` is a fairness-weighted refit with no geometry errors; the blade2d fit remains available as
  `naca65_1210_blade2d()`.
- blade2d's `enforce_c1` makes the thickness slopes at the peak equal in magnitude only; the law keeps a
  corner there (kept for compatibility, documented).
- OCCT periodic interpolation is C1 at the seam (this affects the STEP files of the Python Rotor 67
  scripts); `BRepFilletAPI` fails on the blade/hub edge.

### Known limitations
- **Fillet:** constant radius, hub only (no tip or casing fillet, no variable radius). The hub contour's
  curvature radius must exceed the fillet radius. Tangency to the walls is 0.2–0.4°, not exact, because the
  arcs are interpolated. About 5–10 s per blade.
- **Profile model:** the edges are always rounded, so cusped or sharp trailing edges (NACA 65) can only be
  approximated by a small round edge. The pressure side is shaped through three thickness stations. The
  power-law thickness distribution has a corner at its peak (it only places PS control points, so the PS
  itself stays smooth).
- **Continuity:** the profile is G2 (curvature-continuous) at the four junctions in the default mode;
  parametric C2 of the whole loop requires `C2Parametric`. The 3-D skin interpolates sampled section points,
  so it reproduces the analytic sections to interpolation accuracy (sub-µm), not exactly.
- **Inverse fit:** a local optimiser with a finite-difference Jacobian, so the result depends on the start
  point. No similarity transform is fitted, the PS stations are fixed, and a section takes 5–60 s.
- **Quality check:** the thresholds are heuristic. Circular published LEs give "double peak" warnings
  because a polynomial cannot have constant curvature.
- **Scope:** single blade rows. No splitter or tandem blades, squealer tips or blade-to-blade passage / fluid
  domain export (the Python scripts had one). Spanwise stacking is set by the station placements; there are
  no lean or sweep law curves.
- **Platforms:** macOS is tested in CI on Apple Silicon only. The notebooks need JupyterLab with ipywidgets.
- **Dependencies:** CppAD is not used (the fit's finite differences are fast enough at 24 parameters).

### Next steps
- **ENIGMAgrid coupling:** export section loops in (m′, θ), stream lines and the periodic passage boundaries
  directly in ENIGMAgrid's geometry format, so the mesher's blade-to-blade topologies read BladeGen output.
- **TURBOlenza as a solver option:** design → mesh (ENIGMAgrid) → flow solution (TURBOlenza) from one
  parameter set, then aerodynamic optimisation. Exact gradients (CppAD or an adjoint) become worthwhile at
  that point.
- Geometry: sharp / cusped TE and elliptic LE options, lean / sweep / bow stacking laws, variable-radius and
  tip fillets, passage and fluid-domain export with named faces, splitter blades.
- Inverse design: global start search, fitting the PS stations and a similarity transform, a full 3-D
  surface fit.

### Relation to ENIGMAgrid and TURBOlenza
BladeGen is the geometry stage of a three-part toolchain:

- **BladeGen** produces the blade: sections, NURBS skin, hub and casing, STEP / JSON.
- **ENIGMAgrid** (C++20 structured multi-block mesher) grids the blade-to-blade passage.
- **TURBOlenza** (flow solver, to be added as a solver option to these modules) computes the flow.

All three share the conformal (m′, θ) frame and the hub / casing contours, so a blade designed or
inverse-designed here can be meshed and analysed without a CAD round trip.
