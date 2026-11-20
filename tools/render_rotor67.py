"""Render NASA Rotor 67 as rebuilt by BladeGen (assets/rotor67.png, assets/rotor67_blade.png).

    python tools/render_rotor67.py            (bladegen module on PYTHONPATH: build preset 'python')

Uses the fitted definition data/rotor67/rotor67_bladegen.json (written by apps/rotor67), builds the
C2 blade with the 1.78 mm hub fillet, copies it 22 times around the axis and renders it with VTK
(off-screen; VTK ships with the conda-forge OCCT build).
"""
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

import bladegen as bg

REPO = Path(__file__).resolve().parents[1]
data = bg.read_rotor67(str(REPO / "data" / "rotor67" / "Rotor_14_Sections_and_Flowpath_EN.txt"))
blade = bg.Blade(bg.read_blade_json(str(REPO / "data" / "rotor67" / "rotor67_bladegen.json")))
loops = [s["xyz"] for s in blade.sections(60, 200)]
mid = blade.section(0.5, 60, 200)["xyz"]
inside = 0.5 * (mid[len(mid) // 4] + mid[3 * len(mid) // 4])
sh = bg.build_shapes(loops, data.hub_curve(), data.sections[-1].streamline, inside, data.tip_clearance, 1.78e-3, 1.0)
print(sh.fillet_message)
body = sh.blade_with_fillet if sh.fillet_ok else sh.blade
v, f, _ = body.tessellate(2e-4, 0.25)
v = v * 1000.0
tri = v[f]




def rot_x(p, a):  # rotation that increases theta (y = r sin, z = r cos)
    c, s = np.cos(a), np.sin(a)
    q = p.copy()
    q[..., 1] = c * p[..., 1] + s * p[..., 2]
    q[..., 2] = -s * p[..., 1] + c * p[..., 2]
    return q


def polydata(verts, faces):
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(verts, dtype=float), deep=True))
    cells = np.hstack([np.full((len(faces), 1), 3), faces]).astype(np.int64).ravel()
    ca = vtk.vtkCellArray()
    ca.SetCells(len(faces), numpy_to_vtkIdTypeArray(cells, deep=True))
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetPolys(ca)
    nm = vtk.vtkPolyDataNormals()
    nm.SetInputData(pd)
    nm.SetFeatureAngle(60.0)
    nm.Update()
    return nm.GetOutput()


def actor(pd, rgb, spec=0.25):
    m = vtk.vtkPolyDataMapper()
    m.SetInputData(pd)
    a = vtk.vtkActor()
    a.SetMapper(m)
    pr = a.GetProperty()
    pr.SetColor(*rgb)
    pr.SetSpecular(spec)
    pr.SetSpecularPower(30)
    pr.SetInterpolationToPhong()
    return a


def render(actors, path, size, cam):
    ren = vtk.vtkRenderer()
    ren.SetBackground(1, 1, 1)
    for a in actors:
        ren.AddActor(a)
    win = vtk.vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.AddRenderer(ren)
    win.SetSize(*size)
    win.SetMultiSamples(8)
    c = ren.GetActiveCamera()
    c.SetFocalPoint(0, 0, 0)
    c.SetPosition(*cam["direction"])  # view direction only; ResetCamera frames the scene
    c.SetViewUp(*cam["up"])
    ren.ResetCamera()
    c.Zoom(cam.get("zoom", 1.0))
    ren.ResetCameraClippingRange()
    win.Render()
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(win)
    w2i.Update()
    wr = vtk.vtkPNGWriter()
    wr.SetFileName(str(path))
    wr.SetInputConnection(w2i.GetOutputPort())
    wr.Write()
    print("wrote", path)


# hub: surface of revolution of the flowpath hub over the bladed region
hub = data.hub_curve()
xr = np.array([hub.at_m(m) for m in np.linspace(0, hub.length, 600)]) * 1000
z = np.linspace(-20, 105, 140)
rh = np.interp(z, xr[:, 0], xr[:, 1])
th = np.linspace(0, 2 * np.pi, 361)
Z, T = np.meshgrid(z, th, indexing="ij")
R = np.interp(Z, z, rh) - 0.2  # 0.2 mm inside: no z-fighting with the blade footprints
hv = np.stack([Z, R * np.sin(T), R * np.cos(T)], axis=-1).reshape(-1, 3)
nz, nt = Z.shape
idx = np.arange(nz * nt).reshape(nz, nt)
q = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], axis=-1).reshape(-1, 4)
hf = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])

blades = []
for k in range(data.n_blades):
    blades.append(actor(polydata(rot_x(v, 2 * np.pi * k / data.n_blades), f), (0.30, 0.50, 0.80)))
hub_actor = actor(polydata(hv, hf), (0.70, 0.70, 0.72), 0.1)
render([hub_actor] + blades, REPO / "assets" / "rotor67.png", (1600, 1200),
       {"direction": (-1.0, -0.85, 0.55), "up": (0, 0, 1), "zoom": 1.5})
# one blade with its fillet, close up near the root
render([actor(polydata(v, f), (0.30, 0.50, 0.80))], REPO / "assets" / "rotor67_blade.png", (1200, 1200),
       {"direction": (-1.0, 0.35, 0.3), "up": (0, 0, 1), "zoom": 1.1})
