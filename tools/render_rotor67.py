"""Render NASA Rotor 67 as rebuilt by BladeGen (assets/rotor67.png).

    python tools/render_rotor67.py            (bladegen module on PYTHONPATH: build preset 'python')

Uses the fitted definition data/rotor67/rotor67_bladegen.json (written by apps/rotor67), builds the
C2 blade with the 1.78 mm hub fillet on its hub sector, draws the 14 fitted stations on it and renders
the scene of notebook 04 with VTK, on a transparent background
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
    ren.SetBackgroundAlpha(0.0)  # transparent PNG: blends with light and dark page themes
    ren.UseFXAAOn()
    kit = vtk.vtkLightKit()  # key, fill, back and head lights: even shading like plotly
    kit.SetKeyLightIntensity(0.9)
    kit.AddLightsToRenderer(ren)
    for a in actors:
        ren.AddActor(a)
    win = vtk.vtkRenderWindow()
    win.SetOffScreenRendering(1)
    win.AddRenderer(ren)
    win.SetSize(*size)
    win.SetMultiSamples(0)
    win.SetAlphaBitPlanes(1)
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
    w2i.SetInputBufferTypeToRGBA()
    w2i.ReadFrontBufferOff()
    w2i.Update()
    wr = vtk.vtkPNGWriter()
    wr.SetFileName(str(path))
    wr.SetInputConnection(w2i.GetOutputPort())
    wr.Write()
    print("wrote", path)


def tube(points, rgb, radius):
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(points, dtype=float), deep=True))
    line = vtk.vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(len(points))
    for i in range(len(points)):
        line.GetPointIds().SetId(i, i)
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(line)
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetLines(cells)
    tf = vtk.vtkTubeFilter()
    tf.SetInputData(pd)
    tf.SetRadius(radius)
    tf.SetNumberOfSides(12)
    tf.Update()
    return actor(tf.GetOutput(), rgb, 0.1)


# the same scene as notebook 04: blade (with fillet) on its hub sector, plus the 14 fitted stations
hv, hf, _ = sh.hub.tessellate(2e-4, 0.25)
scene = [actor(polydata(v, f), (0.69, 0.77, 0.87)),          # plotly 'lightsteelblue'
         actor(polydata(hv * 1000.0, hf), (0.80, 0.80, 0.80), 0.1)]
for s in blade.sections(60, 200):
    xyz = s["xyz"] * 1000.0
    scene.append(tube(xyz, (0.10, 0.12, 0.18), 0.22))
render(scene, REPO / "assets" / "rotor67.png", (1600, 1300),
       {"direction": (-0.35, 1.0, 0.45), "up": (0, 0, 1), "zoom": 1.25})
