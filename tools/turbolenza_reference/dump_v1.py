"""
dump_v1.py
==========
Reference data for the C++ port of the TURBOlenza v1 blade-to-blade solver (src/turbolenza).

Runs the v1 Python package (turbolenza.b2b, given by its repository path) on a few fixed cases
and writes every object of the method as JSON: the cascade and its C2 surface curves, the grid
(nodes, wall identity, wall-curve derivatives), the metrics, the coefficients A..E, the boundary
relations, the assembled sparse system F Psi = Q, psi, the flow, the static state, the Kutta
residual and correction, and the pass records. The C++ tests read these files and compare
module by module.

    python tools/turbolenza_reference/dump_v1.py <path-to-turbolenza-v1-repo> [out_dir]

out_dir defaults to tests/turbolenza/data. Every double is written with repr() (17 significant
digits, round-trip exact). Needs the v1 environment (NumPy, SciPy).
"""

from __future__ import annotations

import dataclasses
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np


def jsonable(obj, depth: int = 0):
    """Dataclasses, arrays, scalars, containers -> JSON values. Arrays: {"shape", "data"} with
    data flattened in C order (row major, last index fastest); string arrays as lists."""
    if depth > 12:
        return None
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if math.isfinite(v) else ("nan" if math.isnan(v) else ("inf" if v > 0 else "-inf"))
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind in "US":
            return {"shape": list(obj.shape), "strings": [str(s) for s in obj.ravel()]}
        if obj.dtype.kind == "b":
            return {"shape": list(obj.shape), "data": [bool(b) for b in obj.ravel()]}
        flat = np.asarray(obj, dtype=float).ravel()
        return {"shape": list(obj.shape), "data": [jsonable(float(v)) for v in flat]}
    if hasattr(obj, "tocoo") and hasattr(obj, "shape") and hasattr(obj, "nnz"):  # scipy sparse
        coo = obj.tocoo()
        return {"sparse_shape": list(coo.shape), "row": [int(r) for r in coo.row], "col": [int(c) for c in coo.col],
                "val": [float(v) for v in coo.data]}
    if isinstance(obj, slice):
        return {"start": obj.start, "stop": obj.stop}
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        out = {"__type__": type(obj).__name__}
        for f in dataclasses.fields(obj):
            out[f.name] = jsonable(getattr(obj, f.name), depth + 1)
        return out
    if isinstance(obj, dict) or hasattr(obj, "items"):
        try:
            return {str(k): jsonable(v, depth + 1) for k, v in obj.items()}
        except Exception:
            return None
    if isinstance(obj, (list, tuple)):
        return [jsonable(v, depth + 1) for v in obj]
    return None     # functions, scipy splines, ...: not data


def grid_json(grid):
    """The Grid with its properties (wall identity, wall-curve derivatives, regions)."""
    out = jsonable(grid)
    out["wall_identity"] = jsonable(grid.wall_identity)
    out["region"] = jsonable(np.asarray(grid.region))
    out["dm"] = grid.dm
    out["deta"] = grid.deta
    out["has_wall_curves"] = bool(grid.has_wall_curves)
    if grid.has_wall_curves:
        out["wall_curve_derivative"] = {str(k): {"lower": jsonable(np.asarray(lo)), "upper": jsonable(np.asarray(up))}
                                       for k in (0, 1, 2) for lo, up in [grid.wall_curve_derivative(k)]}
    return out


def cascade_json(cascade):
    out = jsonable(cascade)
    out["metadata"] = {k: jsonable(v) for k, v in cascade.metadata.items() if jsonable(v) is not None}
    return out


class _Trace:
    """Observer: keeps the pass records and the full fields of the first passes."""

    def __init__(self, keep: int) -> None:
        self.keep, self.records, self.fields = keep, [], []

    def start(self, info) -> None:
        self.info = info

    def update(self, record) -> None:
        self.records.append(record)

    def finish(self, solution, error) -> None:
        self.error = None if error is None else f"{type(error).__name__}: {error}"

    def inspect(self, fields) -> None:
        if len(self.fields) < self.keep:
            self.fields.append(fields)


def pass_fields_json(f, solver):
    from turbolenza.b2b.boundary_conditions import build_boundary_conditions
    from turbolenza.b2b.matrix_assembly import assemble_from_coefficients
    from turbolenza.b2b.linear_solver import solve_linear_system
    from turbolenza.b2b import kutta

    bc = build_boundary_conditions(f.grid)
    assembly = assemble_from_coefficients(f.coefficients, f.metrics, bc)
    linear = solve_linear_system(assembly)
    out = {
        "record": jsonable(f.record),
        "grid": grid_json(f.grid),
        "metrics": jsonable(f.metrics),
        "boundary_conditions": jsonable(bc),
        "density_resident": jsonable(np.asarray(f.density_resident)),
        "density_target": jsonable(np.asarray(f.density_target)),
        "density_next": jsonable(np.asarray(f.density_next)),
        "coefficients": jsonable(f.coefficients),
        "assembly": jsonable(assembly),
        "psi": jsonable(np.asarray(f.psi)),
        "psi_resolve": jsonable(np.asarray(linear.psi)),
        "flow": jsonable(f.flow),
        "static_state": jsonable(f.static_state),
        "mach": jsonable(np.asarray(f.mach)),
        "mass_flow_per_depth": f.mass_flow_per_depth,
        "c_ref": f.c_ref,
    }
    try:
        residual = kutta.kutta_residual(f.flow, f.grid)
        out["kutta_residual"] = jsonable(residual)
        out["exit_angle_correction"] = jsonable(kutta.exit_angle_correction(residual, f.flow, f.metrics))
    except Exception as error:                  # noqa: BLE001 - recorded, not hidden
        out["kutta_error"] = f"{type(error).__name__}: {error}"
    return out


def settings(**changes):
    from turbolenza.b2b.controller import CaseSettings
    base = dict(case_name="ref", solver_mode="coupled_kutta", fixed_exit_angle_deg=None, gamma=1.4,
                gas_constant=287.05, inlet_total_pressure=101325.0, inlet_total_temperature=288.15, omega=0.0,
                stream_surface_radius=None, inlet_absolute_tangential_velocity=None, blade_family="NACA65",
                design_lift_coefficient=1.2, thickness_ratio=0.10, meanline_type="A10",
                profile_variant="kovach_sandercock", surface_continuity="C2", n_camber_stations=101, chord=0.1,
                solidity=1.0, stagger_deg=33.0, inlet_mach=0.65, inlet_flow_angle_deg=45.0,
                mass_flow_per_depth_override=None, grid_level="very_coarse", n_passage=None, n_eta=None,
                upstream_extension_chords=1.0, downstream_extension_chords=1.0, density_tolerance=1e-8,
                velocity_tolerance=None, max_density_iterations=12, density_relaxation=0.5, stop_at_sonic=True,
                kutta_tolerance=1e-8, max_kutta_iterations=12, kutta_relaxation=1.0,
                initial_exit_angle_mode="blade_te_angle", initial_exit_angle_deg=None, measure_runtime=False)
    base.update(changes)
    return CaseSettings(**base)


def run_case(name: str, s, keep: int, out_dir: Path) -> None:
    from turbolenza.b2b import controller, solver as sv
    case = controller.build_case(s)
    trace = _Trace(keep)
    run = controller.run_case(case, observers=(trace,))
    sol = run.solution
    data = {
        "name": name,
        "settings": jsonable(s),
        "cascade": cascade_json(case.cascade),
        "config": jsonable(case.config),
        "first_grid": grid_json(case.grid),
        "records": [jsonable(r) for r in trace.records],
        "error": None if run.error is None else f"{type(run.error).__name__}: {run.error}",
        "passes": [pass_fields_json(f, sv) for f in trace.fields],
    }
    if sol is not None:
        from turbolenza.b2b import postprocess as pp
        data["solution"] = {
            "converged": sol.converged, "iterations": sol.iterations, "exit_angle": sol.exit_angle,
            "kutta_residual": sol.kutta_residual, "normalized_kutta_residual": sol.normalized_kutta_residual,
            "mass_flow_per_depth": sol.mass_flow_per_depth, "inlet": jsonable(sol.inlet),
            "density": jsonable(np.asarray(sol.density)), "psi": jsonable(np.asarray(sol.psi)),
            "mach": jsonable(np.asarray(sol.mach)), "flow": jsonable(sol.flow), "static_state": jsonable(sol.static_state),
            "grid": grid_json(sol.grid),
        }
        try:
            data["surfaces"] = jsonable(pp.blade_surfaces(sol))
            data["summary"] = jsonable(pp.summarize(sol))
            data["surface_pressure_coefficient"] = jsonable(pp.surface_pressure_coefficient(sol))
        except Exception as error:              # noqa: BLE001
            data["postprocess_error"] = f"{type(error).__name__}: {error}"
    path = out_dir / f"{name}.json"
    path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    print(f"{path.name}: {len(trace.records)} passes, error={data['error']}, {path.stat().st_size / 1e6:.2f} MB")


def grid_family(out_dir: Path) -> None:
    """The five standard levels and one explicit grid of the reference cascade (no solve)."""
    from turbolenza.b2b import grid as gr
    from turbolenza.b2b.config import GridSettings
    from turbolenza.b2b.geometry import CascadeSpec, make_cascade
    cascade = make_cascade(CascadeSpec(blade_type="NACA65", solidity=1.0, stagger_deg=33.0, cl0=1.2, t_max_over_c=0.10))
    beta1, kte = math.radians(45.0), cascade.kappa_te
    data = {"cascade": cascade_json(cascade), "inlet_angle": beta1, "exit_angle": kte, "levels": {}}
    data["transitions"] = jsonable(gr.passage_transitions(cascade))
    for level in gr.GRID_LEVELS:
        sg = gr.build_grid(cascade, level, inlet_angle=beta1, exit_angle=kte)
        g = sg.grid
        data["levels"][level] = grid_json(g)
        data["levels"][level]["info"] = jsonable(sg.info)
        data["levels"][level]["boundaries"] = jsonable(g.boundaries)
        data["levels"][level]["rebuilt_exit_20deg"] = jsonable(np.asarray(gr.rebuild_exit(g, math.radians(20.0)).nodes))
    explicit = gr.make_grid(cascade, GridSettings(n_passage=45, n_eta=21), beta1, kte)
    data["explicit_45x21"] = grid_json(explicit)
    data["explicit_45x21_boundaries"] = jsonable(explicit.boundaries)
    path = out_dir / "grid_family.json"
    path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    print(f"{path.name}: {path.stat().st_size / 1e6:.2f} MB")


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    v1 = Path(sys.argv[1]).resolve()
    sys.path.insert(0, str(v1 / "src"))
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parents[2] / "tests" / "turbolenza" / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    import subprocess
    commit = subprocess.run(["git", "-C", str(v1), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (out_dir / "SOURCE.txt").write_text(f"turbolenza v1 commit {commit}\nnumpy {np.__version__}\n", encoding="utf-8")
    grid_family(out_dir)
    # M1 0.65 at 45 deg, very_coarse: 12 coupled passes, the first two in full (pass-by-pass port check)
    run_case("m065_very_coarse_12", settings(), keep=2, out_dir=out_dir)
    # M1 0.30: converges on very_coarse (end-to-end check of the converged state)
    run_case("m030_very_coarse_converged", settings(inlet_mach=0.30, density_tolerance=1e-9, kutta_tolerance=1e-9,
                                                    max_density_iterations=400, max_kutta_iterations=400),
             keep=1, out_dir=out_dir)
    # fixed exit angle, explicit grid 45 x 21
    run_case("fixed_exit_explicit_45x21", settings(solver_mode="fixed_exit", grid_level=None, n_passage=45, n_eta=21,
                                                   inlet_mach=0.5, max_density_iterations=60, density_tolerance=1e-9),
             keep=1, out_dir=out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
