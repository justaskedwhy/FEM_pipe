"""Command-line entry point (Layer 5) and pipeline orchestration.

Usage: femsolver --in <deck.inp> --out <solution.vtu> [--verbose]
Exit codes follow doc 01 section 7.1.
"""

import argparse
import logging
import sys

import numpy as np

from . import meta
from .errors import FEMError, InputError, ModelError
from .model_builder import build_model
from .parser import load_inp
from .preprocessor import preprocess
from .vtk_writer import write_vtu_manual

logger = logging.getLogger("fem")


def _fill_solver_input(femcore, prepared: "PreprocessedInput"):
    si = femcore.SolverInput()
    si.coords = np.ascontiguousarray(prepared.coords, dtype=np.float64)
    si.elem_type_name = list(prepared.elem_type_name)
    si.elem_conn_flat = np.ascontiguousarray(prepared.elem_conn_flat, dtype=np.int32)
    si.elem_conn_offsets = np.ascontiguousarray(prepared.elem_conn_offsets, dtype=np.int32)
    si.elem_mat = np.ascontiguousarray(prepared.elem_mat, dtype=np.int32)
    si.mat_E = np.ascontiguousarray(prepared.mat_E, dtype=np.float64)
    si.mat_nu = np.ascontiguousarray(prepared.mat_nu, dtype=np.float64)
    si.mat_thickness = np.ascontiguousarray(prepared.mat_thickness, dtype=np.float64)
    si.mat_plane_mode = np.ascontiguousarray(prepared.mat_plane_mode, dtype=np.int32)
    si.dof_map = np.ascontiguousarray(prepared.dof_map, dtype=np.int32)
    si.n_dofs_total = int(prepared.n_dofs_total)
    si.n_dofs_free = int(prepared.n_dofs_free)
    si.prescribed_dofs = np.ascontiguousarray(prepared.prescribed_dofs, dtype=np.int32)
    si.prescribed_vals = np.ascontiguousarray(prepared.prescribed_vals, dtype=np.float64)
    si.point_load_dofs = np.ascontiguousarray(prepared.point_load_dofs, dtype=np.int32)
    si.point_load_vals = np.ascontiguousarray(prepared.point_load_vals, dtype=np.float64)
    si.pressure_elem = np.ascontiguousarray(prepared.pressure_elem, dtype=np.int32)
    si.pressure_face = np.ascontiguousarray(prepared.pressure_face, dtype=np.int32)
    si.pressure_val = np.ascontiguousarray(prepared.pressure_val, dtype=np.float64)
    si.n_nodes = int(prepared.n_nodes)
    si.n_elem = int(prepared.n_elem)
    si.n_mat = int(prepared.n_mat)
    si.dimension = int(prepared.dimension)
    return si


def solve(inp_path: str, out_path: str, verbose: bool = False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(message)s")

    raw = load_inp(inp_path)
    logger.info("[parser]  Read %d nodes, %d elements from %s",
                len(raw.nodes), len(raw.elements), inp_path)

    model = build_model(raw)
    logger.info("[model]   Dimension = %d, active materials = %d",
                model.dimension, model.n_mat)

    prepared = preprocess(model)
    n_presc = prepared.prescribed_dofs.size
    logger.info("[preproc] DOFs total = %d, free = %d, prescribed = %d",
                prepared.n_dofs_total, prepared.n_dofs_free, n_presc)

    try:
        import femcore
    except ImportError:
        raise FEMError(
            "femcore extension is not built. Run `pip install . -v` first."
        )

    meta.configure(lambda name: femcore.element_meta(name))

    si = _fill_solver_input(femcore, prepared)

    result = femcore.fem_solve(si)

    if not result.converged:
        raise SolverFailure(
            "linear solve failed (LDLT and LU); check boundary conditions "
            "(rigid body modes must be constrained)"
        )

    logger.info("[solve]   SimplicialLDLT solve OK, residual norm = %.3e",
                float(result.residual_norm))

    disp = np.asarray(result.displacement)
    max_disp = float(np.max(np.linalg.norm(disp, axis=1))) if disp.size else 0.0
    logger.info("[post]    Max nodal displacement magnitude |u|_max = %.3e m", max_disp)

    write_vtu_manual(out_path, model, result)
    logger.info("[vtu]     Wrote %s", out_path)


class SolverFailure(FEMError):
    exit_code = 4


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="femsolver", description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="inp_path", required=True, help="input ABAQUS .inp deck")
    parser.add_argument("--out", dest="out_path", required=True, help="output VTK .vtu file")
    parser.add_argument("--verbose", action="store_true", help="debug logging")
    args = parser.parse_args(argv)

    try:
        solve(args.inp_path, args.out_path, args.verbose)
        return 0
    except FEMError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return exc.exit_code
    except RuntimeError as exc:
        print(f"femcore: {exc}", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"internal error: {exc}", file=sys.stderr)
        return 5


if __name__ == "__main__":
    sys.exit(main())