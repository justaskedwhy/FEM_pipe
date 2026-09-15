"""Module M3: Preprocessor / DOF manager.

Converts the FEMModel into flat, contiguous, typed SolverInput buffers.
Global DOF numbering follows doc 03 section 5.2:
    dof(n, d) = n * dim + d   (d in {0,1,2})
inactive directions are stored as -1 in dof_map.
"""

import logging
from dataclasses import dataclass, field

import numpy as np

from .errors import ModelError
from .meta import element_meta
from .model_builder import FEMModel

logger = logging.getLogger("fem.preproc")


@dataclass
class PreprocessedInput:
    coords: np.ndarray = field(
        default_factory=lambda: np.zeros((0, 3), dtype=np.float64)
    )
    elem_type_name: list[str] = field(default_factory=list)
    elem_conn_flat: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    elem_conn_offsets: np.ndarray = field(
        default_factory=lambda: np.zeros(1, dtype=np.int32)
    )
    elem_mat: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    mat_E: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    mat_nu: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    mat_thickness: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    mat_plane_mode: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    dof_map: np.ndarray = field(
        default_factory=lambda: np.full((0, 3), -1, dtype=np.int32)
    )
    n_dofs_total: int = 0
    n_dofs_free: int = 0
    prescribed_dofs: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    prescribed_vals: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    point_load_dofs: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    point_load_vals: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    pressure_elem: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    pressure_face: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    pressure_val: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    n_nodes: int = 0
    n_elem: int = 0
    n_mat: int = 0
    dimension: int = 0


def preprocess(model: FEMModel) -> PreprocessedInput:
    dim = model.dimension
    n_nodes = model.n_nodes

    dof_map = np.full((n_nodes, 3), -1, dtype=np.int32)
    counter = 0
    for n in range(n_nodes):
        for d in range(dim):
            dof_map[n, d] = counter
            counter += 1
    n_dofs_total = counter

    elem_conn_flat = (
        np.concatenate([np.asarray(c, dtype=np.int32) for c in model.elem_conn])
        if model.elem_conn
        else np.zeros(0, dtype=np.int32)
    )
    offsets = np.zeros(model.n_elem + 1, dtype=np.int32)
    run = 0
    for ei, conn in enumerate(model.elem_conn):
        run += len(conn)
        offsets[ei + 1] = run

    node_id_to_idx = {nid: i for i, nid in enumerate(model.node_ids)}
    elem_id_to_idx = {eid: i for i, eid in enumerate(model.elem_ids)}

    if len(node_id_to_idx) != model.n_nodes:
        raise ModelError("Duplicate or missing node ids in FEMModel")
    if len(elem_id_to_idx) != model.n_elem:
        raise ModelError("Duplicate or missing element ids in FEMModel")

    prescribed: dict[int, float] = {}
    for b in model.boundaries:
        if b.node_id not in node_id_to_idx:
            raise ModelError(f"*BOUNDARY references unknown node {b.node_id}")
        n = node_id_to_idx[b.node_id]
        for d in range(b.dof_first - 1, b.dof_last):
            if d >= dim:
                continue
            dof = int(dof_map[n, d])
            prescribed[dof] = b.value
    prescribed_dofs = np.asarray(list(prescribed), dtype=np.int32)
    prescribed_vals = np.asarray([prescribed[d] for d in prescribed], dtype=np.float64)
    n_dofs_free = n_dofs_total - len(prescribed)

    load_dofs: dict[int, float] = {}
    for p in model.point_loads:
        if p.node_id not in node_id_to_idx:
            raise ModelError(f"*CLOAD references unknown node {p.node_id}")
        n = node_id_to_idx[p.node_id]
        d = p.dof - 1
        if d >= dim:
            continue
        dof = int(dof_map[n, d])
        load_dofs.setdefault(dof, 0.0)
        load_dofs[dof] += p.magnitude
    point_load_dofs = np.asarray(list(load_dofs), dtype=np.int32)
    point_load_vals = np.asarray([load_dofs[d] for d in load_dofs], dtype=np.float64)

    pressure = []
    for p in model.pressure_loads:
        if p.elem_id not in elem_id_to_idx:
            raise ModelError(f"*DLOAD references unknown element {p.elem_id}")
        face_num = int(p.face_label[1:])
        elem_idx = elem_id_to_idx[p.elem_id]
        n_faces = element_meta(model.elem_type_name[elem_idx])["n_faces"]
        if not 1 <= face_num <= n_faces:
            raise ModelError(
                f"*DLOAD face {p.face_label} on element {p.elem_id} of type "
                f"'{model.elem_type_name[elem_idx]}' is out of range (1..{n_faces})"
            )
        pressure.append((elem_idx, face_num, p.magnitude))
    pressure_elem = np.asarray([t[0] for t in pressure], dtype=np.int32)
    pressure_face = np.asarray([t[1] for t in pressure], dtype=np.int32)
    pressure_val = np.asarray([t[2] for t in pressure], dtype=np.float64)

    return PreprocessedInput(
        coords=model.coords.astype(np.float64, copy=True),
        elem_type_name=list(model.elem_type_name),
        elem_conn_flat=elem_conn_flat,
        elem_conn_offsets=offsets,
        elem_mat=model.elem_mat.astype(np.int32, copy=True),
        mat_E=model.mat_E.astype(np.float64, copy=True),
        mat_nu=model.mat_nu.astype(np.float64, copy=True),
        mat_thickness=model.mat_thickness.astype(np.float64, copy=True),
        mat_plane_mode=model.mat_plane_mode.astype(np.int32, copy=True),
        dof_map=dof_map,
        n_dofs_total=n_dofs_total,
        n_dofs_free=n_dofs_free,
        prescribed_dofs=prescribed_dofs,
        prescribed_vals=prescribed_vals,
        point_load_dofs=point_load_dofs,
        point_load_vals=point_load_vals,
        pressure_elem=pressure_elem,
        pressure_face=pressure_face,
        pressure_val=pressure_val,
        n_nodes=model.n_nodes,
        n_elem=model.n_elem,
        n_mat=model.n_mat,
        dimension=model.dimension,
    )
