"""Module M2: FEM model builder.

Transforms the unvalidated RawModel from M1 into a validated, 0-indexed
FEMModel. Raises ModelError on semantic violations.
"""

import logging
from dataclasses import dataclass, field

import numpy as np

from .errors import ModelError
from .meta import element_meta
from .parser import RawBoundary, RawModel, RawPointLoad

logger = logging.getLogger("fem.model")

_PREFIX_PLANE_MODE = {"CPS": 0, "CPE": 1, "C3D": 2}


@dataclass
class FEMModel:
    coords: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    elem_type_name: list[str] = field(default_factory=list)
    elem_conn: list[list[int]] = field(default_factory=list)
    elem_mat: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    mat_E: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    mat_nu: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    mat_thickness: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.float64)
    )
    mat_plane_mode: np.ndarray = field(
        default_factory=lambda: np.zeros(0, dtype=np.int32)
    )
    boundaries: list = field(default_factory=list)
    point_loads: list = field(default_factory=list)
    pressure_loads: list = field(default_factory=list)
    node_ids: list[int] = field(default_factory=list)
    elem_ids: list[int] = field(default_factory=list)
    n_nodes: int = 0
    n_elem: int = 0
    n_mat: int = 0
    dimension: int = 0


def _expand_set_boundaries(boundaries, nsets):
    expanded = []
    for b in boundaries:
        if b.node_set is None:
            expanded.append(b)
            continue
        node_ids = nsets.get(b.node_set)
        if not node_ids:
            raise ModelError(
                f"*BOUNDARY references undefined or empty node set '{b.node_set}'"
            )
        for nid in node_ids:
            expanded.append(
                RawBoundary(
                    node_id=nid,
                    dof_first=b.dof_first,
                    dof_last=b.dof_last,
                    value=b.value,
                )
            )
    return expanded


def _expand_set_point_loads(point_loads, nsets):
    expanded = []
    for p in point_loads:
        if p.node_set is None:
            expanded.append(p)
            continue
        node_ids = nsets.get(p.node_set)
        if not node_ids:
            raise ModelError(
                f"*CLOAD references undefined or empty node set '{p.node_set}'"
            )
        for nid in node_ids:
            expanded.append(RawPointLoad(node_id=nid, dof=p.dof, magnitude=p.magnitude))
    return expanded


def _prefix_plane_mode(type_name: str) -> int:
    for prefix, mode in _PREFIX_PLANE_MODE.items():
        if type_name.startswith(prefix):
            return mode
    raise ModelError(f"Element type '{type_name}' has no known analysis mode prefix")


def build_model(raw: RawModel) -> FEMModel:
    if not raw.elements:
        raise ModelError("Model contains no elements")

    element_metas = {}
    for elem in raw.elements.values():
        if elem.type_name not in element_metas:
            try:
                element_metas[elem.type_name] = element_meta(elem.type_name)
            except ModelError:
                raise ModelError(f"Element type '{elem.type_name}' is not registered")

    dimensions = {m["dim"] for m in element_metas.values()}
    if len(dimensions) > 1:
        raise ModelError(
            f"Mixed dimension: elements span dimensions {sorted(dimensions)}"
        )
    dimension = dimensions.pop()

    node_order = list(raw.nodes.values())
    node_id_to_idx = {n.id: i for i, n in enumerate(node_order)}

    coords = np.zeros((len(node_order), 3), dtype=np.float64)
    for i, n in enumerate(node_order):
        coords[i, 0] = n.x
        coords[i, 1] = n.y
        coords[i, 2] = n.z

    elem_order = list(raw.elements.values())
    elem_conn: list[list[int]] = []
    elem_type_name: list[str] = []
    for elem in elem_order:
        meta = element_metas[elem.type_name]
        if len(elem.node_ids) != meta["npe"]:
            raise ModelError(
                f"Element {elem.id} of type '{elem.type_name}' has "
                f"{len(elem.node_ids)} nodes, expected {meta['npe']}"
            )
        conn = []
        for nid in elem.node_ids:
            if nid not in node_id_to_idx:
                raise ModelError(
                    f"Node {nid} referenced by element {elem.id} does not exist"
                )
            conn.append(node_id_to_idx[nid])
        elem_conn.append(conn)
        elem_type_name.append(elem.type_name)

    if not raw.materials:
        raise ModelError("Model contains no materials")
    for name, mat in raw.materials.items():
        if mat.E <= 0.0:
            raise ModelError(f"Material '{name}' must have E > 0")
        if not (-1.0 < mat.nu < 0.5):
            raise ModelError(f"Material '{name}' must satisfy -1 < nu < 0.5")

    raw_mat_names = [name for name in raw.materials]
    mat_name_to_idx = {n: i for i, n in enumerate(raw_mat_names)}
    n_mat = len(raw_mat_names)

    elem_mat = np.zeros(len(elem_order), dtype=np.int32)
    elem_mat[:] = -1
    elem_id_to_idx = {e.id: i for i, e in enumerate(elem_order)}
    for section in raw.sections:
        if section.material_name not in mat_name_to_idx:
            raise ModelError(
                f"*SOLID SECTION references unknown material '{section.material_name}'"
            )
        if section.elset_name not in raw.elsets:
            logger.warning(
                "[model] *SOLID SECTION elset '%s' is empty or undefined",
                section.elset_name,
            )
            continue
        for eid in raw.elsets[section.elset_name]:
            if eid not in elem_id_to_idx:
                raise ModelError(
                    f"*SOLID SECTION references unknown element {eid}"
                    f" in elset '{section.elset_name}'"
                )
            elem_mat[elem_id_to_idx[eid]] = mat_name_to_idx[section.material_name]

    for ei, mat_idx in enumerate(elem_mat):
        if mat_idx < 0:
            raise ModelError(
                f"Element {elem_order[ei].id} has no assigned material (missing *SOLID SECTION)"
            )

    mat_E = np.array([raw.materials[n].E for n in raw_mat_names], dtype=np.float64)
    mat_nu = np.array([raw.materials[n].nu for n in raw_mat_names], dtype=np.float64)
    mat_thickness = np.full(n_mat, 1.0, dtype=np.float64)
    mat_plane_mode = np.full(n_mat, -1, dtype=np.int32)

    for section in raw.sections:
        m = mat_name_to_idx[section.material_name]
        mat_thickness[m] = section.thickness

    for ei, type_name in enumerate(elem_type_name):
        m = int(elem_mat[ei])
        mode = _prefix_plane_mode(type_name)
        if mat_plane_mode[m] < 0:
            mat_plane_mode[m] = mode
        elif mat_plane_mode[m] != mode:
            raise ModelError(
                f"Material '{raw.mat_names[m]}' is used by mixed 2D analysis modes"
            )

    return FEMModel(
        coords=coords,
        elem_type_name=elem_type_name,
        elem_conn=elem_conn,
        elem_mat=elem_mat,
        mat_E=mat_E,
        mat_nu=mat_nu,
        mat_thickness=mat_thickness,
        mat_plane_mode=mat_plane_mode,
        boundaries=_expand_set_boundaries(raw.boundaries, raw.nsets),
        point_loads=_expand_set_point_loads(raw.point_loads, raw.nsets),
        pressure_loads=list(raw.pressure_loads),
        node_ids=[n.id for n in node_order],
        elem_ids=[e.id for e in elem_order],
        n_nodes=len(node_order),
        n_elem=len(elem_order),
        n_mat=n_mat,
        dimension=dimension,
    )
