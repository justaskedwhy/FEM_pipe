import numpy as np
import pytest

from fem.errors import ModelError
from fem.model_builder import build_model
from fem.parser import parse_inp
from fem.preprocessor import preprocess


def test_dof_numbering_2d(minimal_deck):
    pre = preprocess(build_model(parse_inp(minimal_deck)))
    assert pre.n_dofs_total == 8
    assert pre.n_dofs_free == 6
    dof_map = pre.dof_map
    assert dof_map.shape == (4, 3)
    for n in range(4):
        assert dof_map[n, 0] == 2 * n
        assert dof_map[n, 1] == 2 * n + 1
        assert dof_map[n, 2] == -1


def test_dof_numbering_3d():
    deck = (
        "*NODE\n1, 0,0,0\n2, 1,0,0\n3, 1,1,0\n4, 0,1,0\n"
        "5, 0,0,1\n6, 1,0,1\n7, 1,1,1\n8, 0,1,1\n"
        "*ELEMENT, TYPE=C3D8, ELSET=B\n1, 1,2,3,4,5,6,7,8\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n70.0e9, 0.33\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n"
        "*BOUNDARY\n1, 3, 3, 0.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    pre = preprocess(build_model(parse_inp(deck)))
    assert pre.dimension == 3
    assert pre.n_dofs_total == 24
    assert pre.n_dofs_free == 23
    assert pre.dof_map[0, 2] == 2
    assert pre.dof_map[7, 2] == 23


def test_csr_offsets(minimal_deck):
    pre = preprocess(build_model(parse_inp(minimal_deck)))
    assert pre.elem_conn_offsets.tolist() == [0, 4]
    assert pre.elem_conn_flat.tolist() == [0, 1, 2, 3]


def test_2d_z_dof_loads_ignored(minimal_deck):
    deck = minimal_deck + "\n*CLOAD\n4, 3, 5.0\n"
    pre = preprocess(build_model(parse_inp(deck)))
    assert len(pre.point_load_dofs) == 1
    assert pre.point_load_vals[0] == -100.0


def test_pressure_conversion(minimal_deck):
    deck = minimal_deck + "\n*DLOAD\n1, P2, 50.0\n"
    pre = preprocess(build_model(parse_inp(deck)))
    assert pre.pressure_elem.tolist() == [0]
    assert pre.pressure_face.tolist() == [2]
    assert pre.pressure_val.tolist() == [50.0]


def test_pressure_face_out_of_range_raises(minimal_deck):
    deck = minimal_deck + "\n*DLOAD\n1, P5, 50.0\n"
    with pytest.raises(ModelError, match="out of range"):
        preprocess(build_model(parse_inp(deck)))


def test_boundary_unknown_node_raises(minimal_deck):
    deck = minimal_deck.replace("1, 1, 2, 0.0", "99, 1, 2, 0.0")
    with pytest.raises(ModelError, match="unknown node"):
        preprocess(build_model(parse_inp(deck)))
