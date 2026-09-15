import numpy as np
import pytest

from fem.errors import ModelError
from fem.model_builder import build_model
from fem.parser import parse_inp


def test_build_model_minimal(minimal_deck):
    model = build_model(parse_inp(minimal_deck))
    assert model.dimension == 2
    assert model.n_nodes == 4
    assert model.n_elem == 1
    assert model.n_mat == 1
    assert model.coords.shape == (4, 3)
    assert np.allclose(model.coords[:, 2], 0.0)
    assert model.elem_conn == [[0, 1, 2, 3]]
    assert model.mat_thickness[0] == 0.01


def test_missing_node_raises():
    deck = (
        "*NODE\n1, 0.0, 0.0\n2, 1.0, 0.0\n3, 1.0, 1.0\n"
        "*ELEMENT, TYPE=CPS3\n100, 1, 2, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n1.0\n"
    )
    with pytest.raises(ModelError, match="does not exist"):
        build_model(parse_inp(deck))


def test_unregistered_type_raises():
    deck = (
        "*NODE\n1, 0, 0\n"
        "*ELEMENT, TYPE=CPS6\n1, 1, 1, 1, 1, 1, 1\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n1.0\n"
    )
    with pytest.raises(ModelError, match="not registered"):
        build_model(parse_inp(deck))


def test_wrong_node_count_raises():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPS3\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n1.0\n"
    )
    with pytest.raises(ModelError, match="has 4 nodes, expected 3"):
        build_model(parse_inp(deck))


def test_mixed_dimension_raises():
    deck = (
        "*NODE\n1, 0, 0, 0\n"
        "*ELEMENT, TYPE=CPS3\n1, 1, 1, 1\n"
        "*ELEMENT, TYPE=C3D4\n2, 1, 1, 1, 1\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n1.0\n"
    )
    with pytest.raises(ModelError, match="Mixed dimension"):
        build_model(parse_inp(deck))


def test_plane_mode_derived_from_element_family(minimal_deck):
    model = build_model(parse_inp(minimal_deck.replace("CPS4", "CPE4")))
    assert model.mat_plane_mode[0] == 1


def test_unassigned_element_raises():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n"
        "*ELEMENT, TYPE=CPS3, ELSET=A\n1, 1, 2, 3\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=OTHER, MATERIAL=M\n1.0\n"
    )
    with pytest.raises(ModelError, match="no assigned material"):
        build_model(parse_inp(deck))


def test_boundary_node_set_expansion():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*NSET, NSET=FIXED\n1, 4\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*BOUNDARY\nFIXED, 1, 2, 0.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model = build_model(parse_inp(deck))
    expanded = [b for b in model.boundaries if b.node_set is None]
    assert expanded == model.boundaries
    assert sorted(b.node_id for b in model.boundaries) == [1, 4]
    assert all(b.dof_first == 1 and b.dof_last == 2 for b in model.boundaries)


def test_cload_node_set_expansion():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*NSET, NSET=LOADED\n2, 3\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*CLOAD\nLOADED, 2, -50.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model = build_model(parse_inp(deck))
    assert sorted(p.node_id for p in model.point_loads) == [2, 3]
    assert all(p.magnitude == -50.0 for p in model.point_loads)


def test_boundary_unknown_node_set_raises(minimal_deck):
    deck = minimal_deck.replace("1, 1, 2, 0.0", "GHOST, 1, 2, 0.0")
    with pytest.raises(ModelError, match="undefined or empty node set"):
        build_model(parse_inp(deck))


def test_surface_pressure_expansion():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n5, 0.5, 0\n6, 1, 0.5\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n"
        "1, 1, 5, 6, 4\n2, 5, 2, 3, 6\n"
        "*ELSET, ELSET=EDGE\n1\n2\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*SURFACE, TYPE=ELEMENT, NAME=TOPP\nEDGE, S2\n"
        "*DSLOAD\nTOPP, P, -3.5\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model = build_model(parse_inp(deck))
    assert sorted(p.elem_id for p in model.pressure_loads) == [1, 2]
    assert all(p.face_label == "P2" for p in model.pressure_loads)
    assert all(p.magnitude == -3.5 for p in model.pressure_loads)


def test_surface_undefined_surface_raises():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*DSLOAD\nGHOST, P, -3.5\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    with pytest.raises(ModelError, match="undefined surface"):
        build_model(parse_inp(deck))


def test_surface_undefined_elset_raises():
    deck = (
        "*NODE\n1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*SURFACE, TYPE=ELEMENT, NAME=P1\nNOELEMS, S1\n"
        "*DSLOAD\nP1, P, -3.5\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    with pytest.raises(ModelError, match="undefined or empty element set"):
        build_model(parse_inp(deck))
