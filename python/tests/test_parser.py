import logging

import pytest

from fem.errors import InputError
from fem.parser import parse_inp


def test_parses_minimal_deck(minimal_deck):
    model = parse_inp(minimal_deck)
    assert len(model.nodes) == 4
    assert len(model.elements) == 1
    assert model.elements[1].type_name == "CPS4"
    assert model.elements[1].node_ids == [1, 2, 3, 4]
    assert model.materials["Steel"].E == 210.0e9
    assert model.materials["Steel"].nu == 0.3
    assert model.sections[0].thickness == 0.01
    assert model.boundaries[0].node_id == 1
    assert model.point_loads[0].node_id == 4
    assert model.point_loads[0].magnitude == -100.0


def test_comments_and_blanks_ignored():
    deck = "** a comment\n\n\n*NODE\n1, 0.0, 0.0\n"
    model = parse_inp(deck)
    assert 1 in model.nodes


def test_case_insensitive_keywords():
    deck = "*node\n1, 0.0, 0.0\n"
    model = parse_inp(deck)
    assert 1 in model.nodes


def test_unsupported_keyword_skipped_with_warning(caplog):
    deck = (
        "*NODE\n1, 0.0, 0.0\n"
        "*EL PRINT\n1\n"
        "*ELEMENT, TYPE=CPS4\n"
        "1, 1, 1, 1, 1\n"
    )
    with caplog.at_level(logging.WARNING):
        model = parse_inp(deck)
    assert len(model.elements) == 1


def test_malformed_float_raises_input_error():
    with pytest.raises(InputError):
        parse_inp("*NODE\n1, not-a-number, 0.0\n")


def test_duplicate_node_raises():
    with pytest.raises(InputError):
        parse_inp("*NODE\n1, 0.0, 0.0\n1, 1.0, 0.0\n")


def test_boundary_set_target_parsed(minimal_deck):
    model = parse_inp("*BOUNDARY\nNSETA, 1, 2, 0.0\n")
    assert len(model.boundaries) == 1
    assert model.boundaries[0].node_set == "NSETA"
    assert model.boundaries[0].node_id == 0
    assert model.boundaries[0].dof_first == 1
    assert model.boundaries[0].dof_last == 2


def test_cload_set_target_parsed(minimal_deck):
    model = parse_inp("*CLOAD\nLOADNODES, 2, -10.0\n")
    assert len(model.point_loads) == 1
    assert model.point_loads[0].node_set == "LOADNODES"
    assert model.point_loads[0].dof == 2
    assert model.point_loads[0].magnitude == -10.0


def test_standalone_nset_block(minimal_deck):
    deck = "*NSET, NSET=BOTTOM\n1, 2\n3\n" "*NSET, NSET=TOP\n4\n" "*NODE\n5, 0, 0\n"
    model = parse_inp(deck)
    assert model.nsets["BOTTOM"] == [1, 2, 3]
    assert model.nsets["TOP"] == [4]


def test_standalone_elset_block(minimal_deck):
    deck = "*ELSET, ELSET=PARTS\n1\n2, 3\n"
    model = parse_inp(deck)
    assert model.elsets["PARTS"] == [1, 2, 3]


def test_nset_requires_name():
    with pytest.raises(InputError, match="requires NSET="):
        parse_inp("*NSET\n1, 2\n")


def test_dload_face_label_parsed(minimal_deck):
    deck = minimal_deck + "\n*DLOAD\n1, P1, 100.0\n"
    model = parse_inp(deck)
    assert model.pressure_loads[0].face_label == "P1"
    assert model.pressure_loads[0].magnitude == 100.0


def test_generate_form_in_nset():
    deck = "*NSET, NSET=RANGE, GENERATE\n1, 10, 2\n"
    model = parse_inp(deck)
    assert model.nsets["RANGE"] == [1, 3, 5, 7, 9]


def test_generate_form_in_elset_default_step():
    deck = "*ELSET, ELSET=ALL, GENERATE\n1, 5\n"
    model = parse_inp(deck)
    assert model.elsets["ALL"] == [1, 2, 3, 4, 5]


def test_generate_form_invalid_step():
    with pytest.raises(InputError):
        parse_inp("*NSET, NSET=R, GENERATE\n1, 10, 0\n")


def test_elset_references_other_set():
    deck = (
        "*ELSET, ELSET=PART_A\n1, 2, 3\n"
        "*ELSET, ELSET=PART_B\n4\n"
        "*ELSET, ELSET=ALL\nPART_A, PART_B\n"
    )
    model = parse_inp(deck)
    assert model.elsets["ALL"] == [1, 2, 3, 4]


def test_nset_references_undefined_set_raises():
    with pytest.raises(InputError, match="references undefined set"):
        parse_inp("*NSET, NSET=GHOST\nMISSING\n")


def test_surface_and_dsload_parsed():
    deck = (
        "*ELSET, ELSET=EDGE\n1\n2\n"
        "*SURFACE, TYPE=ELEMENT, NAME=PRESS\nEDGE, S2\n"
        "*DSLOAD\nPRESS, P, -5.0\n"
    )
    model = parse_inp(deck)
    assert model.surfaces["PRESS"] == [("EDGE", "S2")]
    assert len(model.pressure_loads) == 1
    assert model.pressure_loads[0].surface_name == "PRESS"
    assert model.pressure_loads[0].magnitude == -5.0


def test_part_instance_scoping():
    deck = (
        "*PART, NAME=MAIN\n"
        "*NODE\n1, 0, 0\n2, 1, 0\n"
        "*ELEMENT, TYPE=CPS3, ELSET=MAIN\n1, 1, 2, 1\n"
        "*END PART\n"
        "*ASSEMBLY, NAME=A\n"
        "*INSTANCE, NAME=INST1, PART=MAIN\n"
        "*END INSTANCE\n"
        "*NSET, NSET=BC, INSTANCE=INST1\n1, 2\n"
        "*END ASSEMBLY\n"
        "*BOUNDARY\nBC, 1, 2, 0.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model = parse_inp(deck)
    assert "part:MAIN:MAIN" in model.elsets
    assert model.nsets["instance:INST1:BC"] == [1, 2]
    assert model.boundaries[0].node_set == "instance:INST1:BC"
