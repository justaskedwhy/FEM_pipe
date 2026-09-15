import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp


def test_assembly_residual_contract():
    for etype in ("CPS4", "CPS3"):
        _, result, _ = solve_inp(cantilever_deck(etype, 6, 2))
        if not result.converged:
            continue
        assert result.residual_norm < 1e-8
        assert np.asarray(result.reaction).shape == (result.displacement.shape[0], 3)


def test_rigid_body_insufficient_bcs_fails():
    deck = (
        "*NODE\n"
        "1, 0,0\n2, 1,0\n3, 1,1\n4, 0,1\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1,2,3,4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*CLOAD\n2, 1, 100.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    _, result, _ = solve_inp(deck)
    assert result.converged is False