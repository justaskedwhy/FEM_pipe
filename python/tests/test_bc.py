import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp


def test_reactions_balance_point_loads():
    deck = (
        "*NODE\n"
        "1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*BOUNDARY\n1, 1, 2, 0.0\n4, 1, 2, 0.0\n"
        "*CLOAD\n2, 1, 100.0\n3, 1, 100.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model, result, _ = solve_inp(deck)
    assert result.converged
    reaction = np.asarray(result.reaction)
    rx_sum = float(reaction[:, 0].sum())
    assert np.isclose(rx_sum, -200.0, rtol=1e-6)
    assert np.isclose(float(reaction[:, 1].sum()), 0.0, atol=1e-8)
    assert result.residual_norm < 1e-8


def test_patch_with_pure_prescribed_displacements():
    deck = (
        "*NODE\n"
        "1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*BOUNDARY\n"
        "1, 1, 2, 0.0\n2, 1, 1, 1e-3\n2, 2, 2, 0.0\n"
        "3, 1, 1, 1e-3\n3, 2, 2, -3e-4\n4, 1, 2, 0.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model, result, _ = solve_inp(deck)
    assert result.converged
    assert result.residual_norm < 1e-8
    disp = np.asarray(result.displacement)
    assert np.allclose(disp[2], [1e-3, -3e-4, 0.0])