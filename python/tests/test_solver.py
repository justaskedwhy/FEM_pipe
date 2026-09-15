import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp


def test_single_element_q4_patch_exact_recovery():
    deck = (
        "*NODE\n"
        "1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n5, 0.5, 0\n"
        "6, 1, 0.5\n7, 0.5, 1\n8, 0, 0.5\n9, 0.5, 0.5\n"
        "*ELEMENT, TYPE=CPS4, ELSET=P\n"
        "1, 1, 5, 9, 8\n2, 5, 2, 6, 9\n3, 8, 9, 7, 4\n4, 9, 6, 3, 7\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*BOUNDARY\n"
    )
    boundary = []
    for n, (x, y) in enumerate(
        [(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0), (1, 0.5), (0.5, 1), (0, 0.5)],
        start=1,
    ):
        ux = 1e-3 * x
        uy = -3e-4 * y
        boundary.append(f"{n}, 1, 1, {ux}")
        boundary.append(f"{n}, 2, 2, {uy}")
    deck += "\n".join(boundary) + "\n*STEP\n*STATIC\n*END STEP\n"

    model, result, _ = solve_inp(deck)
    assert result.converged
    assert result.residual_norm < 1e-8

    disp = np.asarray(result.displacement)
    center = disp[8]
    assert abs(center[0] - 1e-3 * 0.5) < 1e-10
    assert abs(center[1] + 3e-4 * 0.5) < 1e-10

    strain = np.asarray(result.elem_strain)
    assert np.allclose(strain[:, 0], 1e-3, rtol=1e-9)
    assert np.allclose(strain[:, 1], -3e-4, atol=1e-12)
    assert result.elem_von_mises[0] > 0.0


def test_cantilever_q4_converged_and_reasonable():
    model, result, prepared = solve_inp(cantilever_deck("CPS4", 10, 2))
    assert result.converged
    assert result.residual_norm < 1e-8
    disp = np.asarray(result.displacement).reshape(model.n_nodes, 3)
    assert np.max(np.abs(disp[:, 1])) > 1e-4
