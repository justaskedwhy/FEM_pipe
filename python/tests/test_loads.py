import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp


def test_point_load_produces_downward_tip_motion():
    _, result, _ = solve_inp(cantilever_deck("CPS4", 8, 2))
    disp = np.asarray(result.displacement)
    assert disp[:, 1].min() < -1e-4, "tip must deflect downward under -y load"


def test_point_load_on_inactive_2d_z_dof_is_ignored():
    deck = cantilever_deck("CPS4", 4, 2)
    deck = deck.replace("\n*STEP", "\n*CLOAD\n1, 3, 500.0\n*STEP")
    model, result, prepared = solve_inp(deck)
    assert prepared.point_load_dofs.size == 3, "z-dof load must be dropped, y tip loads kept"
    assert result.residual_norm < 1e-8
    assert np.asarray(result.displacement).max() > 0.0