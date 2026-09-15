import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp

DELTA_TIP_BEAM = 100.0 * 1.0**3 / (3.0 * 210e9 * 8.333333e-7)  # 1.9048e-4 m


def _tip_displacement(result, model):
    n = model.n_nodes
    disp = np.asarray(result.displacement).reshape(n, 3)
    best = 0.0
    for i in range(n):
        x = model.coords[i, 0]
        if abs(x - 1.0) < 1e-12:
            d = abs(disp[i, 1])
            if d > best:
                best = d
    return best


def test_swap_q4_to_t3_runs_both_and_converges():
    q4 = solve_inp(cantilever_deck("CPS4", 10, 2))
    t3 = solve_inp(cantilever_deck("CPS3", 10, 2))
    assert q4[1].converged
    assert t3[1].converged

    dq4 = _tip_displacement(q4[1], q4[0])
    dt3 = _tip_displacement(t3[1], t3[0])

    assert 0.0 < dt3 <= dq4, "same mesh: Q4 should be softer than T3"
    assert (
        dq4 <= DELTA_TIP_BEAM
    ), "Q4 must approach analytical tip deflection from below"


def test_convergence_monotone_under_refinement():
    prev = None
    for nx in (10, 20, 40, 80):
        model, result, _ = solve_inp(cantilever_deck("CPS4", nx, 2))
        tip = _tip_displacement(result, model)
        assert tip > 0.0
        if prev is not None:
            assert tip >= prev - 1e-9, "Q4 tip deflection must converge monotonically"
        prev = tip
