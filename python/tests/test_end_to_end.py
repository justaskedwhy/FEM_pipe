import os

import pytest

from conftest import cantilever_deck


def test_cli_end_to_end(tmp_path, monkeypatch):
    pytest.importorskip("femcore")
    import femcore  # ensure imported for meta

    inp = tmp_path / "cantilever.inp"
    inp.write_text(cantilever_deck("CPS4", 8, 2))

    out = tmp_path / "solution.vtu"
    from fem.cli import main

    code = main(["--in", str(inp), "--out", str(out)])
    assert code == 0
    assert out.exists()
    assert out.stat().st_size > 0


def test_cli_parser_error_exit_1(tmp_path, capsys):
    from fem.cli import main

    bad = tmp_path / "bad.inp"
    bad.write_text("*NODE\n1, oops, 0.0\n")
    code = main(["--in", str(bad), "--out", str(tmp_path / "x.vtu")])
    assert code == 1


def test_cli_model_error_exit_2(tmp_path):
    from fem.cli import main

    deck = tmp_path / "model.inp"
    deck.write_text(
        "*NODE\n1,0,0\n"
        "*ELEMENT, TYPE=CPS3\n1, 1, 2, 3\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n1e9, 0.3\n"
        "*SOLID SECTION, ELSET=B, MATERIAL=M\n1.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    code = main(["--in", str(deck), "--out", str(tmp_path / "x.vtu")])
    assert code == 2


_INPUTS = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "tests", "inputs"
)


def test_q4_pressure_deck_end_to_end(tmp_path):
    pytest.importorskip("femcore")
    from fem.cli import main

    out = tmp_path / "q4_pressure.vtu"
    code = main(["--in", os.path.join(_INPUTS, "q4_pressure.inp"), "--out", str(out)])
    assert code == 0
    assert out.exists() and out.stat().st_size > 0


def test_h8_pressure_deck_analytical_stress(tmp_path):
    femcore = pytest.importorskip("femcore")
    import numpy as np

    from fem import meta
    from fem.cli import _fill_solver_input
    from fem.model_builder import build_model
    from fem.parser import load_inp
    from fem.preprocessor import preprocess

    meta.configure(lambda name: femcore.element_meta(name))
    inp = os.path.join(_INPUTS, "h8_pressure.inp")
    model = build_model(load_inp(inp))
    prepared = preprocess(model)
    res = femcore.fem_solve(_fill_solver_input(femcore, prepared))
    assert res.converged
    sigma_zz = np.asarray(res.elem_stress)[0, 2]
    assert np.isclose(sigma_zz, -1.0e6, rtol=1e-6)
