import os

from fem.model_builder import build_model
from fem.parser import parse_inp
from fem.vtk_writer import write_vtu_manual


def test_manual_writer_ascii_structure(minimal_deck, tmp_path):
    model = build_model(parse_inp(minimal_deck))

    class DummyOutput:
        displacement = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -1e-5, 0.0],
            [0.0, -1e-5, 0.0],
        ]
        elem_stress = [[2.0e8, 0.0, 0.0, 0.0, 0.0, 0.0]]
        elem_strain = [[1e-3, -3e-4, -3e-4, 0.0, 0.0, 0.0]]
        elem_von_mises = [2.0e8]
        elem_principal = [[2.0e8, 0.0, 0.0]]

    out = tmp_path / "solution.vtu"
    write_vtu_manual(str(out), model, DummyOutput())

    text = out.read_text()
    assert 'type="UnstructuredGrid"' in text
    assert 'NumberOfPoints="4" NumberOfCells="1"' in text
    assert '<DataArray type="Int32" Name="connectivity" format="ascii">' in text
    assert "\n          4\n" in text
    assert '<DataArray type="UInt8" Name="types" format="ascii">' in text
    assert "9\n" in text
    assert 'Name="strain"' in text
    assert 'Name="von_mises"' in text
    assert "1.00000000e-03" in text
    assert "2.00000000e+08" in text


def test_manual_writer_no_hardcoded_vtk_code(minimal_deck, tmp_path):
    model = build_model(parse_inp(minimal_deck))

    class DummyOutput:
        displacement = [[0.0] * 3] * 4
        elem_stress = [[0.0] * 6]
        elem_strain = [[0.0] * 6]
        elem_von_mises = [0.0]
        elem_principal = [[0.0] * 3]

    out = tmp_path / "solution.vtu"
    write_vtu_manual(str(out), model, DummyOutput())
    assert "9\n" in out.read_text()
