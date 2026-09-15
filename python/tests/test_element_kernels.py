import numpy as np
import pytest

femcore = pytest.importorskip("femcore")


@pytest.fixture(autouse=True)
def live_registry():
    from fem import meta

    meta.ELEMENT_META.clear()
    meta.configure(lambda name: femcore.element_meta(name))
    yield
    meta.configure(None)
    meta.ELEMENT_META.clear()


def test_registered_elements():
    assert sorted(femcore.registered_elements()) == [
        "C3D4", "C3D8", "CPE3", "CPE4", "CPS3", "CPS4"
    ]


def test_metadata_matches_vtk_table():
    assert femcore.element_meta("CPS3")["vtk"] == 5
    assert femcore.element_meta("CPS4")["vtk"] == 9
    assert femcore.element_meta("C3D4")["vtk"] == 10
    assert femcore.element_meta("C3D8")["vtk"] == 12


@pytest.mark.parametrize("etype,coords,thick,plane", [
    ("CPS3", np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0.0]], dtype=float), 0.01, 0),
    ("CPS4", np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float), 0.01, 0),
    ("C3D4", np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float), 1.0, 0),
    ("C3D8", np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float), 1.0, 0),
])
def test_stiffness_symmetric_psd_and_rigid_body(etype, coords, thick, plane):
    Ke = femcore.test_element_matrix(etype, coords, 210.0e9, 0.3, thick, plane)
    assert np.allclose(Ke, Ke.T, rtol=1e-10, atol=1e-8)

    evals = np.linalg.eigvalsh(Ke)
    scale = max(1.0, np.abs(evals).max())
    assert evals.min() >= -1e-12 * scale, "stiffness must be positive semi-definite"

    npe = coords.shape[0]
    dim = femcore.element_meta(etype)["dim"]
    for axis in range(dim):
        block = np.zeros(dim)
        block[axis] = 1.0
        rig = np.tile(block, npe)
        scale = np.abs(Ke).max()
        assert np.max(np.abs(Ke @ rig)) < 100.0 * np.finfo(float).eps * scale, \
            "rigid body translation must produce zero force"


def test_inverted_triangle_does_not_throw_large_area():
    inverted = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0.0]], dtype=float)
    Ke = femcore.test_element_matrix("CPS3", inverted, 210.0e9, 0.3, 0.01, 0)
    assert np.allclose(Ke, Ke.T, rtol=1e-10, atol=1e-8)


def test_degenerate_triangle_throws():
    degenerate = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0.0]], dtype=float)
    with pytest.raises(RuntimeError):
        femcore.test_element_matrix("CPS3", degenerate, 210.0e9, 0.3, 0.01, 0)


Q4_SQUARE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
H8_CUBE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
T4_TET = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)


def test_q4_pressure_force_bottom_edge_is_positive_y():
    fe = femcore.test_element_pressure("CPS4", Q4_SQUARE, 1, 100.0, 1.0)
    nodal = fe.reshape(-1, 2)
    assert np.allclose(nodal[0], [0.0, 50.0])
    assert np.allclose(nodal[1], [0.0, 50.0])
    assert np.allclose(nodal[2], 0.0)
    assert np.allclose(nodal[3], 0.0)
    assert np.isclose(nodal.sum(axis=0)[1], 100.0)


def test_q4_pressure_force_top_edge_is_negative_y():
    fe = femcore.test_element_pressure("CPS4", Q4_SQUARE, 3, 100.0, 1.0)
    nodal = fe.reshape(-1, 2)
    assert np.allclose(nodal[2], [0.0, -50.0])
    assert np.allclose(nodal[3], [0.0, -50.0])


def test_q4_pressure_force_scales_with_thickness():
    fe = femcore.test_element_pressure("CPS4", Q4_SQUARE, 1, 100.0, 0.01)
    fe_t1 = femcore.test_element_pressure("CPS4", Q4_SQUARE, 1, 100.0, 1.0)
    assert np.allclose(fe, fe_t1 * 0.01)


def test_h8_pressure_force_top_face_compresses():
    fe = femcore.test_element_pressure("C3D8", H8_CUBE, 2, 100.0)
    total = fe.reshape(-1, 3).sum(axis=0)
    assert np.allclose(total, [0.0, 0.0, -100.0])
    nodal = fe.reshape(-1, 3)
    assert np.allclose(nodal[4], [0.0, 0.0, -25.0])
    assert np.allclose(nodal[7], [0.0, 0.0, -25.0])


def test_t4_pressure_force_face_split_and_direction():
    fe = femcore.test_element_pressure("C3D4", T4_TET, 4, 100.0)
    nodal = fe.reshape(-1, 3)
    assert np.allclose(nodal[1], np.zeros(3))
    for i in (0, 2, 3):
        assert np.allclose(nodal[i], [100.0 / 6.0, 0.0, 0.0])
    assert np.allclose(nodal.sum(axis=0), [50.0, 0.0, 0.0])


def test_q4_b_at_is_nonzero_correct_shape():
    B = femcore.test_element_b_at("CPS4", Q4_SQUARE, 0.0, 0.0)
    assert B.shape == (3, 8)
    assert np.abs(B).max() > 1e-12