import numpy as np
import pytest

from conftest import cantilever_deck, solve_inp


def test_postprocess_voigt_ordering_and_invariants():
    model, result, _ = solve_inp(cantilever_deck("CPS4", 10, 2))
    stress = np.asarray(result.elem_stress)
    strain = np.asarray(result.elem_strain)
    principal = np.asarray(result.elem_principal)
    vm = np.asarray(result.elem_von_mises)

    assert stress.shape == (20, 6)
    assert strain.shape == (20, 6)
    assert principal.shape == (20, 3)

    sxx, syy, szz, sxy = (stress[:, i] for i in range(4))
    assert np.allclose(szz, 0.0, atol=1e-9), "plane stress sigma_zz must be zero"

    for i in range(2):
        assert np.all(principal[:, i] >= principal[:, i + 1] - 1e-9)

    expected_vm = np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * sxy**2
    )
    assert np.allclose(vm, expected_vm, rtol=1e-9)


def test_plane_strain_has_nonzero_sigma_zz():
    deck = (
        "*NODE\n"
        "1, 0, 0\n2, 1, 0\n3, 1, 1\n4, 0, 1\n"
        "*ELEMENT, TYPE=CPE4, ELSET=P\n1, 1, 2, 3, 4\n"
        "*MATERIAL, NAME=M\n*ELASTIC\n210.0e9, 0.3\n"
        "*SOLID SECTION, ELSET=P, MATERIAL=M\n0.01\n"
        "*BOUNDARY\n"
        "1, 1, 2, 0.0\n2, 1, 1, 1e-3\n2, 2, 2, 0.0\n"
        "3, 1, 1, 1e-3\n3, 2, 2, -3e-4\n4, 1, 2, 0.0\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model, result, _ = solve_inp(deck)
    stress = np.asarray(result.elem_stress)
    strain = np.asarray(result.elem_strain)
    assert np.allclose(
        strain[:, 2], 0.0, atol=1e-12
    ), "plane strain eps_zz must be zero"
    assert np.all(stress[:, 2] > 1e6)


def test_3d_cube_uniaxial_compression():
    deck = (
        "*NODE\n"
        "1, 0,0,0\n2, 1,0,0\n3, 1,1,0\n4, 0,1,0\n"
        "5, 0,0,1\n6, 1,0,1\n7, 1,1,1\n8, 0,1,1\n"
        "*ELEMENT, TYPE=C3D8, ELSET=C\n1, 1,2,3,4,5,6,7,8\n"
        "*MATERIAL, NAME=Al\n*ELASTIC\n70.0e9, 0.33\n"
        "*SOLID SECTION, ELSET=C, MATERIAL=Al\n"
        "*BOUNDARY\n"
        "1, 1, 1, 0.0\n4, 1, 1, 0.0\n5, 1, 1, 0.0\n8, 1, 1, 0.0\n"
        "1, 2, 2, 0.0\n2, 2, 2, 0.0\n5, 2, 2, 0.0\n6, 2, 2, 0.0\n"
        "1, 3, 3, 0.0\n2, 3, 3, 0.0\n3, 3, 3, 0.0\n4, 3, 3, 0.0\n"
        "5, 3, 3, -1e-3\n6, 3, 3, -1e-3\n7, 3, 3, -1e-3\n8, 3, 3, -1e-3\n"
        "*STEP\n*STATIC\n*END STEP\n"
    )
    model, result, _ = solve_inp(deck)
    assert result.converged
    stress = np.asarray(result.elem_stress)
    assert np.isclose(stress[0, 2], -70.0e6, rtol=1e-6)
    assert np.allclose(stress[0, 1], stress[0, 0], rtol=1e-6, atol=1e-6)
    assert np.allclose(stress[0, 0], 0.0, atol=1e-6)
