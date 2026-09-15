import pytest

from fem import meta
from fem.errors import ModelError

FAKE_META = {
    "CPS3": {
        "vtk": 5,
        "npe": 3,
        "dim": 2,
        "dofs": 6,
        "strain_components": 3,
        "n_gauss": 1,
        "n_faces": 3,
        "faces": [0, 1, 2],
        "vtk_permutation": [0, 1, 2],
    },
    "CPE3": {
        "vtk": 5,
        "npe": 3,
        "dim": 2,
        "dofs": 6,
        "strain_components": 3,
        "n_gauss": 1,
        "n_faces": 3,
        "faces": [0, 1, 2],
        "vtk_permutation": [0, 1, 2],
    },
    "CPS4": {
        "vtk": 9,
        "npe": 4,
        "dim": 2,
        "dofs": 8,
        "strain_components": 3,
        "n_gauss": 4,
        "n_faces": 4,
        "faces": [0, 1, 2, 3],
        "vtk_permutation": [0, 1, 2, 3],
    },
    "CPE4": {
        "vtk": 9,
        "npe": 4,
        "dim": 2,
        "dofs": 8,
        "strain_components": 3,
        "n_gauss": 4,
        "n_faces": 4,
        "faces": [0, 1, 2, 3],
        "vtk_permutation": [0, 1, 2, 3],
    },
    "C3D4": {
        "vtk": 10,
        "npe": 4,
        "dim": 3,
        "dofs": 12,
        "strain_components": 6,
        "n_gauss": 1,
        "n_faces": 4,
        "faces": [0, 1, 2, 3],
        "vtk_permutation": [0, 1, 2, 3],
    },
    "C3D8": {
        "vtk": 12,
        "npe": 8,
        "dim": 3,
        "dofs": 24,
        "strain_components": 6,
        "n_gauss": 8,
        "n_faces": 6,
        "faces": list(range(8)),
        "vtk_permutation": list(range(8)),
    },
}

REGISTERED_NAMES = sorted(FAKE_META)


@pytest.fixture(autouse=True)
def fake_registry():
    meta.ELEMENT_META.clear()

    def provider(name):
        try:
            return dict(FAKE_META[name])
        except KeyError:
            raise ModelError(f"Element type '{name}' is not registered")

    meta.configure(provider)
    yield
    meta.configure(None)
    meta.ELEMENT_META.clear()


@pytest.fixture
def minimal_deck():
    return """\
*NODE
1, 0.0, 0.0
2, 1.0, 0.0
3, 1.0, 1.0
4, 0.0, 1.0

*ELEMENT, TYPE=CPS4, ELSET=PLATE
1, 1, 2, 3, 4

*MATERIAL, NAME=Steel
*ELASTIC
210.0e9, 0.3

*SOLID SECTION, ELSET=PLATE, MATERIAL=Steel
0.01

*BOUNDARY
1, 1, 2, 0.0

*CLOAD
4, 2, -100.0

*STEP
*STATIC
*END STEP
"""


def femcore_or_skip():
    femcore = pytest.importorskip("femcore")
    return femcore


def cantilever_deck(elem_type, nx, ny, load=-100.0):
    L, H, t = 1.0, 0.1, 0.01
    dx, dy = L / nx, H / ny
    lines = ["*NODE"]

    def nid(i, j):
        return j * (nx + 1) + i + 1

    for j in range(ny + 1):
        for i in range(nx + 1):
            lines.append(f"{nid(i,j)}, {i*dx:.8f}, {j*dy:.8f}")

    lines.append(f"*ELEMENT, TYPE={elem_type}, ELSET=BEAM")
    e = 0
    for j in range(ny):
        for i in range(nx):
            n00, n10 = nid(i, j), nid(i + 1, j)
            n11, n01 = nid(i + 1, j + 1), nid(i, j + 1)
            if elem_type in ("CPS4", "CPE4"):
                e += 1
                lines.append(f"{e}, {n00}, {n10}, {n11}, {n01}")
            else:
                e += 1
                lines.append(f"{e}, {n00}, {n10}, {n11}")
                e += 1
                lines.append(f"{e}, {n00}, {n11}, {n01}")

    lines += [
        "*MATERIAL, NAME=Steel",
        "*ELASTIC",
        "210.0e9, 0.3",
        "*SOLID SECTION, ELSET=BEAM, MATERIAL=Steel",
        f"{t}",
        "*BOUNDARY",
    ]
    for j in range(ny + 1):
        lines.append(f"{nid(0,j)}, 1, 2, 0.0")
    lines += ["*CLOAD"]
    for j in range(ny + 1):
        lines.append(f"{nid(nx,j)}, 2, {load/(ny+1):.8f}")
    lines += ["*STEP", "*STATIC", "*END STEP"]
    return "\n".join(lines)


def solve_inp(deck):
    femcore = femcore_or_skip()
    from fem.cli import _fill_solver_input
    from fem.model_builder import build_model
    from fem.parser import parse_inp
    from fem.preprocessor import preprocess

    model = build_model(parse_inp(deck))
    prepared = preprocess(model)
    solver_input = _fill_solver_input(femcore, prepared)
    result = femcore.fem_solve(solver_input)
    return model, result, prepared
