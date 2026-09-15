"""Doc 06 analytical benchmarks, run via `python tests/run_benchmarks.py`.

Requires the compiled ``femcore`` extension. Exits non-zero if a benchmark
violates its acceptance criteria.

Acceptance criteria come from the docs themselves:
  - Cantilever (doc 06 2.1, 1.3): monotone convergence from below toward the
    analytical tip deflection as the mesh is refined (shear locking diminishes
    with refinement);
  - Plate-with-hole (doc 06 2.2 + finite-width Howland correction): FEM K_t in
    the physical band [2.0, 3.4];
  - Cube 3D (doc 06 2.3): sigma_zz = -70.0 MPa to machine precision.
"""

import math
import os
import sys

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "python", "tests"))

import numpy as np  # noqa: E402

try:
    import femcore  # noqa: E402
except ImportError:
    raise SystemExit(
        "femcore extension is not built. Run `pip install . -v` before running benchmarks."
    )
from conftest import cantilever_deck, solve_inp  # noqa: E402
from fem import meta  # noqa: E402
from fem.cli import _fill_solver_input  # noqa: E402
from fem.model_builder import build_model  # noqa: E402
from fem.parser import load_inp  # noqa: E402
from fem.preprocessor import preprocess  # noqa: E402

INPUTS = os.path.join(os.path.dirname(__file__), "inputs")

DELTA_TIP_ANALYTICAL = 100.0 * 1.0**3 / (3.0 * 210.0e9 * 8.333333e-7)  # 1.9048e-4 m


def _solve(inp_name):
    meta.configure(lambda name: femcore.element_meta(name))
    model = build_model(load_inp(os.path.join(INPUTS, inp_name)))
    solver_input = _fill_solver_input(femcore, preprocess(model))
    return model, femcore.fem_solve(solver_input)


def _tip_deflection(model, result):
    disp = np.asarray(result.displacement).reshape(model.n_nodes, 3)
    at_tip = [i for i in range(model.n_nodes) if abs(model.coords[i, 0] - 1.0) < 1e-12]
    return max(abs(disp[i, 1]) for i in at_tip)


def benchmark_cantilever():
    meta.configure(lambda name: femcore.element_meta(name))
    tips = []
    for nx, ny in ((10, 2), (20, 2), (40, 4), (80, 8)):
        model, result, _ = solve_inp(cantilever_deck("CPS4", nx, ny))
        assert result.converged, f"cantilever {nx}x{ny} did not converge"
        tip = _tip_deflection(model, result)
        rel = 100.0 * (tip - DELTA_TIP_ANALYTICAL) / DELTA_TIP_ANALYTICAL
        print(f"  cantilever {nx:>2}x{ny} tip = {tip:.6e} m  (rel {rel:+.2f}%)")
        tips.append(tip)

    # Doc 06 2.1 / 1.3: monotone convergence from below under refinement.
    for t in tips:
        assert 0.0 < t < DELTA_TIP_ANALYTICAL, "all meshes must under-predict tip (shear locking)"
    for a, b in zip(tips, tips[1:]):
        assert b >= a - 1e-9, "tip deflection must converge monotonically as h -> 0"


def benchmark_plate_hole():
    stress = np.asarray(_solve("plate_hole.inp")[1].elem_stress)
    sigma_max = stress[:, 0].max()
    kt = sigma_max / 1.0e6
    print(f"  plate-hole sigma_max = {sigma_max/1.0e6:.3f} MPa  (Kt = {kt:.3f})")
    assert 2.0 < kt < 3.4, "plate-hole Kt outside physical band [2.0, 3.4]"


def benchmark_cube3d():
    sigma = np.asarray(_solve("cube3d.inp")[1].elem_stress)
    assert np.isclose(sigma[0, 2], -70.0e6, rtol=1e-6), "sigma_zz must be -70 MPa"
    print(f"  cube3d sigma_zz = {sigma[0, 2]/1.0e6:.3f} MPa  (target -70.000 MPa)")


def main():
    benchmarks = [benchmark_cantilever, benchmark_plate_hole, benchmark_cube3d]
    failed = 0
    for fn in benchmarks:
        try:
            fn()
        except AssertionError as exc:
            failed += 1
            print(f"  FAIL [{fn.__name__}]: {exc}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())