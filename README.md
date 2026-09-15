# fem-project

Hybrid Python–C++ linear elastic finite element solver. Reads ABAQUS `.inp`
decks, solves the static equilibrium system, and writes VTK `.vtu` results.

```
.inp ──► M1 parser ──► M2 model ➤ M3 preprocessor ──► femcore.fem_solve() ──► .vtu
                  Python orchestration          │ C++ core (Eigen)
                                                └─ single SolverInput in, one SolverOutput out
```

## Features

- **Four element families** (ABAQUS names): `CPS3`/`CPE3` (T3), `CPS4`/`CPE4`
  (Q4), `C3D4` (T4), `C3D8` (H8). Kernels self-register via the
  `REGISTER_ELEMENT` macro.
- **Static linear elasticity**: plane stress / plane strain / 3D, isotropic
  materials, nodal `*BOUNDARY` and `*CLOAD`, surface `*DLOAD` pressure.
- **Element sets**: `*NSET`/`*ELSET` blocks and string set targets in
  `*BOUNDARY`/`*CLOAD`.
- **ASCII VTK `.vtu` output** ready for ParaView (PointData displacement +
  reactions, CellData stress/strain/von Mises/principal).
- **Documentation-first**: the numbered specs in `docs/` are authoritative and
  govern implementation (doc 03 §1.1). The build is validated against
  analytical benchmarks and a 16-item consistency audit (doc 07).

## Requirements

| Tool / Library | Version | Purpose |
| :--- | :--- | :--- |
| Python | >= 3.10 | Orchestration & I/O |
| NumPy | >= 1.23 | Contiguous numeric buffers |
| CMake | >= 3.16 | C++ build |
| C++ compiler | C++17 | Numerical core (GCC 9+, Clang 10+, MSVC 2019+) |
| Eigen | >= 3.4 | Sparse/dense linear algebra (header-only) |
| pybind11 | >= 2.11 | Python <-> C++ memory bridge |
| meshio | >= 5.0 | optional `.vtu` exporter (fallback writer) |

## Install

Python orchestration depends only on NumPy; the numerical core is compiled at
install time.

```bash
pip install . -v
```

System prerequisites (Debian/Ubuntu):

```bash
sudo apt install -y build-essential cmake python3-dev libeigen3-dev \
                    pybind11-dev ninja-build
```

## Quick start

```bash
# Solve a validated deck and write a VTU
femsolver --in tests/inputs/cantilever.inp --out solution_cantilever.vtu [--verbose]

# Or through the Makefile
make run-input  IN=tests/inputs/cantilever.inp  OUT=solution_cantilever.vtu
```

Expected console log:

```text
[parser]  Read 33 nodes, 20 elements from tests/inputs/cantilever.inp
[model]   Dimension = 2, active materials = 1
[preproc] DOFs total = 66, free = 60, prescribed = 6
[solve]   SimplicialLDLT solve OK, residual norm = 5.0e-12
[post]    Max nodal displacement magnitude |u|_max = 1.358e-04 m
[vtu]     Wrote solution_cantilever.vtu
```

Open the `.vtu` in ParaView and apply **Warp By Vector** with
`Vectors = displacement` to see the deformed shape.

### Benchmark decks in `tests/inputs/`

| Deck | Element | What it validates |
| :--- | :--- | :--- |
| `cantilever.inp` | CPS4 | tip deflection vs analytical beam |
| `plate_hole.inp` | CPS4 | stress concentration factor near a hole |
| `cube3d.inp` | C3D8 | 3D uniaxial compression, exact σ_zz |
| `q4_pressure.inp` | CPS4 | top-face pressure, exact σ_yy |
| `h8_pressure.inp` | C3D8 | top-face pressure, exact σ_zz |

## CLI reference

```
femsolver --in <deck.inp> --out <solution.vtu> [--verbose]
```

Exit codes (doc 01 §7.1):

| Code | Class |
| :--- | :--- |
| 0 | success |
| 1 | input parse error (`InputError`) |
| 2 | model semantics error (`ModelError`) |
| 3 | numerical failure (reported by femcore) |
| 4 | linear solve did not converge (`SolverFailure`) |
| 5 | internal error |

## Development

See the task list targets (`make help`) or run steps individually:

```bash
make build      # pip install . -v
make test       # pytest python/tests -m "not slow"
make benchmark  # python tests/run_benchmarks.py
make lint       # flake8 + black --check
make format     # black + clang-format in place
make verify     # build + registry sanity + tests + benchmarks
```

### Verifying the element registry

```bash
python -c "import femcore; print(femcore.registered_elements())"
```

Expected:

```text
['C3D4', 'C3D8', 'CPE3', 'CPE4', 'CPS3', 'CPS4']
```

An empty list means the linker stripped the static registrations — make sure
every `cpp/src/elements_*.cpp` is listed in `CMakeLists.txt`.

## Repository layout

```
cpp/            C++ numerical core (M4–M9): kernels, assembler, solver, bindings
python/fem/     Python orchestration (M1–M3, M10, CLI): parser, builder, preprocessor, writers
docs/           Authoritative specifications (01 pipeline … 07 consistency audit)
tests/inputs/   Validated ABAQUS decks; tests/run_benchmarks.py analytical sweep
python/tests/   pytest suite (unit + integration)
.github/        CI workflow (Linux, Windows, lint)
```

## Documentation

The specs in `docs/` are the ground truth:

- `01` master architecture & pipeline, module contracts, error classes
- `02` element kernel layer & formulations (full math + shipped kernel code)
- `03` data structures & memory schemas, `femcore` binding API
- `04` `.inp` subset grammar & `.vtu` XML spec
- `05` build/execution guides & troubleshooting
- `06` testing, ADRs, analytical benchmarks, CI pipeline
- `07` module consistency audit & open items

## License

MIT. See [LICENSE](LICENSE).