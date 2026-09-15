# AGENTS.md

## What this repo is

Spec-first repo for a hybrid **Python–C++ linear elastic FE solver** fed by ABAQUS `.inp` files, producing VTK `.vtu` output. The spec lives in `docs/`; a scaffolded but not-yet-built file system lives at `cpp/`, `python/`, `tests/`. There is no `git` or README yet.

## Ground truth

The seven numbered docs in `docs/` are the authoritative spec. Doc 03 §1.1 states: if code ever disagrees with these docs, **the docs govern**. Doc 06 ADR-0016 confirms documentation-first development.

- `docs/01` — pipeline & architecture: modules M1–M10, invariants I1–I9 (each module output = next module's exact input), error classes/exits
- `docs/02` — element kernel layer: `ElementKernel` interface, `ElementRegistry` + `REGISTER_ELEMENT` macro, math for T3/Q4/T4/H8, full T3 impl + Q4 skeleton
- `docs/03` — data structures + memory schemas: Python dataclasses, C++ structs (`SolverInput`, `SolverOutput`), complete `types.hpp`/`solver.hpp` layouts
- `docs/04` — `.inp` subset (EBNF grammar, supported keywords) + `.vtu` XML spec
- `docs/05` — build/execution guides + troubleshooting
- `docs/06` — testing, 3 analytical benchmarks, 16 ADRs, CI workflow
- `docs/07` — module consistency audit (M1–M10 recheck; 16 open items)

## Contract rules an agent will get wrong

- **One language boundary**: exactly two crossings per run — one `SolverInput` into `fem_solve()`, one `SolverOutput` back. No per-element Python↔C++ chatter.
- **Flat typed arrays across the boundary**: `float64` / `int32` only. Never `int64` or `float32`. All arrays contiguous.
- **Coordinates are always `(n, 3)`** — even 2D meshes store `z=0`. Element `coords` in `ElementContext` is `(npe, 3)`.
- **Element identity by ABAQUS string name** (e.g. `"CPS4"`, `"C3D8"`), never a numeric code. The `ElementCode` enum in doc 03 is bookkeeping only.
- **No `switch`/`if-else` on element type outside an `ElementKernel` subclass** (invariant I1). No hardcoded element metadata outside kernels or the Python `ELEMENT_META` cache (I2).
- **Suffix conventions** (doc 03 §1.2): `_id` = external 1-based user IDs, `_idx` = internal 0-based solver indices, `_flat`/`_offsets` = CSR connectivity, `_vals`/`_dofs` = value/DOF payloads. Ambiguous bare `id`/`index` are prohibited.
- **Adding an element** = new `cpp/src/elements_<name>.cpp` + exactly one line added to `CMakeLists.txt` (I6). Kernels self-register via the macro in the `.cpp` file.
- **VTK output is ASCII** (P7, ADR-0007). VTK cell codes come from `femcore.element_meta(...)`, not a hardcoded Python table.
- **Fail loudly** (P5): validate every contract; throw exceptions. C++ prints nothing (ADR-0014) — diagnostics live in `SolverOutput` + exceptions.

## Implementation notes from the specs

- Target structure is scaffolded: `cpp/src/*.cpp`, `cpp/include/*.hpp`, `cpp/bindings/bindings.cpp`, `python/fem/`, `python/tests/`, `tests/inputs/*.inp`. The C++ core is not yet compiled.
- Doc 02's T3/Q4 code uses `std::abs` on area so inverted elements throw rather than silently corrupting; keep that invariant.
- Solver: `Eigen::SimplicialLDLT` primary, `SparseLU` fallback; BC handler reduces `K_ff u_f = F_f - K_fc u_c`.
- Reference tests for new work: single-element patch test + mandatory Q4↔T3 swap test (doc 06 §1).

## Command conventions (intended, not yet runnable)

The docs specify these but nothing is buildable yet — no commit/verify commands exist until the C++ core compiles (`pip install . -v`):

```
pip install . -v                                  # builds femcore (pybind11)
python -c "import femcore; print(femcore.registered_elements())"   # registry sanity check
pytest python/tests -m "not slow" -v               # test suite
black python/fem && flake8 python/fem --max-line-length=100        # Python style
clang-format -i cpp/include/*.hpp cpp/src/*.cpp cpp/bindings/*.cpp # C++ style
```

Test decks live in `tests/inputs/`; benchmarks in `tests/run_benchmarks.py`.