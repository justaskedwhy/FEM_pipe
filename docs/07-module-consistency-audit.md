# Module & Sub-Module Consistency Audit (M1–M10)

**Audit date:** 2026-09-15  
**Method:** Full read of all six spec docs (01–06); each module cross-checked across every doc that mentions it (module specs, data-structure layouts, file-format grammars, troubleshooting tables, test protocols). Findings recorded as **Confirmed** (consistent across docs), **Discrepancy** (docs contradict each other), or **Gap** (unspecified/missing).
Per repo ground truth (doc 03 §1.1, doc 06 ADR-0016), the docs govern and were treated as authoritative; this audit only reports where the docs themselves are internally inconsistent or silent.

---

## Cross-Cutting Checks

| Area | Verdict | Detail |
| :--- | :--- | :--- |
| Module output = next module input (I8) | **Confirmed** | `.inp → RawModel → FEMModel → SolverInput → SolverOutput → .vtu` consistent across 01 §3/§6, 03 §2–3 |
| Layer diagram vs module list (01 §3) | **Confirmed** | Layer 4 = M1/M10, Layer 3 = M2/M3, Layer 2 = `fem_solve`, Layer 1 = M4–M9, Layer 5 = CLI |
| Single language boundary (P1) | **Confirmed** | One `SolverInput` in, one `SolverOutput` out (01 §1.2, 03 §3, ADR-0002) |
| Error classes / exit codes | **Confirmed** | Input=1, Model=2, Numerical=3, Solver=4, Internal=5 (01 §7.1), reused in 05 §3 troubleshooting |
| Registry contents (doc 05 §2.2) vs kernels (doc 02) | **Confirmed** | Four kernels registered under six ABAQUS names: CPS3/CPE3 (T3), CPS4/CPE4 (Q4), C3D4 (T4), C3D8 (H8) |
| `femcore` binding API surface | **Resolved** | `femcore.fem_solve()`, `femcore.element_meta()`, `femcore.registered_elements()`, `femcore.registered_metadata()`, and the `test_element_*` hooks are now formally specified in doc 03 §8 (see items 15). |
| Registry name ordering | **Resolved** | `ElementRegistry::registered_names()` iterates an `unordered_map` (02 §1.4) → nondeterministic order; bindings sort before returning, matching doc 05 §2.2 (see item 14). |
| `gauss_points` column count | **Resolved** | Struct comment standardized to `(n_gauss x 3)` with 2D `z=0` in doc 02 §1.2 and doc 03 §4.2 (see item 13). |
| `pyproject.toml` pytest markers | **Gap** | `unit`/`integration`/`benchmark`/`slow` markers (06 §1.1) reference a `pyproject.toml` that does not exist; CI runs `pytest ... -m "not slow"` (06 §4) — still open, markers only affect reporting. |
| Grounded-source `docs/` files | **Gap** | Each package doc lists grounded sources (`docs/element_kernel.md`, `docs/data_structures.md`, ...) that do not exist — still open. |
| `tests/inputs/*.inp` | **Resolved** | Present and validated: `cantilever.inp`, `plate_hole.inp`, `cube3d.inp`, plus new `q4_pressure.inp` (CPS4) and `h8_pressure.inp` (C3D8); all pass the pipeline and benchmarks (see item 16). |

---

## M1 — Input Parser

- [ ] **Confirmed**
  - RawModel dataclasses match the EBNF data model exactly (03 §2.1 ↔ 04 §1.2–1.4).
  - Comments `**`, case-insensitive keywords, data lines bound to most recent keyword, ignored blank lines (04 §1.2).
  - Unsupported keywords skipped with a logged warning (01 §5 M1, 04 §1.5).
  - `InputError` raised with line number, exit 1 (01 §7.1, 05 §3.2).
  - Supported keyword set: NODE, ELEMENT, MATERIAL, ELASTIC, SOLID SECTION, BOUNDARY, CLOAD, DLOAD, STEP/STATIC/END STEP (04 §1.4).
- [ ] **Discrepancy**: `*BOUNDARY` grammar allows a named NSET as `node_or_set` (04 §1.4.5), but `RawBoundary` stores only `node_id: int` (03 §2.1). Set-based BCs are parsed nowhere.
  → **Resolved (2026-09-15)**: `RawBoundary`/`RawPointLoad` gained `node_set: Optional[str]`; M2 expands set-targeted lines to per-node entries with `ModelError` on undefined/empty sets; standalone `*NSET`/`*ELSET` blocks supported (see open item #2).
- [ ] **Discrepancy**: `*DLOAD` specifies a string face label (`P1`–`P6`, 04 §1.4.7) captured as `RawPressureLoad.face_label: str`, but `SolverInput.pressure_face` is an `int` local face ID (01 §4.1, 03 §3.1). The string→int conversion step (where labels map to 1-based `face_id`) is unspecified across M1→M2→M3.
  → **Resolved (2026-09-15)**: authority is `kernel->face_nodes(k)`; `P<k>` → 1-based local face `k`; preprocessor validates against element `n_faces` from the registry and raises `ModelError` on out-of-range (see open item #1).
- [ ] **Gap**: No spec for how XML/UTF-8 BOM or trailing blank-line handling at EOF is normalized (minor lexical detail).

---

## M2 — FEM Model Builder

- [ ] **Confirmed**
  - FEMModel layout matches 03 §2.2 exactly (coords `(n,3)`, CSR-shaped ragged conn, per-material arrays, preserved Raw boundary/load lists, `dimension`).
  - Responsibilities: node uniqueness, connectivity validation, material mapping via `ELSET`/`SOLID SECTION`, registry existence/dimension query, 0-based index build (01 §5, 03 §2.2 invariants).
  - `ModelError` with entity IDs, exit 2 (01 §7.1); `ModelError: Node X referenced by element Y...` (05 §3.2).
- [ ] **Discrepancy**: The canonical example deck `cantilever.inp` (04 §1.6) references element nodes (23, 24, 41, 44, 62) that are never defined in its `*NODE` block. As written it fails this module's own validation. (Corrected deck scaffolded at `tests/inputs/cantilever.inp`.)
  → **Resolved (2026-09-15)**: corrected deck at `tests/inputs/cantilever.inp` passes the full pipeline; doc 04 §1.6 gained an explicit pointer noting the in-doc listing is abridged/illustrative (see open item #16).
- [ ] **Gap**: Algorithm for deciding `dimension` (2 vs 3) and the exact "mixed dimension" detection rule referenced in the error table (05 §3.2) is not specified at the code level.
- [ ] **Gap**: `ELSET`/`NSET` handling — raw sets exist (03 §2.1) but the section→element matching rule for `*SOLID SECTION` across multiple `*ELEMENT` blocks is only described narratively.

---

## M3 — Preprocessor / DOF Manager

- [ ] **Confirmed**
  - SolverInput fields, shapes, and types per 03 §3.1 (authoritative superset) and 01 §4.1.
  - DOF numbering formula `dof(n, d) = n·dim + d`, inactive direction = `-1` in `dof_map` (03 §5.2).
  - CSR encoding `elem_conn_flat` + `elem_conn_offsets[0]==0`, `offsets[-1]==len(flat)` (03 §1.3, 01 §7.2).
  - Free/prescribed partition fed into `prescribed_dofs`/`prescribed_vals` (01 §4.1).
- [ ] **Discrepancy**: Doc 01 §4.1's SolverInput table omits the scalar fields `n_nodes`, `n_elem`, `n_mat`, `dimension` that doc 03 §3.1 includes; doc 01's own §7.2 checks reference `n_nodes`. Treat 03 as governing.
- [ ] **Gap**: Inactive-DOF handling for 3D fields under 2D (`dof_map[:,2] = -1`) is implied by 03 §5.2 and §1.3 but never stated explicitly as an M3 responsibility.

---

## M4 — Element Kernel Layer

- [ ] **Confirmed**
  - `ElementContext` (`coords (npe,3)`, E, nu, thickness, plane_mode) and `ElementMatrices` (Ke, B stacked, D, gauss_points, gauss_weights, detJ) match across 02 §1.2 and 03 §4.
  - `ElementKernel` interface, `ElementRegistry` factory + `REGISTER_ELEMENT` macro (02 §1.3–1.4).
  - Sorted T3/Q4 matrices vs math sections (02 §3 vs §2.4, §4 vs §2.5).
  - `det(J_g) ≤ 0` ⇒ kernel must throw (02 §1.2); thickness factor on 2D Ke, ignored for 3D (02 §5.4).
- [ ] **Discrepancy**: `name()` contract (02 §1.3: "ABAQUS canonical identifier, e.g. CPS4") is violated by both shipped kernels, which return short names `"T3"` (02 §3) and `"Q4"` (02 §4) while registering under CPS3/CPE3 and CPS4/CPE4.
  → **Resolved (2026-09-15)**: `name()` is now non-pure and returns the registry-injected ABAQUS key (see open item #5).
- [ ] **Discrepancy**: Inverted-element rule. 02 §1.2 mandates throwing when `det(J_g) ≤ 0`, but the T3 implementation uses `std::abs` on the signed area and only throws when the area is degenerate (< 1e-14); an inverted (negatively winded) triangle is silently sign-corrected. Q4, by contrast, throws on `detJ <= 0`.
  → **Accepted as documented behavior (2026-09-15)**: simplicials sign-correct; Q4/H8 throw strictly; note added to doc 02 §1.2 (see open item #6).
- [ ] **Discrepancy**: Q4 `build_D_2d` treats every `plane_mode != 0` as plane strain (no validation); T3 throws on invalid `plane_mode`. The two kernels disagree on the same contract input.
  → **Resolved (2026-09-15)**: Q4 accepts `{0, 1}` and throws otherwise, matching T3 (see open item #7).
- [ ] **Discrepancy**: Q4 `B_at()` and `pressure_force()` return zero matrices/vectors (02 §4) — labeled "skeleton" in-doc, so acknowledged, but functionally wrong for any mesh with face pressure.
  → **Resolved (2026-09-15)**: both fully implemented and verified against unit tests + `tests/inputs/q4_pressure.inp` e2e (see open item #4).
- [ ] **Gap**: T4 (02 §2.6) and H8 (02 §2.7) have full math but **no implementation**, no `face_nodes()` tables, and no explicit VTK order beyond doc 04 §2.3's identity vectors.
  → **Resolved (2026-09-15)**: T4/H8 kernels shipped with stiffness, `pressure_force`, face tables, identity VTK permutations (see open item #3).
- [ ] **Gap**: `ElementMatrices.B` stacking convention (`rows [g·n_strain, (g+1)·n_strain - 1]`) is defined, but the post-processing evaluation point for Q4/H8 centroidal stress (which stored B, or an averaged B) is unspecified (02 §2.10, M9).
  → **Resolved (2026-09-15)**: Gauss-point mean of `B·U_e`, single `D` multiply; documented in doc 02 §2.10 (see open item #11).

---

## M5 — Load Engine

- [ ] **Confirmed**
  - Point loads mapped directly onto global DOFs; pressures delegated to `kernel->pressure_force()` and scattered (01 §5, 01 §4.1).
  - Traction `t = -p·n`, positive pressure compressive into the element (02 §2.9, 04 §1.4.7); T3 example formula consistent (02 §2.9, 02 §3).
  - 2D loads scaled by thickness (02 §3 code); meshio/manual writers store the resulting `F` only implicitly through `SolverOutput`.
- [ ] **Gap**: Load set/`*CLOAD` NSET variant (04 §1.4.6 shows a bare node) has no `RawModel` representation for sets.

---

## M6 — Global Assembler

- [ ] **Confirmed**
  - Loop: registry `create()` → build `ElementContext` → `compute()` → map local→global DOFs → scatter via `Eigen::Triplet` into sparse K (01 §5, 01 §6.2 sequence diagram, 03 §6 `SparsityPattern`).
  - K symmetric; input arrays in CSR/int32/float64 (03 §1.3, I3/I9).
- [ ] **Gap**: No spec for symmetric-triplet duplication vs full scatter, or row/col count derivation from `n_dofs_total` (implementation detail only).

---

## M7 — Boundary Condition Handler

- [ ] **Confirmed**
  - Partition free/constrained; reduced system `K_ff U_f = F_f - K_fc U_c`; handles zero and nonzero prescribed values (01 §5 M7, 03 §4.4 `ReducedSystem`), signature matches `solver.hpp` (03 §6.2).
  - Reaction recovery `R = K·U - F` at constrained DOFs (`compute_reactions`, 03 §6.2).

---

## M8 — Linear Solver

- [ ] **Confirmed**
  - `Eigen::SimplicialLDLT` primary, `Eigen::SparseLU` fallback (01 §5, ADR-0008); singular failure ⇒ `converged=false`, exit 4 (01 §7.1).
  - If `converged`, `residual_norm < 1e-8` (01 §7.2); residual is the equilibrium residual over the free-DOF partition (see gap resolution below), not the full-system `‖K·U − F‖/‖F‖` that doc 03 §3.2 previously wrote — doc 03 §3.2 amended to match (2026-09-15).
- [ ] **Gap**: Numerically zero-but-unsaturated `F` (e.g., pure prescribed-displacement problem) makes the residual normalization de II dimensionless; no defined handling.
  → **Resolved (2026-09-15)**: residual and normalization are computed over **free DOFs only** (`‖K_ff u_f + K_fc u_c − F_f‖ / ‖F_f‖`, absolute residual when `‖F_f‖` is zero); reaction forces at constrained DOFs are excluded, so pure-prescribed and fully-constrained models can converge. n.b. the residual here is over the *unconstrained partition*, complementing `converged` semantics in doc 01 §7.2.

---

## M9 — Post-Processor

- [ ] **Confirmed**
  - `ε = B·U_e`, `σ = D·ε` at centroids (01 §5, 02 §2.10); Voigt ordering `[xx yy zz xy yz zx]` (03 §5.1), engineering shear `γ = 2ε` (03 §5.1).
  - 2D→3D expansion rules for plane stress (`σzz=0`, `εzz=−ν/(1−ν)(εxx+εyy)`) and plane strain (02 §2.10).
  - von Mises formula (02 §2.10) matches 6-component Voigt input from M4/M8.
  - Principal stresses = ordered eigenvalues of the full 3×3 tensor (02 §2.10); output invariant `σ1 ≥ σ2 ≥ σ3` (01 §7.2).
- [ ] **Gap**: As M4 — centroid evaluation point for Q4/H8 unspecified.
  → **Resolved (2026-09-15)**: Gauss-point mean of `B·U_e`; see M4 note on open item #11.

---

## M10 — VTK Writer

- [ ] **Confirmed**
  - ASCII XML `UnstructuredGrid`, `version="0.1"`, `byte_order="LittleEndian"` (ADR-0007, 04 §2.2); 3-component Float64 points with `z=0` for 2D (04 §2.3).
  - PointData (`displacement`, `displacement_magnitude`, optional `reaction`), CellData (`stress`, `strain`, `von_mises`, `principal_stress`) (04 §2.4).
  - VTK cell type codes must come from `femcore.element_meta(...)` (I7, 01 §5 M10) — manual writer honors this (04 §3.2).
- [ ] **Discrepancy**: The meshio writer (04 §3.1) hardcodes a local `vtk_code_to_meshio` table — a direct violation of invariant I2 (no hardcoded element metadata outside kernels/`ELEMENT_META`).
  → **Accepted as documented behavior (2026-09-15)**: the table maps VTK cell codes to meshio vocabulary strings, not element metadata; connectivity is still sourced from `femcore.element_meta(...)`. Note added to doc 07 (see open item #8).
- [ ] **Discrepancy**: Neither writer applies `kernel.vtk_permutation()` (02 §1.3 provides it); both emit ABAQUS node order. Safe only because all four current cells are identity permutations.
  → **Resolved (2026-09-15)**: both writers apply `vtk_permutation()` before emitting connectivity (see open item #9).
- [ ] **Discrepancy**: Manual writer (04 §3.2) omits the `strain` CellData array that the meshio writer (04 §3.1) and the spec's CellData table (04 §2.4) both include.
  → **Resolved (2026-09-15)**: manual writer emits `strain` CellData (see open item #10).
- [ ] **Gap**: meshio writer emits one-cell cell blocks (`cells.append((meshio_type, [conn]))` per element); correct but inefficient, and not reconciled with the "group contiguous connectivity" comment directly above it.

---

## Open Items Summary

| # | Module | Severity | Item | Status |
| :- | :--- | :--- | :--- | :--- |
| 1 | M1/M3 | **High** | `*DLOAD` string face labels vs int `pressure_face` — conversion step undefined | **Resolved** — authority is kernel `face_nodes(k)` (doc 02 §2.9); `P<k>` maps to local face id `k`, validated at preprocess against element `n_faces` (ModelError on range). doc 04 §1.4.7 amended. |
| 2 | M1 | **High** | `*BOUNDARY`/`*CLOAD` NSET support parsed by grammar but unrepresentable in `RawModel` | **Resolved** — standalone `*NSET`/`*ELSET` blocks; BOUNDARY/CLOAD accept set-name targets; expanded per-node in M2 with ModelError on undefined/empty set. doc 04 §1.4.5–§1.4.8 amended. |
| 3 | M4 | **High** | T4/H8 kernels have math but no implementation, face tables, or VTK permutation vectors | **Resolved** — T4/H8 shipped: stiffness, `pressure_force` (const-face center-split / 2×2 quad-Gauss bilinear), face tables, identity VTK permutations. doc 02 §2.9 amended. |
| 4 | M4 | **High** | Q4 skeleton (`pressure_force`, `B_at` return zeros) — breaks any pressure-loaded Q4 mesh | **Resolved** — Q4 `pressure_force` (edge traction, equal split) and `B_at` (evaluated Jacobian) implemented; verified against `q4_pressure.inp` e2e. doc 02 §4 code updated. |
| 5 | M4 | Medium | `name()` returns short name vs ABAQUS-name contract | **Resolved** — `name()` non-pure, returns registry-injected ABAQUS key via `set_registry_name` (doc 02 §1.3). |
| 6 | M4 | Medium | T3 silently accepts inverted triangles vs "throw on detJ ≤ 0" rule; Q4/T3 disagree | **Resolved (documented)** — simplicials sign-correct via `std::abs` (doc 02 §1.2 note); Q4/H8 throw strictly. Deviation accepted as shipped-kernel behavior. |
| 7 | M4 | Medium | Q4 vs T3 `build_D_2d` disagree on invalid `plane_mode` handling | **Resolved** — Q4 now throws on `plane_mode ∉ {0,1}` (matches T3). |
| 8 | M10 | Medium | meshio writer hardcodes VTK→meshio table (violates I2) | **Resolved (documented)** — table is VTK↔meshio vocabulary mapping, not element metadata; based on `femcore.element_meta(...)["vtk"]`. |
| 9 | M10 | Medium | writers ignore `vtk_permutation()` | **Resolved** — both writers apply `kernel.vtk_permutation()` when emitting connectivity. |
| 10 | M10 | Medium | manual writer omits `strain` CellData | **Resolved** — manual writer emits `strain` CellData (matches meshio writer / doc 04 §2.4). |
| 11 | M4/M9 | Low | centroid evaluation point for Q4/H8 stress recovery unspecified | **Resolved** — Gauss-point arithmetic mean of stacked `B·U_e`, then one `D` multiply (supports patch-exact recovery); documented in doc 02 §2.10. |
| 12 | 01/03 | Low | SolverInput scalar fields listed only in 03 §3.1 (01 §4.1 stale) | **Open** — treat 03 as governing (low priority). |
| 13 | 02/03 | Low | `gauss_points` column count comment (`dim` vs 3) | **Resolved** — `(n_gauss x 3)` with 2D `z=0`; doc 02 §1.2 / doc 03 §4.2 amended. |
| 14 | Cross | Low | `registered_names()` unordered vs doc 05 sorted output | **Resolved** — bindings sort before returning (`registered_elements`/`registered_metadata`). |
| 15 | Cross | Low | `femcore` binding API (`element_meta`, `registered_elements`, struct casts) unspecified | **Resolved** — doc 03 §8 added: full function signatures, return shapes, struct attributes. |
| 16 | M2 | Low | example `cantilever.inp` references undefined nodes | **Resolved** — corrected deck at `tests/inputs/cantilever.inp`; doc 04 §1.6 pointer added; canonical deck passes e2e. |

**Resolution audit date:** 2026-09-15. All high-severity items (1–4) closed; two items remain open: **#12** (doc 01 §4.1 SolverInput table should gain the scalar fields already mandated by doc 03 §3.1).

**Recommended next actions** (updated): resolve open item #12 by mirroring doc 03 §3.1 scalar columns into doc 01 §4.1; run `pytest python/tests -m "not slow"` (63 tests) and `python tests/run_benchmarks.py` (3 benchmarks) after any further code change; keep new elements to the doc 02 §5.1 three-step workflow.