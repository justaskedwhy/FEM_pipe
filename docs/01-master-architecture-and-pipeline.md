# Hybrid Python-C++ Finite Element Solver: Master Architecture & Pipeline Specification

**Document Package:** 1 of 6  
**Version:** 1.0  
**Target Platforms:** Windows 10/11, Linux (Ubuntu 20.04+)  
**Language Standards:** Python 3.10+, C++17  

---

## 1. Executive Summary & Core Architectural Principles

This document provides the complete architectural specification for the hybrid Python–C++ finite element (FE) solver. The system is designed to solve 2D (plane stress and plane strain) and 3D linear elastostatic problems with isotropic materials from ABAQUS `.inp` input files, producing VTK XML (`.vtu`) output files for visualization in ParaView.

### 1.1 Architectural Division of Labor
The core design philosophy enforces a strict division of responsibility:
* **Python Orchestrates:** Python handles file I/O, text parsing, mesh validation, global DOF mapping, post-processing formatting, and output file generation.
* **C++ Calculates:** C++ performs all numerical computations, element matrix generation, global assembly, boundary condition application, linear system solving, and stress/strain recovery.

```
┌─────────────────────────────────────────────────────────────┐
│  Python Layer: Orchestration, I/O, Parsing, Data Munging    │
└──────────────────────────────┬──────────────────────────────┘
                               │ Flat Contiguous NumPy Buffers
═══════════════════════════════▼═══════════════════════════════
                 PYTHON ↔ C++ INTERFACE (pybind11)
═══════════════════════════════╦═══════════════════════════════
                               │ C++ Eigen Primitive Structures
┌──────────────────────────────▼──────────────────────────────┐
│  C++ Layer: Numerical Core, Element Kernels, Assembly, Solve│
└─────────────────────────────────────────────────────────────┘
```

### 1.2 Core Architectural Principles (P1–P7)
* **P1: One Clean Language Boundary.** Exactly two boundary crossings occur per analysis run: one Python-to-C++ transfer (`SolverInput`) and one C++-to-Python transfer (`SolverOutput`). There is zero per-element chatter across languages.
* **P2: Elements Are Plugins.** All element-specific mathematical behavior, shape functions, integration rules, face topologies, and VTK cell codes are contained inside plugin kernels inheriting from `ElementKernel`.
* **P3: Unidirectional Data Flow.** Data flows strictly forward through standard intermediate data models: `.inp` → `RawModel` → `FEMModel` → `SolverInput` → `SolverOutput` → `.vtu`.
* **P4: Flat Numeric Buffers at Boundary.** All physical arrays crossing between Python and C++ are 1D/2D flat, contiguous, typed NumPy/Eigen buffers (`float64` and `int32`).
* **P5: Fail Loudly.** Every module validates its inputs against its contract and throws explicit exceptions upon violation. Silent failures or fallbacks are forbidden.
* **P6: Immutable Global State.** The only global state is the `ElementRegistry` singleton, which is populated at static initialization and immutable during runtime.
* **P7: ASCII Output First.** Output VTK XML files are produced in ASCII for human inspectability and version control diffing.

---

## 2. Binding Architectural Invariants (I1–I9)

To ensure long-term maintainability and modularity, the codebase strictly enforces nine binding invariants:

| # | Invariant Rule | Description |
|---|---|---|
| **I1** | **No Switch Outside Kernels** | No `switch` or `if-else` branching on element types exists outside an `ElementKernel` subclass. |
| **I2** | **No Hardcoded Metadata** | Parameters like `npe`, `dim`, `n_gauss`, `n_faces`, or `vtk_code` are never hardcoded outside kernels or the Python metadata cache. |
| **I3** | **CSR Connectivity Layout** | Mesh connectivity is stored as Compressed Sparse Row (CSR) flat arrays + offsets, avoiding zero-padded rectangular matrices. |
| **I4** | **ABAQUS String Naming** | Elements are identified strictly by their official ABAQUS string names (e.g., `"CPS4"`, `"C3D8"`) up to the C++ core. |
| **I5** | **Kernel as Single Source of Truth** | Shape functions, quadrature rules, strain matrices, face node maps, and VTK ordering are defined solely by the kernel. |
| **I6** | **Single-Line Element Addition** | Adding a new element requires creating its `.cpp` kernel file and adding exactly one line to `CMakeLists.txt`. |
| **I7** | **Metadata Cache Rule** | Python's `ELEMENT_META` dictionary is an ephemeral cache populated by querying the C++ registry via `femcore.element_meta()`. |
| **I8** | **Strict Modular Pipeline Contract** | The output of Module $M_i$ forms the exact input of Module $M_{i+1}$ without skipping stages. |
| **I9** | **Typed Contiguous Arrays** | All array buffers crossing the language boundary are typed (`float64` / `int32`) and contiguous in memory. |

---

## 3. Layered System Architecture

The solver is structured into six horizontal layers, where any layer may only depend on layers directly below it:

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 5 — CLI / Orchestration                                  │
│  python/fem/cli.py                                              │
│  Responsibility: Parse CLI arguments, orchestrate execution.    │
└────────────────────────────────┬────────────────────────────────┘
                                 │
┌────────────────────────────────▼────────────────────────────────┐
│  Layer 4 — I/O Layer                                            │
│  M1 Input Parser (in)    │  M10 VTK Writer (out)                │
│  Responsibility: Read .inp syntax and format .vtu XML.          │
└────────────────────────────────┬────────────────────────────────┘
                                 │
┌────────────────────────────────▼────────────────────────────────┐
│  Layer 3 — Python Data Model                                    │
│  M2 Model Builder        │  M3 Preprocessor / DOF Manager       │
│  Responsibility: Build semantic model, map global DOFs, CSR.    │
└────────────────────────────────┬────────────────────────────────┘
                                 │
═════════════════════════════════▼═════════════════════════════════
                 PYTHON ↔ C++ BOUNDARY (bindings.cpp)
═════════════════════════════════╦═════════════════════════════════
                                 │
┌────────────────────────────────▼────────────────────────────────┐
│  Layer 2 — C++ Numerical Pipeline Orchestrator                  │
│  fem_solve() entry point in cpp/src/solver.cpp                  │
│  Responsibility: Drive element loop, assembly, solve, post-proc.│
└────────────────────────────────┬────────────────────────────────┘
                                 │
┌────────────────────────────────▼────────────────────────────────┐
│  Layer 1 — C++ Numerical Primitives                             │
│  M4 Kernels │ M5 Loads │ M6 Assembler │ M7 BC │ M8 Solver │ M9 Post│
└────────────────────────────────┬────────────────────────────────┘
                                 │
┌────────────────────────────────▼────────────────────────────────┐
│  Layer 0 — External Libraries                                   │
│  Eigen (Sparse/Dense Linear Algebra) │ pybind11 (Language Glue) │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. The Python ↔ C++ Language Boundary

Data exchange across the language boundary is strictly encapsulated in two C++ structs exposed to Python via `pybind11`: `SolverInput` (Python → C++) and `SolverOutput` (C++ → Python).

### 4.1 SolverInput Schema (Python → C++)

| Array Name | Data Type | Buffer Shape | Semantic Description |
|---|---|---|---|
| `coords` | `float64` | `(n_nodes, 3)` | Physical node coordinates (always 3D; $z=0$ for 2D). |
| `elem_type_name` | `std::string` list | `(n_elem,)` | Official ABAQUS string name for each element. |
| `elem_conn_flat` | `int32` | `(sum_npe,)` | Flat 0-based concatenated node connectivity. |
| `elem_conn_offsets` | `int32` | `(n_elem + 1,)` | CSR offset pointers into `elem_conn_flat`. |
| `elem_mat` | `int32` | `(n_elem,)` | Material index assigned to each element. |
| `mat_E` | `float64` | `(n_mat,)` | Young's modulus per material ($Pa$). |
| `mat_nu` | `float64` | `(n_mat,)` | Poisson's ratio per material. |
| `mat_thickness` | `float64` | `(n_mat,)` | Thickness per material (2D only; 1.0 for 3D). |
| `mat_plane_mode` | `int32` | `(n_mat,)` | Analysis mode: `0` = Plane Stress, `1` = Plane Strain, `2` = 3D. |
| `dof_map` | `int32` | `(n_nodes, 3)` | Global DOF index per node/direction (`-1` if inactive). |
| `n_dofs_total` | `int32` | Scalar | Total number of global DOFs. |
| `n_dofs_free` | `int32` | Scalar | Total number of unconstrained global DOFs. |
| `prescribed_dofs` | `int32` | `(n_presc,)` | Global indices of constrained DOFs. |
| `prescribed_vals` | `float64` | `(n_presc,)` | Prescribed displacement magnitudes. |
| `point_load_dofs` | `int32` | `(n_ploads,)` | Global DOF indices receiving point loads. |
| `point_load_vals` | `float64` | `(n_ploads,)` | Point load magnitudes ($N$). |
| `pressure_elem` | `int32` | `(n_pl,)` | Element indices receiving face pressure. |
| `pressure_face` | `int32` | `(n_pl,)` | 1-based local face ID receiving pressure. |
| `pressure_val` | `float64` | `(n_pl,)` | Pressure magnitude ($Pa$; positive = compression). |

### 4.2 SolverOutput Schema (C++ → Python)

| Array Name | Data Type | Buffer Shape | Semantic Description |
|---|---|---|---|
| `displacement` | `float64` | `(n_nodes, 3)` | Calculated nodal displacement vector $(u_x, u_y, u_z)$. |
| `reaction` | `float64` | `(n_nodes, 3)` | Calculated nodal reaction forces at constrained DOFs. |
| `elem_stress` | `float64` | `(n_elem, 6)` | Centroidal stress tensor in Voigt order $[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx}]$. |
| `elem_strain` | `float64` | `(n_elem, 6)` | Centroidal strain tensor in Voigt order $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{zx}]$. |
| `elem_von_mises` | `float64` | `(n_elem,)` | Equivalent von Mises scalar stress field. |
| `elem_principal` | `float64` | `(n_elem, 3)` | Principal stresses sorted in descending order $(\sigma_1 \ge \sigma_2 \ge \sigma_3)$. |
| `converged` | `bool` | Scalar | Boolean flag indicating numerical solution success. |
| `residual_norm` | `float64` | Scalar | Equilibrium residual over the free-DOF partition only: $\|K_{ff}U_f + K_{fc}U_c - F_f\|/\|F_f\|$ (absolute when $\|F_f\|=0$); constrained-DOF reactions excluded. |

---

## 5. Detailed Module Specifications (M1–M10)

### M1 — Input Parser (Python)
* **Purpose:** Reads raw text from an ABAQUS `.inp` file and constructs an unvalidated `RawModel`.
* **Responsibilities:** Tokenizes keywords (`*NODE`, `*ELEMENT`, `*MATERIAL`, `*ELASTIC`, `*SOLID SECTION`, `*BOUNDARY`, `*CLOAD`, `*DLOAD`, `*STEP`, `*STATIC`, `*END STEP`); strips comments (`**`) and blank lines; parses data lines. Unsupported keywords log a warning and skip data blocks.
* **Output:** `RawModel` dataclass.

### M2 — FEM Model Builder (Python)
* **Purpose:** Transforms raw text entities into a validated, semantically consistent `FEMModel`.
* **Responsibilities:** Verifies node uniqueness; validates element connectivity; maps elements to materials via `ELSET` / `SOLID SECTION`; queries C++ registry to verify element existence and dimension; builds 0-based node and element indices.
* **Output:** `FEMModel` dataclass.

### M3 — Preprocessor / DOF Manager (Python)
* **Purpose:** Converts the physical FEM model into flat numeric C++ solver input buffers.
* **Responsibilities:** Assigns global DOFs; constructs `dof_map`; encodes element connectivity into CSR arrays (`elem_conn_flat`, `elem_conn_offsets`); partitions free and prescribed DOFs; populates `SolverInput`.
* **Output:** `SolverInput` object.

### M4 — Element Kernel Layer (C++)
* **Purpose:** Provides pluggable numerical evaluations for element formulations.
* **Responsibilities:** Implements `ElementKernel` interface; evaluates element stiffness $K_e = \int B^T D B \, d\Omega$; calculates face pressure vectors $F_e^{(p)}$; supplies $B$ matrix at arbitrary coordinates; specifies face topology and VTK node orderings.
* **Output:** `ElementMatrices` struct ($K_e, B, D$, Gauss point coords/weights, $\det J$).

### M5 — Load Engine (C++)
* **Purpose:** Converts external loads into global force vector contributions.
* **Responsibilities:** Maps point loads directly to global DOFs; delegates pressure loads to `kernel->pressure_force()` and scatters force vectors to global $F$.
* **Output:** Global load vector $F$.

### M6 — Global Assembler (C++)
* **Purpose:** Assembles element stiffness matrices into the sparse global stiffness matrix $K$.
* **Responsibilities:** Iterates over elements, creates kernels via `ElementRegistry`, builds local `ElementContext`, invokes `compute()`, maps local DOFs to global DOFs, and populates `Eigen::SparseMatrix<double>` via triplet list.
* **Output:** Global sparse symmetric matrix $K$.

### M7 — Boundary Condition Handler (C++)
* **Purpose:** Applies prescribed displacement boundary conditions and reduces the linear system.
* **Responsibilities:** Partitions DOFs into free ($f$) and constrained ($c$); forms reduced system $K_{ff} U_f = F_f - K_{fc} U_c$; handles zero and non-zero prescribed displacements.
* **Output:** `ReducedSystem` struct ($K_{ff}, F_r$, free/constrained DOF indices).

### M8 — Linear Solver (C++)
* **Purpose:** Solves the reduced sparse linear system $K_{ff} U_f = F_r$.
* **Responsibilities:** Factorizes $K_{ff}$ using `Eigen::SimplicialLDLT` (sparse Cholesky); falls back to `Eigen::SparseLU` if LDLT fails; reconstructs full displacement vector $U$; computes the equilibrium residual over the free-DOF partition $\|K_{ff}U_f + K_{fc}U_c - F_f\|/\|F_f\|$ (absolute when $\|F_f\|=0$), excluding constrained-DOF reactions.
* **Output:** Full global displacement vector $U$.

### M9 — Post-Processor (C++)
* **Purpose:** Computes strain, stress, and derived stress invariants.
* **Responsibilities:** Calculates element strain $\varepsilon = B U_e$ and stress $\sigma = D \varepsilon$ at element centroids; expands 2D plane stress/strain to 3D Voigt tensors; evaluates von Mises stress $\sigma_{\text{vm}}$; calculates principal stresses $\sigma_1 \ge \sigma_2 \ge \sigma_3$ via 3x3 eigensolver.
* **Output:** Populated `SolverOutput` fields.

### M10 — VTK Writer (Python)
* **Purpose:** Exports mesh and numerical results to a VTK XML UnstructuredGrid (`.vtu`) file.
* **Responsibilities:** Queries C++ element metadata for VTK cell type codes; formats point coordinates, cell connectivity, offsets, and cell types; writes point data (`displacement`, `displacement_magnitude`) and cell data (`stress`, `strain`, `von_mises`, `principal_stress`) in ASCII XML format.
* **Output:** `.vtu` file on disk.

---

## 6. Execution Workflow & Sequence Diagrams

### 6.1 End-to-End Analysis Workflow
```text
  User      CLI         M1         M2         M3       femcore     M4..M9     M10
   │         │          │          │          │          │           │         │
   │ argv ──►│          │          │          │          │           │         │
   │         │ parse ──►│          │          │          │           │         │
   │         │          │ read ────│          │          │           │         │
   │         │          │ RawModel►│          │          │           │         │
   │         │          │          │ validate │          │           │         │
   │         │          │          │ FEMModel►│          │           │         │
   │         │          │          │          │ DOFs     │           │         │
   │         │          │          │          │ CSR      │           │         │
   │         │          │          │          │ SolverIn►│           │         │
   │         │          │          │          │          │ solve ───►│         │
   │         │          │          │          │          │           │ kernels │
   │         │          │          │          │          │           │ assemble│
   │         │          │          │          │          │           │ BC      │
   │         │          │          │          │          │           │ solve   │
   │         │          │          │          │          │           │ stress  │
   │         │          │          │          │          │◄─SolverOut─│         │
   │         │◄─────────SolverOutput─────────────────────│           │         │
   │         │ write_vtu ────────────────────────────────────────────►│         │
   │         │          │          │          │          │           │  .vtu   │
   │◄─exit ──│          │          │          │          │           │         │
```

### 6.2 Element Assembly Loop (C++ Core)
```text
   Assembler (M6)        Registry (M4)         Kernel (M4)
        │                     │                     │
        │ for e in 0..n-1:    │                     │
        │ create(name[e]) ───►│                     │
        │                     │ factory() ─────────►│ new kernel
        │◄─────────unique_ptr─│                     │
        │                                           │
        │ build ElementContext(e)                   │
        │ compute(ctx) ─────────────────────────────►│
        │                                           │ evaluate shape fns
        │                                           │ compute J, detJ
        │                                           │ build B, D
        │                                           │ Ke = sum B^T D B detJ w
        │◄─────────────────ElementMatrices──────────│
        │                                           │
        │ scatter Ke into global K                  │
        │                                           │
```

---

## 7. Error Handling, Failure Classes, & Runtime Invariants

### 7.1 System Failure Classification

| Error Class | Example Root Cause | System Response | Exit Code |
|---|---|---|---|
| **User Input Error** | Malformed float, invalid `.inp` syntax | Raise `InputError` in Python with line number | Exit `1` |
| **Model Semantic Error** | Element references non-existent node ID | Raise `ModelError` in Python with entity IDs | Exit `2` |
| **Numerical Error** | Negative or zero Jacobian determinant | Throw `std::runtime_error` in C++ → Python `RuntimeError` | Exit `3` |
| **Solver Failure** | Singular matrix (unconstrained rigid modes) | Catch in `fem_solve()`, set `converged=false` | Exit `4` |
| **Internal Violation** | Invariant rule broken (e.g. array mismatch) | Debug `assert()` / throw `FemError` | Exit `5` |

### 7.2 Pre-Solve and Post-Solve Runtime Invariants

**Input Invariants (Checked at `fem_solve()` entry):**
1. `coords.shape == (n_nodes, 3)`
2. `elem_conn_offsets[0] == 0` and `elem_conn_offsets[-1] == len(elem_conn_flat)`
3. All `elem_type_name` entries exist in `ElementRegistry::registered_names()`
4. All material indices obey `0 <= elem_mat[e] < n_mat`
5. $E > 0$ and $-1.0 < \nu < 0.5$ for all materials

**Output Invariants (Checked before returning `SolverOutput`):**
1. `displacement.shape == (n_nodes, 3)`
2. `elem_stress.shape == (n_elem, 6)` and `elem_strain.shape == (n_elem, 6)`
3. `elem_principal[:, 0] >= elem_principal[:, 1] >= elem_principal[:, 2]`
4. If `converged == true`, `residual_norm < 1e-8`
