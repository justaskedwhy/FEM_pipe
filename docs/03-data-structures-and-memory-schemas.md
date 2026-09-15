# Package 3: Data Structures & Memory Schemas Specification
**Document Version:** 1.0  
**Target Files:** `docs/data_structures.md`, `cpp/include/types.hpp`, `cpp/include/solver.hpp`  
**Grounded Sources:** `docs/data_structures.md`, `cpp/include/types.hpp`, `cpp/include/solver.hpp`, `docs/master.md`

---

## 1. System Naming, Shape & Type Conventions

### 1.1 Purpose & Authority
This document serves as the single source of truth for all data structures, memory layouts, and C++ header contracts in the hybrid Python-C++ finite element solver. If any discrepancy exists between code implementations and this specification, this specification governs.

### 1.2 Naming Conventions
To prevent indexing bugs across language boundaries, strict suffix rules are enforced:

| Suffix | Canonical Meaning | Example |
| :--- | :--- | :--- |
| `_id` | External identifier as defined in the ABAQUS `.inp` file (sparse, positive integers). | `node_id`, `elem_id` |
| `_idx` | Internal zero-based dense index used for solver computations ($0 \le \text{idx} < N$). | `node_idx`, `elem_idx` |
| `_map` | Lookup array or table mapping one index space to another. | `dof_map`, `g2r` |
| `_flat` | Concatenated 1D array representing ragged multidimensional data. | `elem_conn_flat` |
| `_offsets` | CSR-style cumulative offset array for indexing into `_flat` arrays. | `elem_conn_offsets` |
| `_vals` | Numerical payload array (e.g., prescribed displacement or load magnitudes). | `prescribed_vals`, `point_load_vals` |
| `_dofs` | Array containing global Degree-of-Freedom indices. | `prescribed_dofs`, `free_dofs` |
| `elem_*` | Per-element dense arrays (indexed by `elem_idx`). | `elem_mat`, `elem_stress` |
| `mat_*` | Per-material property arrays (indexed by `mat_idx`). | `mat_E`, `mat_nu` |

> **Rule:** `_id` values belong to user input; `_idx` values belong to solver computations. Ambiguous identifiers like `id` or `index` without qualifiers are prohibited.

---

### 1.3 Shape & Type Conventions
All numerical arrays crossing the Python ↔ C++ boundary or passed between modules must strictly comply with the following shape and typing invariants:

| Category | Shape Rule | Type | Invariant / Notes |
| :--- | :--- | :--- | :--- |
| **Coordinates** | $(N_{\text{nodes}}, 3)$ | `float64` / `double` | Always 3 columns. For 2D meshes, column index 2 ($z$) contains zeros. |
| **Displacements** | $(N_{\text{nodes}}, 3)$ | `float64` / `double` | Always 3 columns. For 2D meshes, column index 2 ($u_z$) contains zeros. |
| **Reactions** | $(N_{\text{nodes}}, 3)$ | `float64` / `double` | Always 3 columns. Unconstrained DOFs store $0.0$. |
| **Stresses / Strains** | $(N_{\text{elem}}, 6)$ | `float64` / `double` | Voigt 6 components: $[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx}]$. |
| **Connectivity** | CSR (`flat` + `offsets`) | `int32` / `int32_t` | `flat` size is $\sum n_{pe}$; `offsets` size is $N_{\text{elem}} + 1$. Never padded. |
| **DOF Map** | $(N_{\text{nodes}}, 3)$ | `int32` / `int32_t` | Global DOF index per node and direction; `-1` if direction is inactive. |
| **Integers** | Scaled / 1D / 2D | `int32` / `int32_t` | `int64` is never used across the C++ boundary. |
| **Floating-Point** | Scaled / 1D / 2D | `float64` / `double` | `float32` is never used. Precision is double everywhere. |

---

### 1.4 Parallel Index Spaces
The solver maintains 6 distinct index spaces:

| Index Space | Identifier | Type | Density | Description |
| :--- | :--- | :--- | :--- | :--- |
| **External Node ID** | `node_id` | Positive `int32` | Sparse | External ID from `.inp` file (e.g., 100, 200, 300). |
| **Internal Node Index** | `node_idx` | 0-based `int32` | Dense | Row index in `coords` array ($0 \le \text{node\_idx} < N_{\text{nodes}}$). |
| **External Element ID** | `elem_id` | Positive `int32` | Sparse | External ID from `.inp` file (e.g., 1, 5, 12). |
| **Internal Element Index** | `elem_idx` | 0-based `int32` | Dense | Element position in `elem_type_name` ($0 \le \text{elem\_idx} < N_{\text{elem}}$). |
| **Material Index** | `mat_idx` | 0-based `int32` | Dense | Material position in `mat_E` array ($0 \le \text{mat\_idx} < N_{\text{mat}}$). |
| **Global DOF Index** | `dof_idx` | 0-based `int32` | Dense | Position in global stiffness matrix $K$ ($0 \le \text{dof\_idx} < N_{\text{dofs\_total}}$). |

---

## 2. Python Data Structures (M1 & M2)

### 2.1 M1 Raw Model Dataclasses (`python/fem/parser.py`)
`RawModel` represents the direct, unvalidated output of the ABAQUS `.inp` parser (Module M1).

```python
from dataclasses import dataclass, field

@dataclass
class RawNode:
    id: int          # Positive external node ID
    x: float        # x-coordinate
    y: float        # y-coordinate
    z: float = 0.0  # z-coordinate (defaults to 0.0 for 2D)

@dataclass
class RawElement:
    id: int               # Positive external element ID
    type_name: str        # Canonical ABAQUS name (e.g., 'CPS4', 'C3D8')
    node_ids: list[int]   # External node IDs in file order

@dataclass
class RawMaterial:
    name: str  # Material identifier (from *MATERIAL, NAME=)
    E: float   # Young's modulus (Pa)
    nu: float  # Poisson's ratio

@dataclass
class RawSection:
    elset_name: str     # Target element set name
    material_name: str  # Referenced material name
    thickness: float    # Thickness for 2D (ignored/1.0 for 3D)

@dataclass
class RawBoundary:
    node_id: int    # External node ID
    dof_first: int  # 1-based start DOF (1=x, 2=y, 3=z)
    dof_last: int   # 1-based end DOF (inclusive)
    value: float    # Prescribed displacement value

@dataclass
class RawPointLoad:
    node_id: int     # External node ID
    dof: int         # 1-based load DOF (1=x, 2=y, 3=z)
    magnitude: float # Force magnitude (N)

@dataclass
class RawPressureLoad:
    elem_id: int        # External element ID
    face_label: str     # Local face label (e.g., 'P1', 'P2', 'P3')
    magnitude: float    # Pressure magnitude (Pa; positive = compressive into element)

@dataclass
class RawModel:
    nodes: dict[int, RawNode]             = field(default_factory=dict)
    elements: dict[int, RawElement]       = field(default_factory=dict)
    materials: dict[str, RawMaterial]     = field(default_factory=dict)
    sections: list[RawSection]            = field(default_factory=list)
    boundaries: list[RawBoundary]         = field(default_factory=list)
    point_loads: list[RawPointLoad]       = field(default_factory=list)
    pressure_loads: list[RawPressureLoad] = field(default_factory=list)
    nsets: dict[str, list[int]]           = field(default_factory=dict)
    elsets: dict[str, list[int]]          = field(default_factory=dict)
```

#### Invariants for `RawModel`
1. `RawNode.id > 0` and `RawElement.id > 0`.
2. No cross-reference validation is performed in `RawModel` (deferred to Module M2).
3. `1 <= RawBoundary.dof_first <= RawBoundary.dof_last <= 3`.

---

### 2.2 M2 Solver-Ready Model (`python/fem/model_builder.py`)
`FEMModel` represents the validated, 0-indexed mesh produced by Module M2.

```python
import numpy as np
from dataclasses import dataclass

@dataclass
class FEMModel:
    coords:         np.ndarray         # shape (n_nodes, 3), float64
    elem_type_name: list[str]          # length (n_elem,), str
    elem_conn:      list[list[int]]    # ragged list of 0-based node indices
    elem_mat:       np.ndarray         # shape (n_elem,), int32
    mat_E:          np.ndarray         # shape (n_mat,), float64
    mat_nu:         np.ndarray         # shape (n_mat,), float64
    mat_thickness:  np.ndarray         # shape (n_mat,), float64
    mat_plane_mode: np.ndarray         # shape (n_mat,), int32 (0=stress, 1=strain, 2=3D)
    boundaries:     list[RawBoundary]  # Preserved from RawModel
    point_loads:    list[RawPointLoad] # Preserved from RawModel
    pressure_loads: list[RawPressureLoad] # Preserved from RawModel
    n_nodes: int
    n_elem: int
    n_mat: int
    dimension: int                     # 2 or 3
```

#### Invariants for `FEMModel`
1. `coords.shape == (n_nodes, 3)` with `dtype=np.float64`.
2. `all(0 <= mat_idx < n_mat for mat_idx in elem_mat)`.
3. Every referenced node in `elem_conn` exists within $[0, n_{\text{nodes}} - 1]$.
4. `dimension` matches all registered element kernels in `elem_type_name`.

---

## 3. Python ↔ C++ Boundary Transfer Schemas

### 3.1 `SolverInput` (Python → C++ Payload)

`SolverInput` transfers the entire preprocessed model from Python to C++ in a single call.

| Field Name | C++ Type | Python / NumPy Equivalent | Array Shape | Description |
| :--- | :--- | :--- | :--- | :--- |
| `coords` | `Eigen::Matrix<double, Dynamic, 3>` | `np.ndarray(float64)` | $(N_{\text{nodes}}, 3)$ | Node coordinates ($z=0$ for 2D). |
| `elem_type_name` | `std::vector<std::string>` | `list[str]` | $(N_{\text{elem}},)$ | Canonical ABAQUS element names. |
| `elem_conn_flat` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(\sum n_{pe},)$ | Concatenated 0-based node indices. |
| `elem_conn_offsets` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{elem}} + 1,)$ | CSR offsets into `elem_conn_flat`. |
| `elem_mat` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{elem}},)$ | Material index per element. |
| `mat_E` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{mat}},)$ | Young's Modulus per material (Pa). |
| `mat_nu` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{mat}},)$ | Poisson's Ratio per material. |
| `mat_thickness` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{mat}},)$ | Thickness (2D only; 1.0 for 3D). |
| `mat_plane_mode` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{mat}},)$ | `0` = Plane Stress, `1` = Plane Strain, `2` = 3D. |
| `dof_map` | `Eigen::Matrix<int32_t, Dynamic, 3>` | `np.ndarray(int32)` | $(N_{\text{nodes}}, 3)$ | Global DOF index; `-1` if inactive. |
| `n_dofs_total` | `int32_t` | `int` | Scalar | Total number of global DOFs. |
| `n_dofs_free` | `int32_t` | `int` | Scalar | Number of unconstrained DOFs. |
| `prescribed_dofs` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{presc}},)$ | Prescribed global DOF indices. |
| `prescribed_vals` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{presc}},)$ | Prescribed displacement values (m). |
| `point_load_dofs` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{ploads}},)$ | Loaded global DOF indices. |
| `point_load_vals` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{ploads}},)$ | Point load force magnitudes (N). |
| `pressure_elem` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{pl}},)$ | Element index per pressure load. |
| `pressure_face` | `std::vector<int32_t>` | `np.ndarray(int32)` | $(N_{\text{pl}},)$ | Local face ID (1-based). |
| `pressure_val` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{pl}},)$ | Pressure magnitudes (Pa). |
| `n_nodes` | `int32_t` | `int` | Scalar | Total node count. |
| `n_elem` | `int32_t` | `int` | Scalar | Total element count. |
| `n_mat` | `int32_t` | `int` | Scalar | Total material count. |
| `dimension` | `int32_t` | `int` | Scalar | Spatial dimension (`2` or `3`). |

#### C++ Entry Validation Invariants
1. `elem_conn_offsets[0] == 0` and `elem_conn_offsets[n_elem] == elem_conn_flat.size()`.
2. `all(ElementRegistry::instance().has(name))` for all entries in `elem_type_name`.
3. `mat_E[i] > 0.0` and `-1.0 < mat_nu[i] < 0.5` for all materials $i$.

---

### 3.2 `SolverOutput` (C++ → Python Payload)

`SolverOutput` transfers all computed primary and secondary fields from C++ back to Python.

| Field Name | C++ Type | Python / NumPy Equivalent | Array Shape | Description |
| :--- | :--- | :--- | :--- | :--- |
| `displacement` | `Eigen::Matrix<double, Dynamic, 3>` | `np.ndarray(float64)` | $(N_{\text{nodes}}, 3)$ | Primary nodal displacements ($u_x, u_y, u_z$). |
| `reaction` | `Eigen::Matrix<double, Dynamic, 3>` | `np.ndarray(float64)` | $(N_{\text{nodes}}, 3)$ | Nodal reaction forces at constrained DOFs. |
| `elem_stress` | `Eigen::Matrix<double, Dynamic, 6>` | `np.ndarray(float64)` | $(N_{\text{elem}}, 6)$ | Element centroid stress (Voigt 6). |
| `elem_strain` | `Eigen::Matrix<double, Dynamic, 6>` | `np.ndarray(float64)` | $(N_{\text{elem}}, 6)$ | Element centroid strain (Voigt 6). |
| `elem_von_mises` | `Eigen::VectorXd` | `np.ndarray(float64)` | $(N_{\text{elem}},)$ | von Mises equivalent stress. |
| `elem_principal` | `Eigen::Matrix<double, Dynamic, 3>` | `np.ndarray(float64)` | $(N_{\text{elem}}, 3)$ | Principal stresses ($\sigma_1 \ge \sigma_2 \ge \sigma_3$). |
| `converged` | `bool` | `bool` | Scalar | `true` if linear solve succeeded. |
| `residual_norm` | `double` | `float` | Scalar | Equilibrium residual over the **unconstrained (free) partition**: $\frac{\|K_{ff}U_f + K_{fc}U_c - F_f\|}{\|F_f\|}$, or the absolute free residual when $\|F_f\| = 0$. Reaction forces at constrained DOFs are excluded (doc 01 §7.2, M8). |

---

## 4. C++ Internal Numerical Structures

### 4.1 `ElementContext`
Constructed by the element driver loop for each element $e$ and passed directly to `kernel->compute(ctx)`.

```cpp
namespace fem {

struct ElementContext {
    Eigen::Matrix<double, Eigen::Dynamic, 3> coords; // (npe x 3) local node coordinates
    double E;          // Young's Modulus (Pa)
    double nu;         // Poisson's Ratio
    double thickness;  // Element thickness (2D only; 1.0 for 3D)
    int plane_mode;    // 0 = Plane Stress, 1 = Plane Strain, 2 = 3D
};

} // namespace fem
```

### 4.2 `ElementMatrices`
Returned by `kernel->compute(ctx)` containing computed element-level numerical arrays.

```cpp
namespace fem {

struct ElementMatrices {
    Eigen::MatrixXd Ke;           // (ndof_e x ndof_e) dense element stiffness matrix
    Eigen::MatrixXd B;            // ((n_gauss * n_strain) x ndof_e) stacked strain-displacement matrix
    Eigen::MatrixXd D;            // (n_strain x n_strain) constitutive matrix
    Eigen::MatrixXd gauss_points; // (n_gauss x 3) natural coords; 2D stores z=0 (doc 02 §1.2)
    Eigen::VectorXd gauss_weights;// (n_gauss) Gauss point integration weights
    Eigen::VectorXd detJ;         // (n_gauss) Jacobian determinants at Gauss points
};

} // namespace fem
```

### 4.3 `ElementDofMap`
Maps element-local DOF indices ($0 \le i < \text{ndof}_e$) to global system matrix DOFs.

```cpp
namespace fem {

struct ElementDofMap {
    std::vector<int32_t> global_dofs; // size: (npe * dim)
};

} // namespace fem
```

### 4.4 `ReducedSystem`
Generated by Module M7 (BC Handler) after partitioning system matrix $K$ and force vector $F$.

```cpp
namespace fem {

struct ReducedSystem {
    Eigen::SparseMatrix<double> Kr;               // Free-free block of global stiffness matrix
    Eigen::VectorXd             Fr;               // Reduced force RHS vector (F_f - K_fc * U_c)
    std::vector<int32_t>        free_dofs;        // Global DOF indices of free DOFs
    std::vector<int32_t>        constrained_dofs; // Global DOF indices of constrained DOFs
    Eigen::VectorXd             constrained_vals; // Prescribed displacement values
};

} // namespace fem
```

---

## 5. Voigt Notation & DOF Mapping Conventions

### 5.1 Voigt Stress & Strain Ordering

#### 2D Plane Stress / Plane Strain (3 Components)
$$\boldsymbol{\sigma}_{2D} = \begin{bmatrix} \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{xy} \end{bmatrix}, \qquad \boldsymbol{\varepsilon}_{2D} = \begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \gamma_{xy} \end{bmatrix}$$

#### 3D Isotropic Elasticity (6 Components)
$$\boldsymbol{\sigma}_{3D} = \begin{bmatrix} \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \\ \sigma_{xy} \\ \sigma_{yz} \\ \sigma_{zx} \end{bmatrix}, \qquad \boldsymbol{\varepsilon}_{3D} = \begin{bmatrix} \varepsilon_{xx} \\ \varepsilon_{yy} \\ \varepsilon_{zz} \\ \gamma_{xy} \\ \gamma_{yz} \\ \gamma_{zx} \end{bmatrix}$$

* **Engineering Shear Rule**: Shear strain components in Voigt vectors are engineering shear strains:
  $$\gamma_{xy} = 2\varepsilon_{xy}, \quad \gamma_{yz} = 2\varepsilon_{yz}, \quad \gamma_{zx} = 2\varepsilon_{zx}$$

---

### 5.2 Global DOF Numbering Formula

DOFs are assigned in ascending node index order (`node_idx` $0 \le n < N_{\text{nodes}}$):

$$\text{dof}(n, d) = n \cdot \text{dim} + d \quad (d \in \{0, 1, 2\})$$

* **2D Model ($\text{dim} = 2$)**: Node $n$ owns DOFs $2n$ ($u_x$) and $2n+1$ ($u_y$). `dof_map(n, 2)` contains `-1`.
* **3D Model ($\text{dim} = 3$)**: Node $n$ owns DOFs $3n$ ($u_x$), $3n+1$ ($u_y$), and $3n+2$ ($u_z$).

---

## 6. Complete C++ Header Code

### 6.1 Header File 1: `cpp/include/types.hpp`

```cpp
// cpp/include/types.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace fem {

// =============================================================
// Internal Element Codes (Bookkeeping reference only)
// Canonical identifier is the ABAQUS string name.
// =============================================================
enum class ElementCode : int32_t {
    UNKNOWN = 0,
    T3      = 1,   // 2D, 3 nodes
    Q4      = 2,   // 2D, 4 nodes
    T4      = 3,   // 3D, 4 nodes
    H8      = 4    // 3D, 8 nodes
};

// =============================================================
// Analysis Mode
// =============================================================
enum class PlaneMode : int32_t {
    PLANE_STRESS = 0,
    PLANE_STRAIN = 1,
    THREE_D      = 2
};

// =============================================================
// SolverInput — The single Python -> C++ struct payload.
// =============================================================
struct SolverInput {
    // --- Geometry ---
    Eigen::Matrix<double, Eigen::Dynamic, 3> coords;   // (n_nodes, 3)

    // --- Element Connectivity (CSR Layout) ---
    std::vector<std::string>  elem_type_name;          // (n_elem,)
    std::vector<int32_t>      elem_conn_flat;          // (sum_npe,)
    std::vector<int32_t>      elem_conn_offsets;       // (n_elem + 1,)

    // --- Material Assignment ---
    std::vector<int32_t>      elem_mat;                // (n_elem,)
    Eigen::VectorXd           mat_E;                   // (n_mat,)
    Eigen::VectorXd           mat_nu;                  // (n_mat,)
    Eigen::VectorXd           mat_thickness;           // (n_mat,)
    std::vector<int32_t>      mat_plane_mode;          // (n_mat,)

    // --- DOF Map ---
    Eigen::Matrix<int32_t, Eigen::Dynamic, 3> dof_map; // (n_nodes, 3), -1 if inactive
    int32_t                   n_dofs_total = 0;
    int32_t                   n_dofs_free  = 0;

    // --- Boundary Conditions ---
    std::vector<int32_t>      prescribed_dofs;         // (n_presc,)
    Eigen::VectorXd           prescribed_vals;         // (n_presc,)

    // --- Point Loads ---
    std::vector<int32_t>      point_load_dofs;         // (n_ploads,)
    Eigen::VectorXd           point_load_vals;         // (n_ploads,)

    // --- Pressure Loads ---
    std::vector<int32_t>      pressure_elem;           // (n_pl,)
    std::vector<int32_t>      pressure_face;           // (n_pl,) 1-based
    Eigen::VectorXd           pressure_val;            // (n_pl,)

    // --- Metadata Scalars ---
    int32_t                   n_nodes   = 0;
    int32_t                   n_elem    = 0;
    int32_t                   n_mat     = 0;
    int32_t                   dimension = 0;           // 2 or 3
};

// =============================================================
// SolverOutput — The single C++ -> Python struct payload.
// =============================================================
struct SolverOutput {
    Eigen::Matrix<double, Eigen::Dynamic, 3> displacement;   // (n_nodes, 3)
    Eigen::Matrix<double, Eigen::Dynamic, 3> reaction;       // (n_nodes, 3)
    Eigen::Matrix<double, Eigen::Dynamic, 6> elem_stress;    // (n_elem, 6) Voigt
    Eigen::Matrix<double, Eigen::Dynamic, 6> elem_strain;    // (n_elem, 6) Voigt
    Eigen::VectorXd                          elem_von_mises;  // (n_elem,)
    Eigen::Matrix<double, Eigen::Dynamic, 3> elem_principal;  // (n_elem, 3)
    bool                                     converged = false;
    double                                   residual_norm = 0.0;
};

// =============================================================
// Internal Assembly Helpers
// =============================================================
struct ElementDofMap {
    std::vector<int32_t> global_dofs; // Global DOFs for local element DOFs
};

struct SparsityPattern {
    std::vector<Eigen::Triplet<double>> triplets;
    int rows = 0;
    int cols = 0;
};

// =============================================================
// Exception Helpers
// =============================================================
struct FemError : public std::runtime_error {
    explicit FemError(const std::string& msg)
        : std::runtime_error("[femcore] " + msg) {}
};

} // namespace fem
```

---

### 6.2 Header File 2: `cpp/include/solver.hpp`

```cpp
// cpp/include/solver.hpp
#pragma once

#include "types.hpp"
#include "element_registry.hpp"
#include <Eigen/Sparse>

namespace fem {

// =============================================================
// Top-Level Solver Entry Point (Exposed via pybind11)
// =============================================================
SolverOutput fem_solve(const SolverInput& in);

// =============================================================
// Stage Computation Functions
// =============================================================

// M4 + M6: Assemble global sparse stiffness matrix K
Eigen::SparseMatrix<double> assemble_stiffness(
    const SolverInput& in,
    const ElementRegistry& registry);

// M5: Assemble global force vector F
Eigen::VectorXd assemble_forces(
    const SolverInput& in,
    const ElementRegistry& registry);

// M7: Partition system by applying prescribed displacement BCs
ReducedSystem apply_boundary_conditions(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd&             F,
    const SolverInput&                 in);

// M8: Solve reduced system K_r * U_r = F_r
Eigen::VectorXd solve_linear_system(
    const ReducedSystem& reduced);

// Reconstruct full U vector from free and constrained components
Eigen::VectorXd reconstruct_full_displacement(
    const Eigen::VectorXd& Ur,
    const ReducedSystem&   reduced,
    int32_t                n_dofs_total);

// M9: Post-process strains, stresses, and von Mises invariants
void compute_element_results(
    const SolverInput&      in,
    const Eigen::VectorXd&  U,
    const ElementRegistry&  registry,
    SolverOutput&           out);

// Compute reaction forces R = K * U - F at constrained DOFs
Eigen::VectorXd compute_reactions(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd&             U,
    const Eigen::VectorXd&             F);

// =============================================================
// Context & Mapping Builders
// =============================================================

// Extract local ElementContext for element index e
ElementContext build_context(
    const SolverInput& in,
    int32_t            e);

// Extract local ElementDofMap for element index e
ElementDofMap build_element_dof_map(
    const SolverInput& in,
    int32_t            e);

} // namespace fem
```

---

## 7. Lifetime, Ownership & Memory Safety Rules

To guarantee memory safety and zero memory leaks across the C++ / Python boundary, the ownership and lifetime model is strictly defined:

| Data Structure | Owner Entity | Creation Context | Destruction Context |
| :--- | :--- | :--- | :--- |
| `RawModel` | Python GC | Created in `parse_inp()` (M1) | Destroyed after `build_model()` (M2) completes. |
| `FEMModel` | Python GC | Created in `build_model()` (M2) | Destroyed after `preprocess()` (M3) completes. |
| `SolverInput` | Python / C++ | Created in `preprocess()` (M3) | Passed by `const&` into C++; freed after `solve()` returns. |
| `ElementKernel` | C++ `std::unique_ptr` | Instantiated in M5/M6/M9 driver loops | Destroyed at the end of each element iteration. |
| `ElementContext` | C++ Stack | Built in driver loop | Discarded after `kernel->compute()` returns. |
| `ElementMatrices` | C++ Stack | Returned by `kernel->compute()` | Copied into triplets, discarded after element loop iteration. |
| `SolverOutput` | C++ / Python | Created inside `fem_solve()` | Copied into Python runtime, managed by Python GC. |

#### Memory Invariant Guarantee
No pointer or reference to C++ stack memory outlives the execution of `fem_solve()`. All crossing data is either value-copied or transferred via pybind11-managed NumPy buffers.

---

## 8. `femcore` Binding API (pybind11 Module Surface)

The C++ core is exposed to Python via a single `femcore` module. This is the authoritative API surface:

| Function | Signature | Returns / Behavior |
| :--- | :--- | :--- |
| `fem_solve` | `(SolverInput) -> SolverOutput` | Runs `fem_solve` (§6.2); maps over the structs below. |
| `element_meta` | `(name: str) -> dict` | `KeyError("[femcore] Unregistered element: <name>")` on unknown name; else dict: `vtk` (int VTK cell code), `npe`, `dim`, `dofs`, `strain_components`, `n_gauss`, `n_faces`, `faces` (= `face_nodes(1)` local connectivity), `vtk_permutation` (list). |
| `registered_elements` | `() -> list[str]` | **Sorted** list of registered ABAQUS names (registry itself uses `unordered_map`; the binding sorts to make output deterministic). |
| `registered_metadata` | `() -> dict[str, dict]` | `{ABAQUS name: element_meta(name)}` for every registered element, keyed in sorted order. |
| `test_element_matrix` | `(name, coords, E, nu, thickness, plane_mode=0) -> (ndof,ndof)` | Test hook: raw element stiffness for a single element. |
| `test_element_b_at` | `(name, coords, xi, eta, zeta=0, E=0, nu=0, thickness=1, plane_mode=0) -> ndarray` | Test hook: `B` at natural coordinates. |
| `test_element_pressure` | `(name, coords, face_id, p, thickness=1) -> (ndof,)` | Test hook: equivalent nodal pressure vector for one 1-based local face. |

**Class attributes** (`def_readwrite`, NumPy-backed):
* `SolverInput` — exactly the fields of doc 03 §3.1 (`coords`, `elem_type_name`, `elem_conn_flat`, `elem_conn_offsets`, `elem_mat`, `mat_E`, `mat_nu`, `mat_thickness`, `mat_plane_mode`, `dof_map`, `n_dofs_total`, `n_dofs_free`, `prescribed_dofs`, `prescribed_vals`, `point_load_dofs`, `point_load_vals`, `pressure_elem`, `pressure_face`, `pressure_val`, `n_nodes`, `n_elem`, `n_mat`, `dimension`).
* `SolverOutput` — `displacement (n,3)`, `reaction (n,3)`, `elem_stress (n_elem,6)`, `elem_strain (n_elem,6)`, `elem_von_mises (n_elem,)`, `elem_principal (n_elem,3)`, `converged (bool)`, `residual_norm (float)`. All arrays are contiguous `float64`; the engineering-shear Voigt order `[xx yy zz xy yz zx]` (doc 03 §5.1) applies to `elem_stress`/`elem_strain`.

This section closes audit item 15 (unspecified binding API) from doc 07.
