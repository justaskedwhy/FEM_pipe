# Package 2: Element Kernel Layer & Formulations Specification
**Document Version:** 1.0  
**Target Files:** `cpp/include/element_kernel.hpp`, `cpp/include/element_registry.hpp`, `cpp/src/elements_*.cpp`  
**Grounded Sources:** `docs/element_kernel.md`, `docs/formulation.md`, `docs/adding_an_element.md`

---

## 1. Element Kernel Architecture & Interface Specification

### 1.1 Purpose & Swappability Contract
The **Element Kernel Layer** provides a pluggable, polymorphic mechanism that decouples finite element mathematical formulations from the surrounding solver orchestration pipeline. It eliminates hardcoded branching (e.g., `switch(elem_type)` statements) across assembly, load evaluation, post-processing, and visualization.

#### Architectural Guarantee
Only **two components** in the entire codebase possess element-specific knowledge:
1. **`ElementKernel` Subclasses (C++)**: Contain shape functions, quadrature rules, strain-displacement matrices ($B$), constitutive matrices ($D$), face topologies, and VTK cell mappings.
2. **`ELEMENT_META` Cache (Python)**: A thin metadata dictionary queried from the C++ registry at import time.

All other solver modules (Assembler M6, Load Engine M5, BC Handler M7, Linear Solver M8, VTK Writer M10) operate strictly via element-agnostic interfaces.

---

### 1.2 Data Transfer Structs (`ElementContext` & `ElementMatrices`)

To achieve element swappability, the driver loop passes identical input data structures and receives uniform output data structures regardless of the underlying element geometry (2D vs. 3D, triangle vs. quad vs. hex).

#### `ElementContext` Struct
```cpp
namespace fem {

struct ElementContext {
    Eigen::Matrix<double, Eigen::Dynamic, 3> coords; // (npe x 3) node physical coordinates
    double E;             // Young's modulus (Pa)
    double nu;            // Poisson's ratio
    double thickness;     // Out-of-plane thickness (2D only; 1.0 for 3D)
    int    plane_mode;    // 0 = plane stress, 1 = plane strain, 2 = 3D
};

} // namespace fem
```

##### Rules & Invariants for `ElementContext`
* **3-Column Coordinate Matrix**: `coords` ALWAYS has 3 columns, even for 2D elements (where the 3rd column contains zeros). This guarantees uniform matrix operations across 2D and 3D kernels.
* **Element-Local Ordering**: `coords.row(i)` corresponds strictly to local node index `i` ($0 \le i < n_{pe}$).
* **3D Mode Defaults**: `thickness` and `plane_mode` are ignored by 3D kernels.

#### `ElementMatrices` Struct
```cpp
namespace fem {

struct ElementMatrices {
    Eigen::MatrixXd Ke;              // (ndof_e x ndof_e) dense element stiffness matrix
    Eigen::MatrixXd B;               // ((n_gauss * n_strain) x ndof_e) stacked strain-displacement matrix
    Eigen::MatrixXd D;               // (n_strain x n_strain) elastic constitutive matrix
    Eigen::MatrixXd gauss_points;    // (n_gauss x 3) natural coords; 2D stores z=0 (see note below)
    Eigen::VectorXd gauss_weights;   // (n_gauss) Gauss integration weights
    Eigen::VectorXd detJ;            // (n_gauss) Jacobian determinants at Gauss points
};

} // namespace fem
```

##### Rules & Invariants for `ElementMatrices`
* **Dense Stiffness Matrix ($K_e$)**: Element matrices are dense and symmetric positive semi-definite (positive definite once constrained).
* **Stacked $B$ Matrix Layout**: $B$ is stacked row-wise across Gauss points. Rows $[g \cdot n_{\text{strain}}, (g+1) \cdot n_{\text{strain}} - 1]$ store the strain-displacement matrix for Gauss point $g$. This allows post-processing and stress recovery to reuse evaluated $B$ matrices without recomputing shape function derivatives.
* **Single $D$ Storage**: For linear elasticity, $D$ is constant across all Gauss points within an element and is stored once.
* **Positive Jacobian Constraint**: $\det(\mathbf{J}_g) > 0$ for all Gauss points $g$. If $\det(\mathbf{J}_g) \le 0$, the element is distorted, inverted, or degenerate, and the kernel MUST throw a `std::runtime_error`. **Implementation note (shipped constant-strain kernels)**: the verbatim T3 code (§3, and T4 similarly) computes the signed area/determinant and uses `std::abs`, so a negatively wound element is silently sign-corrected and only a truly degenerate one (area < 1e-14) throws. Q4 and H8 honor the strict rule and throw on any $\det(\mathbf{J}_g) \le 0$. Both behaviors are accepted; the strict throw is the default for new quadratic/hex kernels.

---

### 1.3 The `ElementKernel` Abstract Interface

```cpp
namespace fem {

class ElementKernel {
public:
    virtual ~ElementKernel() = default;

    // ---------------------------------------------------------------
    // 1. Metadata Methods (Cheap, Constant Queries)
    // ---------------------------------------------------------------
    virtual std::string name()                 const; // ABAQUS canonical identifier (e.g., "CPS4"); see note below
    virtual int dim()                          const = 0; // Spatial dimension: 2 or 3
    virtual int nodes_per_element()            const = 0; // npe: 3, 4, or 8
    virtual int dofs_per_element()             const = 0; // ndof_e = npe * dim
    virtual int n_strain_components()          const = 0; // 3 (2D) or 6 (3D)
    virtual int n_gauss_points()               const = 0; // Number of integration points
    virtual int vtk_cell_type()                const = 0; // VTK XML cell type code

    // ---------------------------------------------------------------
    // 2. Core Stiffness & Integration Computation
    // ---------------------------------------------------------------
    virtual ElementMatrices compute(const ElementContext& ctx) const = 0;

    // ---------------------------------------------------------------
    // 3. Equivalent Pressure Force Vector Integration
    // ---------------------------------------------------------------
    virtual Eigen::VectorXd pressure_force(
        const ElementContext& ctx,
        int face_id,
        double pressure) const = 0;

    // ---------------------------------------------------------------
    // 4. Arbitrary Natural Coordinate Post-Processing Hook
    // ---------------------------------------------------------------
    virtual Eigen::MatrixXd B_at(
        const ElementContext& ctx,
        const Eigen::VectorXd& xi) const = 0;

    // ---------------------------------------------------------------
    // 5. Face Topology Queries
    // ---------------------------------------------------------------
    virtual int n_faces() const = 0;
    virtual std::vector<int> face_nodes(int face_id) const = 0; // 1-based face_id -> 0-based node vector

    // ---------------------------------------------------------------
    // 6. VTK Connectivity Permutation
    // ---------------------------------------------------------------
    virtual std::vector<int> vtk_permutation() const = 0;
};

} // namespace fem
```

**`name()` resolution**: `name()` is **non-pure** in the base class. It returns the ABAQUS canonical identifier set by the registry at factory time — `REGISTER_ELEMENT` constructs the kernel and injects the registration key via `set_registry_name()`. This is what lets one kernel class serve multiple ABAQUS names (e.g. `Q4Kernel` registered as both `CPS4` and `CPE4`); `name()` therefore reports the concrete registration key, not a class-level short name.

---

### 1.4 Static Registration & Macro System (`ElementRegistry`)

The element registry uses factory functions and C++ static initialization to allow new kernels to self-register without central source code modifications.

```cpp
// cpp/include/element_registry.hpp
#pragma once
#include "element_kernel.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <vector>

namespace fem {

class ElementRegistry {
public:
    using Factory = std::function<std::unique_ptr<ElementKernel>()>;

    static ElementRegistry& instance() {
        static ElementRegistry instance_;
        return instance_;
    }

    void register_kernel(const std::string& abaqus_name, Factory f) {
        factories_[abaqus_name] = f;
    }

    std::unique_ptr<ElementKernel> create(const std::string& abaqus_name) const {
        auto it = factories_.find(abaqus_name);
        if (it == factories_.end()) {
            throw std::runtime_error("[ElementRegistry] Unregistered element: " + abaqus_name);
        }
        return it->second();
    }

    bool has(const std::string& abaqus_name) const {
        return factories_.find(abaqus_name) != factories_.end();
    }

    std::vector<std::string> registered_names() const {
        std::vector<std::string> names;
        for (const auto& [name, _] : factories_) {
            names.push_back(name);
        }
        return names;
    }

private:
    std::unordered_map<std::string, Factory> factories_;
};

#define REGISTER_ELEMENT(ABAQUS_NAME, CLASS_NAME)                              \
    namespace {                                                                \
        const bool _registered_##CLASS_NAME = []() {                           \
            fem::ElementRegistry::instance().register_kernel(                  \
                ABAQUS_NAME,                                                   \
                []() { return std::make_unique<CLASS_NAME>(); });              \
            return true;                                                       \
        }();                                                                   \
    }

} // namespace fem
```

#### Registration Macro Usage
At the end of any element implementation file (e.g., `cpp/src/elements_t3.cpp`), add:
```cpp
REGISTER_ELEMENT("CPS3", fem::T3Kernel);
REGISTER_ELEMENT("CPE3", fem::T3Kernel);
```
This guarantees self-registration during executable initialization without touching any central factory files.

---

## 2. Finite Element Formulations & Mathematical Derivations

### 2.1 Governing Continuum Equations & Weak Form

#### Strong Form (Linear Elastostatics)
For a domain $\Omega \subset \mathbb{R}^d$ ($d \in \{2, 3\}$) bounded by $\Gamma = \Gamma_u \cup \Gamma_t$:
$$\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \mathbf{0} \quad \text{in } \Omega$$
$$\boldsymbol{\sigma} = \mathbf{D} : \boldsymbol{\varepsilon}, \qquad \boldsymbol{\varepsilon} = \tfrac{1}{2}\left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T\right)$$
$$\mathbf{u} = \bar{\mathbf{u}} \quad \text{on } \Gamma_u, \qquad \boldsymbol{\sigma} \cdot \mathbf{n} = \bar{\mathbf{t}} \quad \text{on } \Gamma_t$$

#### Weak Form (Principle of Virtual Work)
Find displacement field $\mathbf{u} \in \mathcal{V}$ such that for all test functions $\mathbf{v} \in \mathcal{V}_0$:
$$\int_{\Omega} \boldsymbol{\varepsilon}(\mathbf{v})^T \mathbf{D} \, \boldsymbol{\varepsilon}(\mathbf{u}) \, d\Omega = \int_{\Omega} \mathbf{v}^T \mathbf{b} \, d\Omega + \int_{\Gamma_t} \mathbf{v}^T \bar{\mathbf{t}} \, d\Gamma$$

#### Discretized Equilibrium Equations
Introducing shape function approximations $\mathbf{u}(\mathbf{x}) = \sum_{a=1}^{n_{pe}} N_a(\mathbf{x}) \mathbf{u}_a$:
$$\mathbf{K} \mathbf{U} = \mathbf{F}$$
$$\mathbf{K} = \bigwedge_e \mathbf{K}_e = \bigwedge_e \int_{\Omega_e} \mathbf{B}^T \mathbf{D} \mathbf{B} \, d\Omega$$
$$\mathbf{F} = \bigwedge_e \left( \int_{\Omega_e} \mathbf{N}^T \mathbf{b} \, d\Omega + \int_{\Gamma_e} \mathbf{N}^T \bar{\mathbf{t}} \, d\Gamma \right)$$

---

### 2.2 Isoparametric Kinematics & Mapping

Each physical element $\Omega_e$ is mapped from a standardized reference domain $\hat{\Omega}$ parameterized by natural coordinates $\boldsymbol{\xi} = (\xi, \eta, \zeta)$:

$$\mathbf{x}(\boldsymbol{\xi}) = \sum_{a=1}^{n_{pe}} N_a(\boldsymbol{\xi}) \, \mathbf{x}_a, \qquad \mathbf{u}(\boldsymbol{\xi}) = \sum_{a=1}^{n_{pe}} N_a(\boldsymbol{\xi}) \, \mathbf{u}_a$$

#### Jacobian Matrix
$$\mathbf{J}(\boldsymbol{\xi}) = \frac{\partial \mathbf{x}}{\partial \boldsymbol{\xi}} = \sum_{a=1}^{n_{pe}} \mathbf{x}_a \otimes \frac{\partial N_a}{\partial \boldsymbol{\xi}}$$

In 2D matrix notation ($2 \times 2$):
$$\mathbf{J} = \begin{bmatrix} \frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\ \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta} \end{bmatrix} = \begin{bmatrix} \sum_a \frac{\partial N_a}{\partial \xi} x_a & \sum_a \frac{\partial N_a}{\partial \xi} y_a \\ \sum_a \frac{\partial N_a}{\partial \eta} x_a & \sum_a \frac{\partial N_a}{\partial \eta} y_a \end{bmatrix}$$

#### Spatial Derivatives & Differential Volume
$$\begin{bmatrix} \frac{\partial N_a}{\partial x} \\ \frac{\partial N_a}{\partial y} \end{bmatrix} = \mathbf{J}^{-T} \begin{bmatrix} \frac{\partial N_a}{\partial \xi} \\ \frac{\partial N_a}{\partial \eta} \end{bmatrix}, \qquad d\Omega = \det(\mathbf{J}) \, d\hat{\Omega}$$

---

### 2.3 Constitutive Tensors (Voigt Matrix $D$)

#### 1. 2D Plane Stress ($\sigma_{zz} = 0$, $\varepsilon_{zz} = -\frac{\nu}{1-\nu}(\varepsilon_{xx} + \varepsilon_{yy})$)
$$\mathbf{D}_{\text{stress}} = \frac{E}{1-\nu^2} \begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{bmatrix}$$

#### 2. 2D Plane Strain ($\varepsilon_{zz} = 0$, $\sigma_{zz} = \nu(\sigma_{xx} + \sigma_{yy})$)
$$\mathbf{D}_{\text{strain}} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix} 1-\nu & \nu & 0 \\ \nu & 1-\nu & 0 \\ 0 & 0 & \frac{1-2\nu}{2} \end{bmatrix}$$

#### 3. 3D Isotropic Linear Elasticity
$$\mathbf{D}_{\text{3D}} = \frac{E}{(1+\nu)(1-2\nu)} \begin{bmatrix} 
1-\nu & \nu & \nu & 0 & 0 & 0 \\ 
\nu & 1-\nu & \nu & 0 & 0 & 0 \\ 
\nu & \nu & 1-\nu & 0 & 0 & 0 \\ 
0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\ 
0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\ 
0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}$$

---

### 2.4 Element 1: T3 (Linear Triangle - 2D)

#### Reference Domain & Shape Functions
Unit reference triangle $0 \le \xi, \eta$ and $\xi + \eta \le 1$. Shape functions match area coordinates $L_1 = 1 - \xi - \eta$, $L_2 = \xi$, $L_3 = \eta$:
$$N_1 = 1 - \xi - \eta, \quad N_2 = \xi, \quad N_3 = \eta$$

#### Natural Derivatives (Constant)
$$\frac{\partial \mathbf{N}}{\partial \boldsymbol{\xi}} = \begin{bmatrix} -1 & 1 & 0 \\ -1 & 0 & 1 \end{bmatrix}$$

#### Jacobian & Area Evaluation
$$\mathbf{J} = \begin{bmatrix} x_2 - x_1 & y_2 - y_1 \\ x_3 - x_1 & y_3 - y_1 \end{bmatrix}, \qquad \det(\mathbf{J}) = 2A = (x_2-x_1)(y_3-y_1) - (x_3-x_1)(y_2-y_1)$$

#### Cartesian Derivatives
$$\frac{\partial N_1}{\partial x} = \frac{y_2 - y_3}{2A}, \quad \frac{\partial N_1}{\partial y} = \frac{x_3 - x_2}{2A}$$
$$\frac{\partial N_2}{\partial x} = \frac{y_3 - y_1}{2A}, \quad \frac{\partial N_2}{\partial y} = \frac{x_1 - x_3}{2A}$$
$$\frac{\partial N_3}{\partial x} = \frac{y_1 - y_2}{2A}, \quad \frac{\partial N_3}{\partial y} = \frac{x_2 - x_1}{2A}$$

#### Strain-Displacement Matrix $\mathbf{B}$ ($3 \times 6$)
$$\mathbf{B} = \begin{bmatrix}
N_{1,x} & 0 & N_{2,x} & 0 & N_{3,x} & 0 \\
0 & N_{1,y} & 0 & N_{2,y} & 0 & N_{3,y} \\
N_{1,y} & N_{1,x} & N_{2,y} & N_{2,x} & N_{3,y} & N_{3,x}
\end{bmatrix}$$

#### Stiffness Matrix Evaluation
Since $\mathbf{B}$ is constant throughout the element, integration reduces to scaling by physical area $A$ and thickness $t$:
$$\mathbf{K}_e = A \cdot t \cdot \mathbf{B}^T \mathbf{D} \mathbf{B}$$

---

### 2.5 Element 2: Q4 (Bilinear Quad - 2D)

#### Reference Domain & Shape Functions
Square reference domain $(\xi, \eta) \in [-1, 1]^2$. Nodes at $(\xi_a, \eta_a) \in \{(-1,-1), (1,-1), (1,1), (-1,1)\}$:
$$N_a(\xi, \eta) = \tfrac{1}{4}(1 + \xi_a \xi)(1 + \eta_a \eta)$$

Expanded:
$$N_1 = \tfrac{1}{4}(1-\xi)(1-\eta), \quad N_2 = \tfrac{1}{4}(1+\xi)(1-\eta)$$
$$N_3 = \tfrac{1}{4}(1+\xi)(1+\eta), \quad N_4 = \tfrac{1}{4}(1-\xi)(1+\eta)$$

#### Natural Derivatives
$$\frac{\partial N_a}{\partial \xi} = \tfrac{1}{4} \xi_a (1 + \eta_a \eta), \qquad \frac{\partial N_a}{\partial \eta} = \tfrac{1}{4} \eta_a (1 + \xi_a \xi)$$

#### Numerical Integration ($2 \times 2$ Gauss Rule)
Integration points at $(\pm 1/\sqrt{3}, \pm 1/\sqrt{3})$ with weights $w_g = 1.0$:
$$\mathbf{K}_e = t \sum_{g=1}^{4} \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \det(\mathbf{J}_g) w_g$$

---

### 2.6 Element 3: T4 (Linear Tetrahedron - 3D)

#### Reference Domain & Shape Functions
Unit reference tetrahedron $\xi, \eta, \zeta \ge 0$, $\xi + \eta + \zeta \le 1$. Volume coordinates:
$$N_1 = 1 - \xi - \eta - \zeta, \quad N_2 = \xi, \quad N_3 = \eta, \quad N_4 = \zeta$$

#### Jacobian Matrix (Constant)
$$\mathbf{J} = \begin{bmatrix} 
x_2 - x_1 & y_2 - y_1 & z_2 - z_1 \\
x_3 - x_1 & y_3 - y_1 & z_3 - z_1 \\
x_4 - x_1 & y_4 - y_1 & z_4 - z_1 
\end{bmatrix}, \qquad \det(\mathbf{J}) = 6V$$

#### Strain-Displacement Matrix $\mathbf{B}$ ($6 \times 12$)
$$\mathbf{B}_a = \begin{bmatrix}
N_{a,x} & 0 & 0 \\
0 & N_{a,y} & 0 \\
0 & 0 & N_{a,z} \\
N_{a,y} & N_{a,x} & 0 \\
0 & N_{a,z} & N_{a,y} \\
N_{a,z} & 0 & N_{a,x}
\end{bmatrix} \quad (a = 1, \dots, 4)$$

#### Stiffness Matrix Evaluation
Single Gauss point at tetrahedron centroid $(\frac{1}{4}, \frac{1}{4}, \frac{1}{4})$ with weight $V = \frac{\det(\mathbf{J})}{6}$:
$$\mathbf{K}_e = V \cdot \mathbf{B}^T \mathbf{D} \mathbf{B}$$

---

### 2.7 Element 4: H8 (Trilinear Hexahedron - 3D)

#### Reference Domain & Shape Functions
Reference cube $(\xi, \eta, \zeta) \in [-1, 1]^3$. Corner nodes $(\xi_a, \eta_a, \zeta_a) \in \{-1, +1\}^3$:
$$N_a(\xi, \eta, \zeta) = \tfrac{1}{8}(1 + \xi_a \xi)(1 + \eta_a \eta)(1 + \zeta_a \zeta)$$

#### Natural Derivatives
$$\frac{\partial N_a}{\partial \xi} = \tfrac{1}{8} \xi_a (1 + \eta_a \eta)(1 + \zeta_a \zeta)$$
$$\frac{\partial N_a}{\partial \eta} = \tfrac{1}{8} \eta_a (1 + \xi_a \xi)(1 + \zeta_a \zeta)$$
$$\frac{\partial N_a}{\partial \zeta} = \tfrac{1}{8} \zeta_a (1 + \xi_a \xi)(1 + \eta_a \eta)$$

#### Numerical Integration ($2 \times 2 \times 2$ Gauss Rule)
Eight Gauss integration points at $(\pm 1/\sqrt{3}, \pm 1/\sqrt{3}, \pm 1/\sqrt{3})$ with unit weights:
$$\mathbf{K}_e = \sum_{g=1}^{8} \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \det(\mathbf{J}_g)$$

---

### 2.8 Gauss Quadrature Rules & Integration Tables

| Dimension | Rule / Topology | Point $g$ | Coordinates $(\xi, \eta, \zeta)$ | Weight $w_g$ |
| :--- | :--- | :--- | :--- | :--- |
| **1D** | 2-point Gauss | 1<br>2 | $-\frac{1}{\sqrt{3}}$<br>$+\frac{1}{\sqrt{3}}$ | 1.0<br>1.0 |
| **2D Triangle** | 1-point centroid | 1 | $(\frac{1}{3}, \frac{1}{3}, 0)$ | 1.0 (absorbed in $A$) |
| **2D Quad** | $2 \times 2$ tensor product | 1<br>2<br>3<br>4 | $(-\\frac{1}{\\sqrt{3}}, -\\frac{1}{\\sqrt{3}})$<br>$(+\\frac{1}{\\sqrt{3}}, -\\frac{1}{\\sqrt{3}})$<br>$(+\\frac{1}{\\sqrt{3}}, +\\frac{1}{\\sqrt{3}})$<br>$(-\\frac{1}{\\sqrt{3}}, +\\frac{1}{\\sqrt{3}})$ | 1.0<br>1.0<br>1.0<br>1.0 |
| **3D Tet** | 1-point centroid | 1 | $(\frac{1}{4}, \frac{1}{4}, \frac{1}{4})$ | $\frac{1}{6}$ (absorbed in $V$) |
| **3D Hex** | $2 \times 2 \times 2$ tensor product | $1 \dots 8$ | $(\pm \frac{1}{\sqrt{3}}, \pm \frac{1}{\sqrt{3}}, \pm \frac{1}{\sqrt{3}})$ | 1.0 each |

---

### 2.9 Consistent Surface Pressure Load Integration

For element boundary face $\Gamma_e$ under normal pressure $p$ (where positive $p$ acts compressive **into** the element):
$$\mathbf{t} = -p \, \mathbf{n}$$
$$\mathbf{F}_e^{(p)} = \int_{\Gamma_e} \mathbf{N}^T \mathbf{t} \, d\Gamma$$

#### T3 Straight Edge Integration Example
For face edge of length $L = \sqrt{\Delta x^2 + \Delta y^2}$ with outward normal $\mathbf{n} = (\frac{\Delta y}{L}, -\frac{\Delta x}{L})$:
$$\mathbf{F}_e^{(p)} = \frac{p L}{2} \begin{bmatrix} -n_x \\ -n_y \\ -n_x \\ -n_y \end{bmatrix} \quad \text{distributed equally to the two edge nodes.}$$

#### Shipped Kernel Implementations
* **T3 / Q4** (2D, face = edge): edge normal $\mathbf{n} = (\frac{\Delta y}{L}, -\frac{\Delta x}{L})$ from the CCW-wound edge, traction $\mathbf{t} = -p\,\mathbf{n}$, equivalent nodal forces $f_{\text{node}} = \frac{p L}{2} \cdot t$ ($t$ = thickness) split equally to the two edge-end nodes. Q4 faces are pairs `{0,1},{1,2},{2,3},{3,0}`.
* **T4** (3D, const-pressure face): face area $A$ from the cross-product magnitude, outward normal re-oriented away from the element centroid (4-node average), equal split $A\,p/3$ to the three face nodes. Local faces: `{0,1,2},{0,1,3},{1,2,3},{0,2,3}`.
* **H8** (3D, bilinear face): $2\times2$ Gauss quadrature over the quad face, $\mathbf{n} = \mathbf{r}_\xi \times \mathbf{r}_\eta$ re-oriented outward from the element centroid, $\mathbf{F}_e = \sum_g N_g^T \mathbf{t}\,\det(\mathbf{J}_g)$. Local faces (ABAQUS P1–P6 order): `{0,1,2,3},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}`.
* Positive pressure $p > 0$ always produces compressive (inward) equivalent forces; degenerate faces (zero length/area) throw.

---

### 2.10 Stress & Strain Recovery (Voigt, von Mises, Principal Stresses)

Given element nodal displacements $\mathbf{U}_e$:
$$\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{U}_e, \qquad \boldsymbol{\sigma} = \mathbf{D} \boldsymbol{\varepsilon}$$

**Evaluation point (all shipped kernels)**: element stress/strain are reported as the **Gauss-point arithmetic mean** — the post-processor averages $\boldsymbol{\varepsilon}_g = \mathbf{B}_g \mathbf{U}_e$ over the stacked $\mathbf{B}$ rows and applies $\boldsymbol{\sigma} = \mathbf{D} \bar{\boldsymbol{\varepsilon}}$. For constant-strain elements (T3, T4) this is exact at every point; for Q4/H8 it is the centroidal-equivalent (mean) recovery used by `compute_element_results`. `B_at()` exists for arbitrary natural-coordinate evaluation (e.g. basis verification).

#### Voigt Array Expansion (2D to 3D for VTK Output)
* **Plane Stress**: $\sigma_{zz} = 0, \quad \varepsilon_{zz} = -\frac{\nu}{1-\nu}(\varepsilon_{xx} + \varepsilon_{yy})$
* **Plane Strain**: $\varepsilon_{zz} = 0, \quad \sigma_{zz} = \nu(\sigma_{xx} + \sigma_{yy})$

#### von Mises Equivalent Stress ($\sigma_{\text{vm}}$)
$$\sigma_{\text{vm}} = \sqrt{\tfrac{1}{2}\left[ (\sigma_{xx}-\sigma_{yy})^2 + (\sigma_{yy}-\sigma_{zz})^2 + (\sigma_{zz}-\sigma_{xx})^2 \right] + 3(\sigma_{xy}^2 + \sigma_{yz}^2 + \sigma_{zx}^2)}$$

#### Principal Stresses ($\sigma_1 \ge \sigma_2 \ge \sigma_3$)
Computed as the eigenvalues of the $3 \times 3$ symmetric stress tensor $\boldsymbol{\sigma}$.

---

## 3. Complete C++ Implementation: T3 Kernel

**File Path:** `cpp/src/elements_t3.cpp`

```cpp
#include "element_kernel.hpp"
#include "element_registry.hpp"
#include <stdexcept>
#include <cmath>

namespace fem {

class T3Kernel : public ElementKernel {
public:
    // ---------------------------------------------------------------
    // Metadata
    // ---------------------------------------------------------------
    int dim()                        const override { return 2; }
    int nodes_per_element()          const override { return 3; }
    int dofs_per_element()           const override { return 6; }
    int n_strain_components()        const override { return 3; }
    int n_gauss_points()             const override { return 1; }
    int vtk_cell_type()              const override { return 5; /* VTK_TRIANGLE */ }

    // ---------------------------------------------------------------
    // Core Stiffness Computation
    // ---------------------------------------------------------------
    ElementMatrices compute(const ElementContext& ctx) const override {
        const double x1 = ctx.coords(0,0), y1 = ctx.coords(0,1);
        const double x2 = ctx.coords(1,0), y2 = ctx.coords(1,1);
        const double x3 = ctx.coords(2,0), y3 = ctx.coords(2,1);

        // Signed area 2A = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
        const double twoA = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
        const double A = 0.5 * std::abs(twoA);
        if (A < 1e-14) {
            throw std::runtime_error("T3Kernel: Degenerate or inverted triangle element.");
        }

        const double sign = (twoA > 0.0) ? 1.0 : -1.0;
        const double dNdx[3] = {
            sign * (y2 - y3) / (2.0 * A),
            sign * (y3 - y1) / (2.0 * A),
            sign * (y1 - y2) / (2.0 * A)
        };
        const double dNdy[3] = {
            sign * (x3 - x2) / (2.0 * A),
            sign * (x1 - x3) / (2.0 * A),
            sign * (x2 - x1) / (2.0 * A)
        };

        // Construct B matrix (3 x 6)
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 6);
        for (int a = 0; a < 3; ++a) {
            B(0, 2*a + 0) = dNdx[a];
            B(1, 2*a + 1) = dNdy[a];
            B(2, 2*a + 0) = dNdy[a];
            B(2, 2*a + 1) = dNdx[a];
        }

        // Construct D matrix (3 x 3)
        Eigen::MatrixXd D = build_D_2d(ctx.E, ctx.nu, ctx.plane_mode);

        // Ke = A * thickness * B^T * D * B
        const double factor = A * ctx.thickness;
        Eigen::MatrixXd Ke = factor * B.transpose() * D * B;

        ElementMatrices out;
        out.Ke = Ke;
        out.B = B;
        out.D = D;
        out.gauss_points = Eigen::MatrixXd::Zero(1, 3);
        out.gauss_points(0,0) = 1.0 / 3.0;
        out.gauss_points(0,1) = 1.0 / 3.0;
        out.gauss_weights = Eigen::VectorXd::Constant(1, 1.0);
        out.detJ = Eigen::VectorXd::Constant(1, 2.0 * A);
        return out;
    }

    // ---------------------------------------------------------------
    // Equivalent Pressure Vector Integration
    // ---------------------------------------------------------------
    Eigen::VectorXd pressure_force(
        const ElementContext& ctx, int face_id, double p) const override 
    {
        if (face_id < 1 || face_id > 3) {
            throw std::runtime_error("T3Kernel: face_id out of bounds (must be 1, 2, or 3)");
        }

        static const int edge[3][2] = {{0,1}, {1,2}, {2,0}};
        const int n0 = edge[face_id - 1][0];
        const int n1 = edge[face_id - 1][1];

        const double x0 = ctx.coords(n0, 0), y0 = ctx.coords(n0, 1);
        const double x1 = ctx.coords(n1, 0), y1 = ctx.coords(n1, 1);
        const double dx = x1 - x0, dy = y1 - y0;
        const double L  = std::sqrt(dx*dx + dy*dy);
        if (L < 1e-14) {
            throw std::runtime_error("T3Kernel: Edge length is zero.");
        }

        // Outward normal vector
        const double nx =  dy / L;
        const double ny = -dx / L;

        // Traction vector t = -p * n
        const double tx = -p * nx;
        const double ty = -p * ny;

        Eigen::VectorXd fe = Eigen::VectorXd::Zero(6);
        const double f_node_x = tx * L / 2.0 * ctx.thickness;
        const double f_node_y = ty * L / 2.0 * ctx.thickness;

        fe(2*n0 + 0) = f_node_x;
        fe(2*n0 + 1) = f_node_y;
        fe(2*n1 + 0) = f_node_x;
        fe(2*n1 + 1) = f_node_y;
        return fe;
    }

    // ---------------------------------------------------------------
    // Arbitrary Natural Coordinate B Matrix (Constant for T3)
    // ---------------------------------------------------------------
    Eigen::MatrixXd B_at(
        const ElementContext& ctx,
        const Eigen::VectorXd& /*xi*/) const override 
    {
        return compute(ctx).B;
    }

    // ---------------------------------------------------------------
    // Topology & VTK
    // ---------------------------------------------------------------
    int n_faces() const override { return 3; }

    std::vector<int> face_nodes(int face_id) const override {
        static const std::vector<std::vector<int>> f = {{0,1}, {1,2}, {2,0}};
        if (face_id < 1 || face_id > 3) {
            throw std::runtime_error("T3Kernel: Invalid face_id");
        }
        return f[face_id - 1];
    }

    std::vector<int> vtk_permutation() const override { return {0, 1, 2}; }

private:
    static Eigen::MatrixXd build_D_2d(double E, double nu, int plane_mode) {
        Eigen::MatrixXd D(3, 3);
        if (plane_mode == 0) { // Plane stress
            const double c = E / (1.0 - nu*nu);
            D << 1.0, nu,  0.0,
                 nu,  1.0, 0.0,
                 0.0, 0.0, (1.0 - nu) / 2.0;
            D *= c;
        } else if (plane_mode == 1) { // Plane strain
            const double c = E / ((1.0 + nu) * (1.0 - 2.0*nu));
            D << 1.0 - nu, nu,        0.0,
                 nu,       1.0 - nu,  0.0,
                 0.0,      0.0,       (1.0 - 2.0*nu) / 2.0;
            D *= c;
        } else {
            throw std::runtime_error("T3Kernel: Invalid plane_mode (must be 0 or 1)");
        }
        return D;
    }
};

} // namespace fem

// Static self-registration
REGISTER_ELEMENT("CPS3", fem::T3Kernel);
REGISTER_ELEMENT("CPE3", fem::T3Kernel);
```

---

## 4. Skeleton Implementation: Q4 Kernel

**File Path:** `cpp/src/elements_q4.cpp`

```cpp
#include "element_kernel.hpp"
#include "element_registry.hpp"
#include <stdexcept>
#include <cmath>

namespace fem {

class Q4Kernel : public ElementKernel {
public:
    int dim()                        const override { return 2; }
    int nodes_per_element()          const override { return 4; }
    int dofs_per_element()           const override { return 8; }
    int n_strain_components()        const override { return 3; }
    int n_gauss_points()             const override { return 4; }
    int vtk_cell_type()              const override { return 9; /* VTK_QUAD */ }

    ElementMatrices compute(const ElementContext& ctx) const override {
        static const double gp = 1.0 / std::sqrt(3.0);
        static const double gauss[4][2] = {{-gp,-gp}, {gp,-gp}, {gp,gp}, {-gp,gp}};
        static const double node_ref[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};

        Eigen::MatrixXd D = build_D_2d(ctx.E, ctx.nu, ctx.plane_mode);
        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(8, 8);
        Eigen::MatrixXd B_stacked(12, 8);
        Eigen::VectorXd detJ_vec(4);
        Eigen::MatrixXd gp_coords(4, 3);
        gp_coords.setZero();

        for (int g = 0; g < 4; ++g) {
            const double xi = gauss[g][0], eta = gauss[g][1];
            gp_coords(g, 0) = xi; gp_coords(g, 1) = eta;

            // Shape function natural derivatives (4 x 2)
            Eigen::Matrix<double, 4, 2> dN_dxi;
            for (int a = 0; a < 4; ++a) {
                dN_dxi(a, 0) = 0.25 * node_ref[a][0] * (1.0 + node_ref[a][1] * eta);
                dN_dxi(a, 1) = 0.25 * node_ref[a][1] * (1.0 + node_ref[a][0] * xi);
            }

            // Jacobian matrix J = dN_dxi^T * coords_2d
            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            for (int a = 0; a < 4; ++a) {
                J += dN_dxi.row(a).transpose() * ctx.coords.row(a).head<2>();
            }

            const double detJ = J.determinant();
            if (detJ <= 0.0) {
                throw std::runtime_error("Q4Kernel: Degenerate Jacobian at Gauss point.");
            }
            detJ_vec[g] = detJ;

            // Spatial derivatives dN_dx = dN_dxi * J^-1
            Eigen::Matrix2d Jinv = J.inverse();
            Eigen::Matrix<double, 4, 2> dN_dx = dN_dxi * Jinv.transpose();

            // Construct B matrix (3 x 8)
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
            for (int a = 0; a < 4; ++a) {
                B(0, 2*a + 0) = dN_dx(a, 0);
                B(1, 2*a + 1) = dN_dx(a, 1);
                B(2, 2*a + 0) = dN_dx(a, 1);
                B(2, 2*a + 1) = dN_dx(a, 0);
            }

            Ke += B.transpose() * D * B * detJ * ctx.thickness;
            B_stacked.block<3, 8>(3*g, 0) = B;
        }

        ElementMatrices out;
        out.Ke = Ke;
        out.B = B_stacked;
        out.D = D;
        out.gauss_points = gp_coords;
        out.gauss_weights = Eigen::VectorXd::Constant(4, 1.0);
        out.detJ = detJ_vec;
        return out;
    }

    Eigen::VectorXd pressure_force(
        const ElementContext& ctx, int face_id, double p) const override 
    {
        if (face_id < 1 || face_id > n_faces()) {
            throw std::runtime_error("Q4Kernel: face_id out of bounds");
        }
        static const int edge[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};
        const int n0 = edge[face_id - 1][0];
        const int n1 = edge[face_id - 1][1];

        const double x0 = ctx.coords(n0, 0), y0 = ctx.coords(n0, 1);
        const double x1 = ctx.coords(n1, 0), y1 = ctx.coords(n1, 1);
        const double dx = x1 - x0, dy = y1 - y0;
        const double L  = std::sqrt(dx*dx + dy*dy);
        if (L < 1e-14) {
            throw std::runtime_error("Q4Kernel: Edge length is zero.");
        }

        // Outward normal vector for a CCW-wound edge
        const double nx =  dy / L;
        const double ny = -dx / L;

        // Traction vector t = -p * n
        const double tx = -p * nx;
        const double ty = -p * ny;

        Eigen::VectorXd fe = Eigen::VectorXd::Zero(8);
        const double f_node_x = tx * L / 2.0 * ctx.thickness;
        const double f_node_y = ty * L / 2.0 * ctx.thickness;

        fe(2*n0 + 0) = f_node_x;
        fe(2*n0 + 1) = f_node_y;
        fe(2*n1 + 0) = f_node_x;
        fe(2*n1 + 1) = f_node_y;
        return fe;
    }

    Eigen::MatrixXd B_at(
        const ElementContext& ctx,
        const Eigen::VectorXd& xi) const override 
    {
        static const double node_ref[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};
        const double xi_v = xi(0), eta = xi(1);

        Eigen::Matrix<double, 4, 2> dN_dxi;
        for (int a = 0; a < 4; ++a) {
            dN_dxi(a, 0) = 0.25 * node_ref[a][0] * (1.0 + node_ref[a][1] * eta);
            dN_dxi(a, 1) = 0.25 * node_ref[a][1] * (1.0 + node_ref[a][0] * xi_v);
        }

        Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
        for (int a = 0; a < 4; ++a) {
            J += dN_dxi.row(a).transpose() * ctx.coords.row(a).head<2>();
        }

        const double detJ = J.determinant();
        if (detJ <= 0.0) {
            throw std::runtime_error("Q4Kernel: Degenerate Jacobian at evaluation point.");
        }

        Eigen::Matrix2d Jinv = J.inverse();
        Eigen::Matrix<double, 4, 2> dN_dx = dN_dxi * Jinv.transpose();

        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
        for (int a = 0; a < 4; ++a) {
            B(0, 2*a + 0) = dN_dx(a, 0);
            B(1, 2*a + 1) = dN_dx(a, 1);
            B(2, 2*a + 0) = dN_dx(a, 1);
            B(2, 2*a + 1) = dN_dx(a, 0);
        }
        return B;
    }

    int n_faces() const override { return 4; }

    std::vector<int> face_nodes(int face_id) const override {
        static const std::vector<std::vector<int>> f = {{0,1}, {1,2}, {2,3}, {3,0}};
        return f.at(face_id - 1);
    }

    std::vector<int> vtk_permutation() const override { return {0, 1, 2, 3}; }

private:
    static Eigen::MatrixXd build_D_2d(double E, double nu, int plane_mode) {
        Eigen::MatrixXd D(3, 3);
        if (plane_mode == 0) {
            const double c = E / (1.0 - nu*nu);
            D << 1.0, nu, 0.0, nu, 1.0, 0.0, 0.0, 0.0, (1.0-nu)/2.0;
            D *= c;
        } else if (plane_mode == 1) {
            const double c = E / ((1.0+nu)*(1.0-2.0*nu));
            D << 1.0-nu, nu, 0.0, nu, 1.0-nu, 0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0;
            D *= c;
        } else {
            throw std::runtime_error("Q4Kernel: Invalid plane_mode (must be 0 or 1)");
        }
        return D;
    }
};

} // namespace fem

REGISTER_ELEMENT("CPS4", fem::Q4Kernel);
REGISTER_ELEMENT("CPE4", fem::Q4Kernel);
```

---

## 5. Practical Guide: How to Add a New Element

### 5.1 Three-Step Workflow

1. **Write `cpp/src/elements_<name>.cpp`**: Implement `ElementKernel` and add `REGISTER_ELEMENT` calls at the bottom.
2. **Add Source to `CMakeLists.txt`**: Add `cpp/src/elements_<name>.cpp` to `pybind11_add_module(femcore ...)`.
3. **Verify with Benchmark Tests**: Run single-element patch tests and the mandatory Q4/T3 element swap test.

---

### 5.2 Method-by-Method Implementation Guidelines

#### 1. `compute(const ElementContext& ctx)`
* Loop over all Gauss points $g = 1 \dots n_{\text{gauss}}$.
* Evaluate shape functions $N_a(\boldsymbol{\xi}_g)$ and natural derivatives $\frac{\partial N_a}{\partial \boldsymbol{\xi}}$.
* Compute physical Jacobian $\mathbf{J} = \sum_a \mathbf{x}_a \otimes \frac{\partial N_a}{\partial \boldsymbol{\xi}}$.
* Check $\det(\mathbf{J}) > 0$. Throw `std::runtime_error` if non-positive.
* Compute spatial derivatives $\nabla N_a = \mathbf{J}^{-T} \hat{\nabla} N_a$.
* Construct strain-displacement matrix $\mathbf{B}_g$.
* Accumulate stiffness: $\mathbf{K}_e += \mathbf{B}_g^T \mathbf{D} \mathbf{B}_g \det(\mathbf{J}_g) w_g \cdot (\text{thickness})$.

#### 2. `pressure_force(const ElementContext& ctx, int face_id, double p)`
* Look up local face nodes via `face_nodes(face_id)`.
* Compute face outward normal $\mathbf{n}$.
* Integrate surface traction $\mathbf{t} = -p \mathbf{n}$ across the face using 1D or 2D Gauss integration.
* Distribute loads to corresponding nodal DOFs.

#### 3. `B_at(const ElementContext& ctx, const Eigen::VectorXd& xi)`
* For constant-strain elements (T3, T4), return the constant $\mathbf{B}$ matrix.
* For higher-order elements (Q4, H8), evaluate $\mathbf{B}(\boldsymbol{\xi})$ at the specified natural coordinate.

---

### 5.3 Verification & Quality Control Checklist

| Test Suite | Verification Procedure | Acceptance Criteria |
| :--- | :--- | :--- |
| **Patch Test** | Apply linear displacement BCs $u = \mathbf{A}\mathbf{x} + \mathbf{b}$ to boundary nodes | Interior nodes recover exact displacements; element strain $\boldsymbol{\varepsilon}$ is constant to $< 10^{-10}$ |
| **Rigid Body Translation** | Apply uniform translation vector $\mathbf{u}_{\text{rigid}}$ to all nodes | $\| \mathbf{K}_e \mathbf{u}_{\text{rigid}} \| < 10^{-12}$ |
| **Rigid Body Rotation** | Apply small-angle rigid rotation to all nodes | $\| \mathbf{K}_e \mathbf{u}_{\text{rot}} \| < 10^{-12}$ |
| **Symmetry Test** | Check matrix symmetry of computed stiffness matrix | $\frac{\| \mathbf{K}_e - \mathbf{K}_e^T \|}{\| \mathbf{K}_e \|} < 10^{-12}$ |
| **Positive Semi-Definiteness** | Compute eigenvalues of unconstrained $\mathbf{K}_e$ | All eigenvalues $\lambda_i \ge -10^{-10}$ |
| **Pressure Resultant** | Apply uniform pressure load $p$ to face $\Gamma_f$ | $\sum \mathbf{F}_{\text{nodal}} = p \cdot \text{Area}(\Gamma_f) \cdot (-\mathbf{n})$ |
| **Swap Test** | Run cantilever benchmark with new element vs. baseline element | Displacement converges monotonically to analytical solution under refinement |

---

### 5.4 Common Pitfalls & Diagnostic Table

| Symptom | Primary Cause | Fix / Resolution |
| :--- | :--- | :--- |
| **Silent numerical corruption** | Unchecked $\det(\mathbf{J}) \le 0$ on distorted elements | Check $\det(\mathbf{J}) > 0$ at every Gauss point; throw exception if distorted |
| **Incorrect pressure direction** | Incorrect face normal orientation | Follow 1-based face ID conventions and ensure positive $p$ produces compressive normal forces |
| **Wrong 2D stiffness scale** | Missing `thickness` factor | Multiply 2D $\mathbf{K}_e$ by `ctx.thickness` |
| **Stress post-processing errors** | $B$ matrix not stacked properly in `ElementMatrices` | Ensure rows $[g \cdot n_{\text{strain}}, (g+1) \cdot n_{\text{strain}} - 1]$ store Gauss point $g$'s matrix |
| **Coordinate array out-of-bounds** | Assuming `coords` has 2 columns in 2D | Always query `.col(0)`, `.col(1)`, `.col(2)` since `coords` is fixed at 3 columns |
