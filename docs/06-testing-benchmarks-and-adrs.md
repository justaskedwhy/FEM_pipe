# Package 6: Testing, Benchmarks & Architectural Decision Records Specification
**Document Version:** 1.0  
**Target Files:** `docs/testing.md`, `docs/design_decisions.md`, `docs/benchmarks.md`, `.github/workflows/ci.yml`  
**Grounded Sources:** `docs/testing.md`, `docs/design_decisions.md`, `docs/benchmarks.md`, `.github/workflows/ci.yml`

---

## 1. Testing Framework & Verification Strategy

### 1.1 Test Hierarchy & Pytest Markers
The test suite is written entirely in Python using `pytest`, exercising the C++ numerical core through `femcore` bindings. This establishes a single test runner and unified coverage reporting.

```text
python/tests/
├── conftest.py                # Shared fixtures and mesh generators
├── test_parser.py             # M1 unit tests
├── test_model_builder.py      # M2 unit tests
├── test_preprocessor.py       # M3 unit tests
├── test_element_kernels.py    # M4 kernel unit tests & patch tests
├── test_loads.py              # M5 load engine unit tests
├── test_assembler.py          # M6 assembler unit tests
├── test_bc.py                 # M7 boundary condition handler unit tests
├── test_solver.py             # M8 linear solver unit tests
├── test_postprocess.py        # M9 post-processor unit tests
├── test_vtk_writer.py         # M10 VTK writer unit tests
├── test_swap_q4_t3.py         # Mandatory Q4 vs T3 element swap integration test
└── test_end_to_end.py         # CLI end-to-end execution tests
```

#### Pytest Markers (`pyproject.toml`)
* `@pytest.mark.unit`: Fast, isolated component tests ($< 1\text{s}$).
* `@pytest.mark.integration`: Cross-module workflow tests ($1\text{s} - 10\text{s}$).
* `@pytest.mark.benchmark`: Analytical validation benchmarks ($10\text{s} - 120\text{s}$).
* `@pytest.mark.slow`: Heavy mesh convergence runs ($> 5\text{min}$, nightly CI).

---

### 1.2 The Single-Element Patch Test Protocol
Every `ElementKernel` subclass must pass the patch test before being declared valid.

#### Procedure
1. Create a 1-element or 4-element mesh of the candidate element type.
2. Apply prescribed displacement boundary conditions matching an exact linear field:
   $$\mathbf{u}(\mathbf{x}) = \mathbf{A} \mathbf{x} + \mathbf{b}$$
3. **Acceptance Criteria**:
   * Internal free nodes recover exact displacements to $< 10^{-10}$.
   * Element strain $\boldsymbol{\varepsilon} = \mathbf{B} \mathbf{U}_e$ and stress $\boldsymbol{\sigma} = \mathbf{D} \boldsymbol{\varepsilon}$ are uniform and match analytical values to $< 10^{-10}$.
   * Rigid body translations ($\mathbf{u} = \text{const}$) and small rotations produce zero element energy:
     $$\| \mathbf{K}_e \mathbf{u}_{\text{rigid}} \| / \| \mathbf{K}_e \| < 10^{-12}$$

---

### 1.3 The Mandatory Element Swap Test (`test_swap_q4_t3.py`)
To prove that the `ElementKernel` swappability contract (Invariant **I1**) holds:
1. Solve the 2D Cantilever Beam problem using bilinear quads (`CPS4`).
2. Swap the mesh to linear triangles (`CPS3`) by editing **only** the `*ELEMENT` keyword in the `.inp` deck.
3. Solve without modifying any Python orchestration or C++ numerical code.
4. **Acceptance Criteria**:
   * Both element models execute cleanly without code changes.
   * Both element types converge monotonically to the analytical tip deflection under mesh refinement.
   * The $L^2$ error norm $\| \mathbf{u}_{\text{Q4}} - \mathbf{u}_{\text{T3}} \|$ decreases monotonically as $h \to 0$.

---

## 2. Analytical Benchmarks & Convergence Suite

### 2.1 Benchmark 1: Cantilever Beam (2D Plane Stress)

#### Problem Specification
* Geometry: Length $L = 1.0\text{ m}$, Height $H = 0.1\text{ m}$, Thickness $t = 0.01\text{ m}$.
* Material: Steel ($E = 210\text{ GPa}$, $\nu = 0.3$).
* Boundary: Fixed at $x = 0$ ($u_x = u_y = 0$).
* Loading: Downward tip load $P = -100\text{ N}$ distributed on $x = L$.

#### Analytical Solution (Euler-Bernoulli)
$$I = \frac{t H^3}{12} = \frac{0.01 \times 0.1^3}{12} = 8.3333 \times 10^{-7}\text{ m}^4$$
$$\delta_{\text{tip}} = \frac{P L^3}{3 E I} = \frac{100 \times 1.0^3}{3 \times 210 \times 10^9 \times 8.3333 \times 10^{-7}} = 1.9048 \times 10^{-4}\text{ m}$$

#### Refinement Results Template
| Mesh ID | Grid Size ($n_x \times n_y$) | Elements | Total DOFs | FEM Tip Deflection (m) | Relative Error (%) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| M1 | $10 \times 2$ | 20 (Q4) | 66 | $1.3545 \times 10^{-4}$ | $-28.89\%$ |
| M2 | $20 \times 2$ | 40 (Q4) | 126 | $1.7054 \times 10^{-4}$ | $-10.47\%$ |
| M3 | $40 \times 4$ | 160 (Q4) | 410 | $1.8577 \times 10^{-4}$ | $-2.47\%$ |
| M4 | $80 \times 8$ | 640 (Q4) | 1458 | $1.9010 \times 10^{-4}$ | $-0.20\%$ |

*Monotone convergence from below is observed as shear locking diminishes with refinement.*

> **Template-history note**: the values above are produced by the pipeline's shipped Q4 kernel (doc 02 verbatim full-integration 2×2 Gauss bilinear quad with shear locking). An earlier revision of this table listed M1–M3 as (`1.810e-4`, `1.847e-4`, `1.888e-4`) — those were generated with a reduced-integration or otherwise non-locking formulation and were inconsistent with the doc-governed kernel. `tests/run_benchmarks.py` enforces the correctness invariants that matter here: every mesh under-predicts the analytical tip and convergence is monotone from below (doc 06 §1.3).

---

### 2.2 Benchmark 2: Plate with a Hole (2D Plane Stress)

#### Problem Specification
* Geometry: Half-width $W = 2.0\text{ m}$, Half-height $H = 2.0\text{ m}$, Hole radius $a = 0.5\text{ m}$, Thickness $t = 1.0\text{ m}$.
* Material: Aluminum ($E = 200\text{ GPa}$, $\nu = 0.3$).
* Boundary: Symmetry BCs at $x = 0$ ($u_x = 0$) and $y = 0$ ($u_y = 0$).
* Loading: Uniaxial tensile stress $\sigma_0 = 1.0\text{ MPa}$ on right edge ($x = W$).

#### Analytical Solution (Kirsch & Howland Correction)
Infinite plate Kirsch solution yields stress concentration factor $K_t = \frac{\sigma_{\max}}{\sigma_0} = 3.0$. For finite width ratio $W/a = 4.0$, Howland's correction yields $K_t \approx 3.0 - 3.1$.

#### Refinement Results Template
| Mesh ID | Element Size at Hole ($h$) | Elements | Total DOFs | FEM $\sigma_{\max}$ (MPa) | Computed $K_t$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| M1 | $a / 4$ | 320 (T3) | 380 | $2.78$ | $2.78$ |
| M2 | $a / 8$ | 1280 (T3) | 1420 | $2.94$ | $2.94$ |
| M3 | $a / 16$ | 5120 (T3) | 5400 | $3.03$ | $3.03$ |

---

### 2.3 Benchmark 3: 3D Block under Uniaxial Compression

#### Problem Specification
* Geometry: Unit cube ($1.0 \times 1.0 \times 1.0\text{ m}$).
* Material: Aluminum ($E = 70\text{ GPa}$, $\nu = 0.33$).
* Loading: Strain-controlled prescribed displacement $u_z = -0.001\text{ m}$ on top face ($z = 1$). Fixed bottom face ($z = 0, u_z = 0$). Rigid body modes constrained.

#### Analytical Solution
$$\sigma_{zz} = E \cdot \varepsilon_{zz} = 70 \times 10^9 \times (-0.001) = -70.0\text{ MPa}$$
$$\varepsilon_{xx} = \varepsilon_{yy} = -\nu \varepsilon_{zz} = 0.33 \times 0.001 = +3.3 \times 10^{-4}$$

#### Verification
Even a $2 \times 2 \times 2$ mesh of H8 or T4 elements reproduces $\sigma_{zz} = -70.0\text{ MPa}$ and $\varepsilon_{xx} = 3.3 \times 10^{-4}$ to machine precision ($< 10^{-12}$ error).

---

## 3. Architectural Decision Records (ADRs)

| ADR ID | Decision Title | Status | Summary & Rationale |
| :--- | :--- | :--- | :--- |
| **ADR-0001** | Hybrid Python + C++ | Accepted | Python handles I/O, parsing, and orchestration; C++ executes numerical loops for 10–100$\times$ performance. |
| **ADR-0002** | Single Language Boundary | Accepted | Exactly two boundary crossings (`SolverInput` in, `SolverOutput` out). Eliminates per-element callback overhead. |
| **ADR-0003** | Element Kernel Interface | Accepted | Decouples element formulations from solver logic via polymorphic `ElementKernel` and static registration macros. |
| **ADR-0004** | CSR Connectivity Layout | Accepted | Connectivity is stored as contiguous 1D array + offsets array (`int32`), eliminating zero-padding in mixed meshes. |
| **ADR-0005** | Eigen for Linear Algebra | Accepted | Uses header-only Eigen library for dense element matrices, sparse assembly, and `SimplicialLDLT` system solves. |
| **ADR-0006** | pybind11 for Glue | Accepted | Header-only bindings providing zero-copy NumPy array conversions and native Eigen interop. |
| **ADR-0007** | ASCII VTK Output | Accepted | VTK XML `.vtu` written in ASCII for V1.0 to enable human readability, text diffing, and easy ParaView debugging. |
| **ADR-0008** | Sparse LDLT Primary Solver | Accepted | `Eigen::SimplicialLDLT` primary solver exploits stiffness matrix symmetry; falls back to `Eigen::SparseLU`. |
| **ADR-0009** | ABAQUS Name Identification| Accepted | ABAQUS canonical string (`CPS4`, `C3D8`) is the single element identifier across Python and C++ layers. |
| **ADR-0010** | Stateless Element Kernels | Accepted | `compute()` calls are thread-safe and stateless, receiving all parameters via `ElementContext`. |
| **ADR-0011** | Registry Global State Only | Accepted | `ElementRegistry` is the only global singleton, populated during static initialization and immutable during solves. |
| **ADR-0012** | Fixed Quadrature Rules | Accepted | Quadrature is fixed per element type (1-point for T3/T4, $2\times2$/$2\times2\times2$ for Q4/H8). |
| **ADR-0013** | Input Order Preservation | Accepted | Preserves input node and element order throughout the pipeline to simplify output verification and debugging. |
| **ADR-0014** | No C++ Logging Library | Accepted | C++ core returns numerical diagnostics in `SolverOutput` and throws exceptions rather than printing to stdout/stderr. |
| **ADR-0015** | Python-Based Test Suite | Accepted | All unit, integration, and benchmark tests are written in Python (`pytest`) exercising compiled `femcore`. |
| **ADR-0016** | Documentation-First | Accepted | Comprehensive specification documents are written prior to code implementation to freeze interfaces. |

---

## 4. Automated CI Pipeline (`.github/workflows/ci.yml`)

```yaml
name: CI

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]

jobs:
  linux:
    name: Linux (Ubuntu 22.04, Python ${{ matrix.python-version }})
    runs-on: ubuntu-22.04
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.10", "3.11", "3.12"]

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install System Dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y build-essential cmake ninja-build libeigen3-dev pybind11-dev

      - name: Install Python Dependencies
        run: |
          python -m pip install --upgrade pip
          pip install numpy pytest pytest-cov meshio matplotlib

      - name: Build & Install femcore
        run: pip install . -v

      - name: Verify C++ Module Import
        run: python -c "import femcore; print(femcore.registered_elements())"

      - name: Run Test Suite
        run: pytest python/tests -m "not slow" -v --cov=fem --cov-report=xml

      - name: Run Analytical Benchmarks
        run: python tests/run_benchmarks.py

  windows:
    name: Windows (Server 2022, Python ${{ matrix.python-version }})
    runs-on: windows-2022
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.10", "3.11", "3.12"]

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Set up MSVC Toolchain
        uses: ilammy/msvc-dev-cmd@v1

      - name: Install vcpkg Dependencies
        run: |
          git clone https://github.com/microsoft/vcpkg $env:VCPKG_ROOT
          & "$env:VCPKG_ROOT\bootstrap-vcpkg.bat"
          & "$env:VCPKG_ROOT\vcpkg.exe" install eigen3 pybind11 --triplet x64-windows

      - name: Build & Install
        run: |
          $env:CMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"
          pip install . -v

      - name: Run Unit Tests
        run: pytest python/tests -m "not slow" -v

  lint:
    name: Code Style & Quality Checks
    runs-on: ubuntu-22.04
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: Install Linters
        run: pip install flake8 black
      - name: Python Format Check
        run: black --check --diff python/fem
      - name: Python Lint Check
        run: flake8 python/fem --max-line-length=100
```
