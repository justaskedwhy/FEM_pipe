# Package 5: Build, Execution & Troubleshooting Guides Specification
**Document Version:** 1.0  
**Target Files:** `README.md`, `docs/beginner_guide.md`, `docs/troubleshooting.md`, `docs/contributing.md`  
**Grounded Sources:** `docs/beginner_guide.md`, `docs/troubleshooting.md`, `docs/contributing.md`

---

## 1. Prerequisites & Environment Setup

### 1.1 Requirements Matrix
| Tool / Library | Version Requirement | Purpose | Availability / License |
| :--- | :--- | :--- | :--- |
| **Python** | $\ge 3.10$ | Orchestration & I/O (Modules M1, M2, M3, M10) | PSF License |
| **NumPy** | $\ge 1.23$ | Contiguous numeric buffers | BSD License |
| **CMake** | $\ge 3.16$ | Cross-platform C++ build system | BSD License |
| **C++ Compiler** | C++17 Standard | C++ Numerical Core compilation | MSVC 2019/2022, GCC 9+, Clang 10+ |
| **Eigen** | $\ge 3.4$ | Sparse/dense linear algebra (header-only) | MPL2 License |
| **pybind11** | $\ge 2.11$ | Python $\leftrightarrow$ C++ memory bridge | BSD License |
| **meshio** | $\ge 5.0$ | VTK XML `.vtu` exporter | MIT License |
| **pytest** | $\ge 7.0$ | Automated test execution | MIT License |

---

### 1.2 Linux Installation (Ubuntu 20.04 / 22.04+)

```bash
# 1. Update system package index and install C++ toolchain & headers
sudo apt update
sudo apt install -y build-essential cmake python3-dev python3-pip \
                    libeigen3-dev pybind11-dev ninja-build

# 2. Install Python dependencies
pip install --upgrade pip
pip install numpy pytest meshio matplotlib
```

---

### 1.3 Windows Installation (Visual Studio & `vcpkg`)

1. **Install Python 3.10+**: Download from `python.org` and ensure **"Add Python to PATH"** is checked during installation.
2. **Install Visual Studio 2019 or 2022**: Include the **"Desktop development with C++"** workload.
3. **Install CMake**: Download installer from `cmake.org` or install via winget:
   ```cmd
   winget install Kitware.CMake
   ```
4. **Install Header Libraries via `vcpkg`**:
   ```cmd
   git clone https://github.com/microsoft/vcpkg %USERPROFILE%\vcpkg
   cd %USERPROFILE%\vcpkg
   .\bootstrap-vcpkg.bat
   .\vcpkg install eigen3 pybind11 --triplet x64-windows
   .\vcpkg integrate install
   ```

---

## 2. Compilation, Verification & Execution

### 2.1 Source Retrieval & Package Build

```bash
# Clone source repository
git clone <repo-url> fem-project
cd fem-project

# Build C++ extension module (femcore) and install Python package (fem)
pip install . -v
```

---

### 2.2 Post-Build Verification Protocol

Verify that static element initialization macros successfully populated the `ElementRegistry`:

```bash
python -c "import femcore; print(femcore.registered_elements())"
```

#### Expected Output
```text
['C3D4', 'C3D8', 'C3D8R', 'CPE3', 'CPE4', 'CPS3', 'CPS4', 'CPS4R']
```
*If `femcore.registered_elements()` returns an empty list, static initializers were stripped by the linker. Ensure every `elements_*.cpp` file is listed under `pybind11_add_module(femcore ...)` in `CMakeLists.txt`.*

---

### 2.3 Command-Line Interface (CLI) Usage

```bash
# General syntax
femsolver --in <input_deck.inp> --out <output_mesh.vtu> [--verbose]

# Example 1: Cantilever Beam benchmark
femsolver --in tests/inputs/cantilever.inp --out solution_cantilever.vtu

# Example 2: Plate with a Hole benchmark
femsolver --in tests/inputs/plate_hole.inp --out solution_plate.vtu
```

#### Expected Console Log
```text
[parser]  Read 121 nodes, 100 elements from tests/inputs/cantilever.inp
[model]   Dimension = 2, active materials = 1
[preproc] DOFs total = 242, free = 238, prescribed = 4
[solve]   Assembled K (sparse 242x242, nnz=2116)
[solve]   SimplicialLDLT solve OK, residual norm = 1.24e-15
[post]    Max nodal displacement magnitude |u|_max = 1.847e-03 m
[vtu]     Wrote solution_cantilever.vtu
```

---

## 3. Diagnostic & Troubleshooting Matrix

### 3.1 Build & Installation Failures

| Symptom / Error | Primary Cause | Resolution |
| :--- | :--- | :--- |
| `ImportError: No module named femcore` | C++ extension module not compiled | Run `pip install . -v` from repository root. Check CMake compiler detection logs. |
| `CMake error: Could not find Eigen3` | Missing `Eigen3Config.cmake` | **Linux**: `sudo apt install libeigen3-dev`<br>**Windows**: Pass `-DCMAKE_TOOLCHAIN_FILE=%USERPROFILE%/vcpkg/scripts/buildsystems/vcpkg.cmake` |
| `MSVC: unresolved external symbol` | Source file omitted from target | Ensure all `cpp/src/*.cpp` files are explicitly listed in `CMakeLists.txt`. |
| `registered_elements()` returns `[]` | Static registration objects stripped | Verify `REGISTER_ELEMENT` macros are present at the bottom of element `.cpp` files and linked into `femcore`. |

---

### 3.2 Input Deck Parsing & Preprocessing Errors

| Symptom / Error | Primary Cause | Resolution |
| :--- | :--- | :--- |
| `InputError: Malformed float at line N` | Non-standard float format | Check line N for comma decimals (`1,5` vs `1.5`) or missing comma-separated fields. |
| `ModelError: Node X referenced by element Y does not exist` | Missing node definition | Ensure all node IDs referenced in `*ELEMENT` blocks exist in `*NODE` blocks. |
| `ModelError: Element type 'CPS6' is not registered` | Unsupported element keyword | Replace with supported elements (`CPS3`, `CPS4`) or implement the requested kernel. |
| `ModelError: Mixed dimension: 2D element in 3D model` | Inconsistent mesh dimensions | Ensure 2D elements (`CPS*`, `CPE*`) and 3D elements (`C3D*`) are not mixed in a single deck. |

---

### 3.3 Numerical & Linear Solver Failures

| Symptom / Error | Primary Cause | Resolution |
| :--- | :--- | :--- |
| `femcore: negative Jacobian in element N` | Inverted or degenerate element | Check node ordering in `.inp` file. 2D quads must be counterclockwise; 3D hexes must follow right-hand rule. |
| `femcore: linear solve failed (LDLT and LU)` | Singular stiffness matrix $\mathbf{K}$ | Insufficient Dirichlet boundary conditions. Apply BCs to eliminate rigid body modes (3 in 2D, 6 in 3D). |
| `residual_norm > 1e-8` | Matrix ill-conditioning or extreme load scales | Check material properties ($E > 0, -1 < \nu < 0.5$) and avoid unphysically massive nodal forces. |

---

### 3.4 Post-Processing & Visualization Issues

| Symptom / Error | Primary Cause | Resolution |
| :--- | :--- | :--- |
| ParaView shows empty window | Blank or invalid `.vtu` XML | Inspect `.vtu` file size (>0 bytes) and verify XML header formatting. |
| ParaView shows mesh with no deformation | Unwarped displacement field | Apply **Warp By Vector** filter in ParaView with `Vectors = displacement` and set scale factor. |
| Displacements are unphysically large | Inconsistent unit system | Ensure $E$ ($Pa$), geometry ($m$), and forces ($N$) are in a consistent unit system (e.g., SI). |

---

## 4. Contributing & Development Guidelines

### 4.1 Git Branching Model
* `main`: Production-ready, fully tested releases.
* `develop`: Integration branch for upcoming features.
* `feature/<feature-name>`: Topic branch for new functionality.
* `bugfix/<issue-name>`: Target fix for specific reported issues.

---

### 4.2 Formatting & Linting Standards
Developers must apply automated formatting prior to committing code:

#### Python Code Standards
* **Formatter**: `black` (line length 88).
* **Linter**: `flake8` (`--max-line-length=100`).
* **Command**:
  ```bash
  black python/fem python/tests
  flake8 python/fem --max-line-length=100
  ```

#### C++ Code Standards
* **Style**: `clang-format` using K&R 4-space indentation.
* **Header Guard**: `#pragma once` on all public headers.
* **Command**:
  ```bash
  clang-format -i cpp/include/*.hpp cpp/src/*.cpp cpp/bindings/*.cpp
  ```

---

### 4.3 Pull Request Checklist
Before submitting a pull request, verify:
- [ ] Code compiles cleanly on both Linux (GCC) and Windows (MSVC).
- [ ] All unit and integration tests pass (`pytest python/tests -v`).
- [ ] Code formatting meets `black` and `clang-format` rules.
- [ ] New functionality includes dedicated unit tests in `python/tests/`.
- [ ] Relevant documentation files (`docs/*.md`) have been updated.
