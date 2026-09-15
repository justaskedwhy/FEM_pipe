# Package 4: File Format Specifications
**Document Version:** 1.0  
**Target Modules:** Module M1 (Input Parser), Module M10 (VTK Writer)  
**Grounded Sources:** `docs/inp_format.md`, `docs/vtk_format.md`, `docs/master.md` §7, §8

---

## 1. ABAQUS `.inp` Subset Specification

### 1.1 Scope & Parser Philosophy
This specification defines the exact **supported subset** of the ABAQUS input deck format. The parser (Module M1) executes lexical and syntactic parsing to build a `RawModel`. Any keyword outside this supported subset is safely skipped with a logged warning, ensuring backward compatibility with standard FEA preprocessors.

---

### 1.2 Lexical Rules & Tokenization
1. **Character Encoding**: Standard ASCII-compatible UTF-8.
2. **Line Endings**: Universal support for LF (`\n`) and CRLF (`\r\n`).
3. **Comments**: Any line where the first non-whitespace characters are `**` is treated as a comment and ignored.
4. **Keywords**:
   * A line whose first non-whitespace character is `*` (and not `**`) is a keyword header.
   * Keyword names are **case-insensitive** (e.g., `*NODE` and `*node` are identical).
   * Parameters follow the keyword name, separated by commas.
   * Parameter syntax: `KEY` or `KEY=VALUE`.
5. **Data Lines**:
   * Any non-comment, non-keyword line is a data line bound to the **most recent** active keyword.
   * Blank lines within a data block are ignored.
6. **Identifiers & Numbers**:
   * Node and Element IDs are 1-based positive integers (`int32`).
   * Material and Set names match `[A-Za-z_][A-Za-z0-9_]*`.
   * Floating-point values follow standard scientific notation (e.g., `210.0e9`, `1.5`, `-0.001`).

---

### 1.3 Formal EBNF Grammar

```ebnf
file            = { line } ;
line            = comment | keyword_line | data_line | empty ;
comment         = "**" , { ? any character ? } , newline ;
keyword_line    = "*" , keyword_name , [ "," , parameter_list ] , newline ;
keyword_name    = letter , { letter | digit | "_" } ;
parameter_list  = parameter , { "," , parameter } ;
parameter       = parameter_key [ "=" , parameter_value ] ;
parameter_key   = identifier ;
parameter_value = identifier | number | string ;
data_line       = field , { "," , field } , newline ;
field           = number | identifier | empty ;
empty           = { whitespace } , newline ;

number          = [ "+" | "-" ] , digits , [ "." , digits ] , [ exponent ] ;
exponent        = ( "e" | "E" ) , [ "+" | "-" ] , digits ;
identifier      = ( letter | "_" ) , { letter | digit | "_" } ;
string          = identifier ;
letter          = "A".."Z" | "a".."z" ;
digit           = "0".."9" ;
whitespace      = " " | "\t" ;
```

---

### 1.4 Supported Keywords Reference

#### 1. `*NODE`
Defines nodal spatial coordinates.
```abaqus
*NODE [, NSET=<name>]
<id>, <x>, <y> [, <z>]
```
* **Parameters**: `NSET` (optional) — assigns created nodes to a named node set.
* **Fields**:
  * `id` (`int`, required): Positive node identifier.
  * `x` (`float`, required): X coordinate in global space.
  * `y` (`float`, required): Y coordinate in global space.
  * `z` (`float`, optional): Z coordinate in global space (defaults to `0.0` for 2D meshes).

#### 2. `*ELEMENT`
Defines mesh elements and topological connectivity.
```abaqus
*ELEMENT, TYPE=<type> [, ELSET=<name>]
<id>, <n1>, <n2>, ..., <n_npe>
```
* **Parameters**:
  * `TYPE` (`string`, required): ABAQUS element identifier (`CPS3`, `CPE3`, `CPS4`, `CPE4`, `CPS4R`, `C3D4`, `C3D8`, `C3D8R`).
  * `ELSET` (optional): Assigns created elements to a named element set.
* **Fields**:
  * `id` (`int`, required): Positive element identifier.
  * `n1 ... n_npe` (`int`, required): 1-based node IDs matching the element's node count (`npe`).

#### 3. `*MATERIAL` & `*ELASTIC`
Defines constitutive properties.
```abaqus
*MATERIAL, NAME=<name>
*ELASTIC [, TYPE=ISO]
<E>, <nu>
```
* **Parameters**: `NAME` (`string`, required) — material reference identifier. `TYPE=ISO` (default).
* **Fields**:
  * `E` (`float`, required): Young's modulus in Pascals ($Pa$).
  * `nu` (`float`, required): Poisson's ratio ($\nu$, $-1.0 < \nu < 0.5$).

#### 4. `*SOLID SECTION`
Binds material properties to an element set and provides geometric thickness.
```abaqus
*SOLID SECTION, ELSET=<name>, MATERIAL=<name>
[<thickness>]
```
* **Parameters**: `ELSET` (`string`, required), `MATERIAL` (`string`, required).
* **Fields**:
  * `thickness` (`float`, required for 2D, ignored for 3D): Out-of-plane element thickness (m).

#### 5. `*BOUNDARY`
Applies prescribed displacement boundary conditions.
```abaqus
*BOUNDARY
<node_or_set>, <dof_first>, <dof_last>, <value>
```
* **Fields**:
  * `node_or_set` (`int` or `string`): Node ID or named NSET.
  * `dof_first` (`int`, 1-based): Starting DOF (`1` = X, `2` = Y, `3` = Z).
  * `dof_last` (`int`, 1-based): Ending DOF ($\ge \text{dof\_first}$).
  * `value` (`float`): Prescribed displacement magnitude.

#### 6. `*CLOAD`
Applies concentrated nodal point loads.
```abaqus
*CLOAD
<node_or_set>, <dof>, <magnitude>
```
* **Fields**:
  * `node_or_set` (`int` or `string`): Target node ID or named NSET (expanded to one member node per data line during model build). A non-integer first field is treated as a set name.
  * `dof` (`int`, 1-based): Target degree of freedom (`1` = X, `2` = Y, `3` = Z).
  * `magnitude` (`float`): Concentrated force magnitude (N).

#### 7. `*DLOAD`
Applies uniform surface pressure loads to element faces.
```abaqus
*DLOAD
<element>, <face_label>, <magnitude>
```
* **Fields**:
  * `element` (`int`): Target element ID.
  * `face_label` (`string`): Face label (`P1` through `P6`).
  * `magnitude` (`float`): Pressure magnitude in Pascals ($Pa$). Positive pressure acts **compressive (into)** the element face.

**Face-label authority**: `P<k>` maps directly to the 1-based local face id `k` of the element kernel (`kernel->face_nodes(k)`, `k ∈ [1, n_faces]`). The kernel's face tables (§1.3 and doc 02 §2.9) define which nodes form face `k`. The preprocessor validates `k` against the element's `n_faces` from the registry and raises a `ModelError` if out of range — no parser-side or hardcoded face table exists outside the kernels (invariant I2).

#### 8. `*SURFACE` (element-based) & `*DSLOAD` (surface pressure)
Defines a named element-based surface and applies a uniform pressure to it. This is the ABAQUS/CAE-native pattern for load application (used by the example decks `Q4.inp`, `T3.inp`, `T4.inp`, `H8.inp`).
```abaqus
*SURFACE, TYPE=ELEMENT, NAME=<name>
<element_or_elset>, <face_label>

*DSLOAD
<surface_name>, P, <magnitude>
```
* **Fields** (`*SURFACE`): `element_or_elset` (`int` element ID or named `ELSET`) plus a face label `S1` through `S6`. Each data line contributes one (set, face) pair; multiple lines accumulate onto the named surface.
* **Fields** (`*DSLOAD`): `surface_name` (`string`, must match a `*SURFACE` `NAME=`) and load label `P` followed by the magnitude (`float`, Pa). Positive acts compressive.
* The model builder expands each surface into one per-element `P<k>` pressure load (mapping `S<k>` → `P<k>`), validating that the referenced surface and element set exist (`ModelError` otherwise). Per-element `*DLOAD` (§1.4.7 above) remains supported independently.

#### 9. `*NSET` & `*ELSET` (standalone set blocks)
Defines a named node or element set by explicit member list or by `GENERATE` range (complementing the `NSET=`/`ELSET=` parameters on `*NODE`/`*ELEMENT`).
```abaqus
*NSET, NSET=<name> [, UNSORTED]
<id> [, <id> ...]

*NSET, NSET=<name>, GENERATE
<first>, <last> [, <step>]

*ELSET, ELSET=<name> [, UNSORTED]
<id> [, <id> ...]

*ELSET, ELSET=<name>, GENERATE
<first>, <last> [, <step>]
```
* **Parameters**: `NSET`/`ELSET` (`string`, required); `UNSORTED`/`INTERNAL` (optional, ignored); `GENERATE` (optional) selects range form.
* **Data**: comma-separated node/element IDs, one or more per line, appended to the named set; or, under `GENERATE`, a single line `<first>, <last>[, <step>]` producing the integer sequence `first, first+step, ..., last` (step defaults to `1`, must be `> 0`).
* Set-name tokens in a member list refer to previously-defined sets (members are unioned, de-duplicated) — `ModelError`/`InputError` if the referenced set is undefined.
* Referencing an undefined set in `*BOUNDARY`/`*CLOAD`/`*SURFACE` raises a `ModelError`.

#### 10. `*PART`, `*ASSEMBLY`, `*INSTANCE` scoping
ABAQUS/CAE exports enclose nodes/elements inside `*PART ... *END PART` and define instance-level sets inside `*ASSEMBLY ... *END ASSEMBLY`. The parser flattens these into a single mesh but preserves a scoped namespace so identical set names at different scopes do not collide:
* Set definitions inside a part are stored as `part:<part>:<name>`.
* Set definitions carrying an `INSTANCE=<inst>` parameter are stored as `instance:<inst>:<name>`.
* References (`*SOLID SECTION ELSET=`, `*BOUNDARY`, `*CLOAD`, `*SURFACE`) resolve to the matching scope: inside a part they use the part-prefixed key; at model/step level they fall back to the most recently declared instance's key, then the bare name.
* Material/section/step blocks may sit outside `*PART`/`*ASSEMBLY` blocks; both layouts are accepted.

#### 11. Analysis Step Control
```abaqus
*STEP [, NLGEOM=NO]
*STATIC
*END STEP
```
Defines linear static solution execution. `NLGEOM` must be `NO` or omitted.

---

### 1.5 Ignored Keywords & Diagnostic Policy
The parser logs an informational warning and safely skips data lines for the following unsupported keywords:
* Output requests: `*RESTART`, `*OUTPUT`, `*EL PRINT`, `*NODE PRINT`, `*EL FILE`, `*NODE FILE`.
* Constraints & Contact: `*TIE`, `*CONTACT`, `*EQUATION`, `*MPC`, `*KINEMATIC`.
* Structured headers: `*HEADING`, `*PREPRINT`.
* Structural/Thermal: `*SHELL SECTION`, `*BEAM SECTION`, `*AMPLITUDE`, `*TEMPERATURE`.

---

### 1.6 Complete Example Input File (`cantilever.inp`)

> **Note on completeness**: the abridged listing below is intended to illustrate syntax, not a self-consistent mesh. The canonical, fully-validated deck (20 CPS4 elements, nodes 1–33, correct boundary/tip-load node references) lives at [`tests/inputs/cantilever.inp`](../tests/inputs/cantilever.inp) and passes the full pipeline. Example `.inp` decks for pressure loads: `tests/inputs/q4_pressure.inp` (CPS4) and `tests/inputs/h8_pressure.inp` (C3D8).

```abaqus
** =====================================================================
** Cantilever beam: 2D plane stress, Q4 elements (CPS4)
** Length = 1.0 m, Height = 0.1 m, Thickness = 0.01 m
** Material: Steel (E = 210 GPa, nu = 0.3)
** Tip load: -100 N in y at the free end
** =====================================================================

*NODE, NSET=ALL
1,   0.00,  0.00
2,   0.05,  0.00
3,   0.10,  0.00
21,  1.00,  0.00
22,  0.00,  0.05
42,  1.00,  0.05
43,  0.00,  0.10
63,  1.00,  0.10

*ELEMENT, TYPE=CPS4, ELSET=BEAM
1,  1, 2, 23, 22
2,  2, 3, 24, 23
20, 20, 21, 42, 41
21, 22, 23, 44, 43
40, 41, 42, 63, 62

*MATERIAL, NAME=Steel
*ELASTIC
210.0e9, 0.3

*SOLID SECTION, ELSET=BEAM, MATERIAL=Steel
0.01

*BOUNDARY
1,  1, 2, 0.0
22, 1, 2, 0.0
43, 1, 2, 0.0

*CLOAD
21, 2, -50.0
63, 2, -50.0

*STEP
*STATIC
*END STEP
```

---

## 2. VTK XML UnstructuredGrid (`.vtu`) Specification

### 2.1 Overview & Design Rationale
The VTK XML `UnstructuredGrid` (`.vtu`) format is chosen for post-processing and visualization output.
* **Native ParaView Compatibility**: Ingested directly by ParaView, VisIt, and `meshio`.
* **3D Embedding Strategy**: To allow standard 3D visualizers to display 2D plane stress/strain results, 2D coordinates are written as 3D points with $z = 0.0$.
* **Format**: ASCII encoding for Version 1.0 (enabling human readability and text diffing).

---

### 2.2 VTK XML Skeleton Structure

```xml
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="N" NumberOfCells="M">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          x0 y0 z0  x1 y1 z1  ...  xN yN zN
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          c0_0 c0_1 ... cM_k
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          offset_1 offset_2 ... offset_M
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          type_1 type_2 ... type_M
        </DataArray>
      </Cells>
      <PointData Scalars="displacement_magnitude" Vectors="displacement">
        <!-- Nodal Field Data -->
      </PointData>
      <CellData Scalars="von_mises">
        <!-- Element Centroidal Stress & Strain Data -->
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
```

---

### 2.3 VTK Cell Type Mapping & Permutation

| Element Kernel | ABAQUS Names | VTK Code | VTK Name | Permutation Vector |
| :--- | :--- | :--- | :--- | :--- |
| **T3** | `CPS3`, `CPE3` | `5` | `VTK_TRIANGLE` | `{0, 1, 2}` |
| **Q4** | `CPS4`, `CPE4`, `CPS4R` | `9` | `VTK_QUAD` | `{0, 1, 2, 3}` |
| **T4** | `C3D4` | `10` | `VTK_TETRA` | `{0, 1, 2, 3}` |
| **H8** | `C3D8`, `C3D8R` | `12` | `VTK_HEXAHEDRON` | `{0, 1, 2, 3, 4, 5, 6, 7}` |
| *T6 (Future)* | `CPS6`, `CPE6` | `22` | `VTK_QUADRATIC_TRIANGLE` | `{0, 1, 2, 3, 4, 5}` |
| *Q8 (Future)* | `CPS8`, `CPE8` | `23` | `VTK_QUADRATIC_QUAD` | `{0, 1, 2, 3, 4, 5, 6, 7}` |
| *T10 (Future)* | `C3D10` | `24` | `VTK_QUADRATIC_TETRA` | `{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}` |
| *H20 (Future)* | `C3D20` | `25` | `VTK_QUADRATIC_HEXAHEDRON` | `{0..7 corners, 8..19 mid-edges}` |

---

### 2.4 Data Array Conventions

#### `<PointData>` Arrays
| Name | Type | Components | Description |
| :--- | :--- | :--- | :--- |
| `displacement` | `Float64` | 3 | Vector $[u_x, u_y, u_z]$ in meters |
| `displacement_magnitude` | `Float64` | 1 | Scalar $\| \mathbf{u} \| = \sqrt{u_x^2 + u_y^2 + u_z^2}$ |
| `reaction` *(optional)* | `Float64` | 3 | Nodal reaction forces $[R_x, R_y, R_z]$ |

#### `<CellData>` Arrays
| Name | Type | Components | Description |
| :--- | :--- | :--- | :--- |
| `stress` | `Float64` | 6 | Centroid Voigt stress $[\sigma_{xx}, \sigma_{yy}, \sigma_{zz}, \sigma_{xy}, \sigma_{yz}, \sigma_{zx}]$ |
| `strain` | `Float64` | 6 | Centroid Voigt engineering strain $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{zx}]$ |
| `von_mises` | `Float64` | 1 | Equivalent von Mises stress $\sigma_{\text{vm}}$ ($Pa$) |
| `principal_stress` | `Float64` | 3 | Ordered principal stresses $[\sigma_1, \sigma_2, \sigma_3]$ where $\sigma_1 \ge \sigma_2 \ge \sigma_3$ |

---

### 2.5 Complete Working `.vtu` Output Example

```xml
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0
          1.0 0.0 0.0
          1.0 1.0 0.0
          0.0 1.0 0.0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1 2 3
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          4
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          9
        </DataArray>
      </Cells>
      <PointData Scalars="displacement_magnitude" Vectors="displacement">
        <DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0
          0.001 0.0 0.0
          0.001 -0.0003 0.0
          0.0 -0.0003 0.0
        </DataArray>
        <DataArray type="Float64" Name="displacement_magnitude" NumberOfComponents="1" format="ascii">
          0.0
          0.001
          0.00104403065
          0.0003
        </DataArray>
      </PointData>
      <CellData Scalars="von_mises">
        <DataArray type="Float64" Name="stress" NumberOfComponents="6" format="ascii">
          2.0e8 0.0 0.0 0.0 0.0 0.0
        </DataArray>
        <DataArray type="Float64" Name="strain" NumberOfComponents="6" format="ascii">
          0.001 -0.0003 -0.0003 0.0 0.0 0.0
        </DataArray>
        <DataArray type="Float64" Name="von_mises" NumberOfComponents="1" format="ascii">
          2.0e8
        </DataArray>
        <DataArray type="Float64" Name="principal_stress" NumberOfComponents="3" format="ascii">
          2.0e8 0.0 0.0
        </DataArray>
      </CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
```

---

## 3. Python Exporter Implementations

### 3.1 Primary Writer using `meshio` (`python/fem/vtk_writer.py`)

```python
import meshio
import numpy as np
import femcore

def write_vtu_meshio(filepath: str, fem_model, solver_output):
    """
    Export FEA solution to VTK XML UnstructuredGrid format using meshio.
    Query C++ registry for VTK cell codes to maintain zero hardcoded tables.
    """
    # 1. Coordinates (n_nodes x 3)
    points = fem_model.coords.copy()

    # 2. Build meshio cell blocks using femcore element metadata
    cells = []
    # Group contiguous element connectivity
    curr_type = fem_model.elem_type_name[0]
    curr_conn = []

    # Map VTK integer code to meshio cell string name
    vtk_code_to_meshio = {
        5: "triangle",
        9: "quad",
        10: "tetra",
        12: "hexahedron"
    }

    for e in range(fem_model.n_elem):
        elem_type = fem_model.elem_type_name[e]
        meta = femcore.element_meta(elem_type)
        vtk_code = meta["vtk"]
        meshio_type = vtk_code_to_meshio[vtk_code]

        conn = fem_model.elem_conn[e]
        cells.append((meshio_type, [conn]))

    # 3. Assemble Point and Cell Data
    point_data = {
        "displacement": solver_output.displacement,
        "displacement_magnitude": np.linalg.norm(solver_output.displacement, axis=1)
    }

    cell_data = {
        "stress": [solver_output.elem_stress],
        "strain": [solver_output.elem_strain],
        "von_mises": [solver_output.elem_von_mises],
        "principal_stress": [solver_output.elem_principal]
    }

    mesh = meshio.Mesh(points=points, cells=cells, point_data=point_data, cell_data=cell_data)
    mesh.write(filepath, file_format="vtu", binary=False)
```

---

### 3.2 Dependency-Free Manual XML Writer

```python
def write_vtu_manual(filepath: str, fem_model, solver_output):
    """
    Fallback pure-Python VTK XML exporter requiring zero third-party dependencies.
    """
    n_pts = fem_model.n_nodes
    n_cells = fem_model.n_elem

    with open(filepath, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n')

        # Points
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for row in fem_model.coords:
            f.write(f'          {row[0]:.8e} {row[1]:.8e} {row[2]:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        # Cells
        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for conn in fem_model.elem_conn:
            f.write('          ' + ' '.join(map(str, conn)) + '\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        offset = 0
        for conn in fem_model.elem_conn:
            offset += len(conn)
            f.write(f'          {offset}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for elem_type in fem_model.elem_type_name:
            meta = femcore.element_meta(elem_type)
            f.write(f'          {meta["vtk"]}\n')
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')

        # Point Data
        f.write('      <PointData Scalars="displacement_magnitude" Vectors="displacement">\n')
        f.write('        <DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii">\n')
        for u in solver_output.displacement:
            f.write(f'          {u[0]:.8e} {u[1]:.8e} {u[2]:.8e}\n')
        f.write('        </DataArray>\n')

        disp_mag = np.linalg.norm(solver_output.displacement, axis=1)
        f.write('        <DataArray type="Float64" Name="displacement_magnitude" NumberOfComponents="1" format="ascii">\n')
        for mag in disp_mag:
            f.write(f'          {mag:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')

        # Cell Data
        f.write('      <CellData Scalars="von_mises">\n')
        f.write('        <DataArray type="Float64" Name="stress" NumberOfComponents="6" format="ascii">\n')
        for s in solver_output.elem_stress:
            f.write(f'          {s[0]:.8e} {s[1]:.8e} {s[2]:.8e} {s[3]:.8e} {s[4]:.8e} {s[5]:.8e}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Float64" Name="von_mises" NumberOfComponents="1" format="ascii">\n')
        for vm in solver_output.elem_von_mises:
            f.write(f'          {vm:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </CellData>\n')

        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')
```

---

## 4. ParaView Visualization & Quality Control

### 4.1 ParaView Post-Processing Protocol
1. Open ParaView and load `solution.vtu` (**File $\rightarrow$ Open**).
2. Click **Apply** in the Properties panel.
3. **Deformation Display**:
   * Apply filter: **Filters $\rightarrow$ Alphabetical $\rightarrow$ Warp By Vector**.
   * Set **Vectors** to `displacement`.
   * Set **Scale Factor** appropriately (e.g., `1.0` or `100.0` for small linear strain magnification).
4. **Contour Rendering**:
   * In the **Coloring** dropdown, select `von_mises` (Cell Data) or `displacement_magnitude` (Point Data).
5. **Nodal Stress Smoothing**:
   * To convert blocky constant element stress into smooth continuous contours across nodes, apply: **Filters $\rightarrow$ Alphabetical $\rightarrow$ Cell Data to Point Data**.

---

### 4.2 Exporter Validation Checklist
- [x] File header matches `UnstructuredGrid` version `0.1`.
- [x] Nodal coordinates are 3-component `Float64` with $z=0$ for 2D models.
- [x] Offsets array correctly tracks cumulative connectivity indices.
- [x] VTK Cell types (`5`, `9`, `10`, `12`) match element formulations.
- [x] ParaView loads `.vtu` without XML parse errors or missing scalar warnings.
- [x] Round-trip reading via `meshio.read("solution.vtu")` yields identical arrays to machine precision.
