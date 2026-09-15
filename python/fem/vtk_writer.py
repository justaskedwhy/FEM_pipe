"""Module M10: VTK XML UnstructuredGrid (.vtu) writers.

Two exporters per doc 04 section 3:
  * write_vtu_manual  - dependency-free ASCII XML writer (fallback)
  * write_vtu_meshio  - exporter via the meshio library

VTK cell type codes come from the C++ registry through meta.element_meta()
(invariant I7). The meshio glue table below maps VTK *cell codes* to meshio
cell names; it is VTK/meshio vocabulary, not element metadata, so it does not
violate invariant I2.
"""

import logging

import numpy as np

from .meta import element_meta

logger = logging.getLogger("fem.vtu")

VTK_TO_MESHIO = {
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
}


def _permuted_connectivity(model) -> list[list[int]]:
    out = []
    for t, conn in zip(model.elem_type_name, model.elem_conn):
        perm = element_meta(t)["vtk_permutation"]
        out.append([conn[i] for i in perm])
    return out


def write_vtu_manual(filepath: str, fem_model, solver_output) -> None:
    """Dependency-free ASCII VTK XML exporter (doc 04 section 3.2)."""
    n_pts = fem_model.n_nodes
    n_cells = fem_model.n_elem
    conn = _permuted_connectivity(fem_model)
    stress = np.asarray(solver_output.elem_stress)
    strain = np.asarray(solver_output.elem_strain)
    disp = np.asarray(solver_output.displacement)
    vm = np.asarray(solver_output.elem_von_mises)

    with open(filepath, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">\n')

        f.write('      <Points>\n')
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for row in fem_model.coords:
            f.write(f'          {row[0]:.8e} {row[1]:.8e} {row[2]:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        f.write('      <Cells>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for c in conn:
            f.write('          ' + ' '.join(map(str, c)) + '\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        offset = 0
        for c in conn:
            offset += len(c)
            f.write(f'          {offset}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for elem_type in fem_model.elem_type_name:
            f.write(f'          {element_meta(elem_type)["vtk"]}\n')
        f.write('        </DataArray>\n')
        f.write('      </Cells>\n')

        f.write('      <PointData Scalars="displacement_magnitude" Vectors="displacement">\n')
        f.write('        <DataArray type="Float64" Name="displacement" '
                'NumberOfComponents="3" format="ascii">\n')
        for u in disp:
            f.write(f'          {u[0]:.8e} {u[1]:.8e} {u[2]:.8e}\n')
        f.write('        </DataArray>\n')

        disp_mag = np.linalg.norm(disp, axis=1)
        f.write('        <DataArray type="Float64" Name="displacement_magnitude" '
                'NumberOfComponents="1" format="ascii">\n')
        for mag in disp_mag:
            f.write(f'          {mag:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </PointData>\n')

        f.write('      <CellData Scalars="von_mises">\n')
        f.write('        <DataArray type="Float64" Name="stress" '
                'NumberOfComponents="6" format="ascii">\n')
        for s in stress:
            f.write(f'          {s[0]:.8e} {s[1]:.8e} {s[2]:.8e} {s[3]:.8e} {s[4]:.8e} '
                    f'{s[5]:.8e}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Float64" Name="strain" '
                'NumberOfComponents="6" format="ascii">\n')
        for s in strain:
            f.write(f'          {s[0]:.8e} {s[1]:.8e} {s[2]:.8e} {s[3]:.8e} {s[4]:.8e} '
                    f'{s[5]:.8e}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Float64" Name="von_mises" '
                'NumberOfComponents="1" format="ascii">\n')
        for v in vm:
            f.write(f'          {v:.8e}\n')
        f.write('        </DataArray>\n')

        f.write('        <DataArray type="Float64" Name="principal_stress" '
                'NumberOfComponents="3" format="ascii">\n')
        for p in np.asarray(solver_output.elem_principal):
            f.write(f'          {p[0]:.8e} {p[1]:.8e} {p[2]:.8e}\n')
        f.write('        </DataArray>\n')
        f.write('      </CellData>\n')

        f.write('    </Piece>\n')
        f.write('  </UnstructuredGrid>\n')
        f.write('</VTKFile>\n')


def write_vtu_meshio(filepath: str, fem_model, solver_output) -> None:
    """Export via the meshio library (doc 04 section 3.1)."""
    import meshio

    points = np.array(fem_model.coords, dtype=np.float64)
    conn = _permuted_connectivity(fem_model)

    cells = []
    for t, c in zip(fem_model.elem_type_name, conn):
        vtk_code = element_meta(t)["vtk"]
        meshio_type = VTK_TO_MESHIO[vtk_code]
        cells.append((meshio_type, np.asarray(c, dtype=np.int64).reshape(1, -1)))

    point_data = {
        "displacement": np.asarray(solver_output.displacement, dtype=np.float64),
        "displacement_magnitude": np.linalg.norm(solver_output.displacement, axis=1),
    }
    cell_data = {
        "stress": [np.asarray(solver_output.elem_stress, dtype=np.float64)],
        "strain": [np.asarray(solver_output.elem_strain, dtype=np.float64)],
        "von_mises": [np.asarray(solver_output.elem_von_mises, dtype=np.float64)],
        "principal_stress": [np.asarray(solver_output.elem_principal, dtype=np.float64)],
    }

    mesh = meshio.Mesh(points=points, cells=cells, point_data=point_data, cell_data=cell_data)
    mesh.write(filepath, file_format="vtu", binary=False)