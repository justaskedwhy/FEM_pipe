#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "solver.hpp"
#include "element_registry.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

py::dict element_meta(const std::string& name) {
    const fem::ElementRegistry& registry = fem::ElementRegistry::instance();
    if (!registry.has(name)) {
        throw py::key_error("[femcore] Unregistered element: " + name);
    }
    auto kernel = registry.create(name);
    py::dict meta;
    meta["vtk"] = kernel->vtk_cell_type();
    meta["npe"] = kernel->nodes_per_element();
    meta["dim"] = kernel->dim();
    meta["dofs"] = kernel->dofs_per_element();
    meta["strain_components"] = kernel->n_strain_components();
    meta["n_gauss"] = kernel->n_gauss_points();
    meta["n_faces"] = kernel->n_faces();
    meta["faces"] = kernel->face_nodes(1);
    meta["vtk_permutation"] = kernel->vtk_permutation();
    return meta;
}

py::list registered_elements() {
    std::vector<std::string> names =
        fem::ElementRegistry::instance().registered_names();
    std::sort(names.begin(), names.end());
    py::list out;
    for (const auto& n : names) out.append(n);
    return out;
}

py::dict registered_metadata() {
    std::vector<std::string> names =
        fem::ElementRegistry::instance().registered_names();
    std::sort(names.begin(), names.end());
    py::dict out;
    for (const auto& n : names) out[n.c_str()] = element_meta(n);
    return out;
}

Eigen::MatrixXd test_element_b_at(const std::string& name,
                                  const Eigen::MatrixXd& coords,
                                  double xi, double eta, double zeta,
                                  double E, double nu,
                                  double thickness, int plane_mode) {
    const fem::ElementRegistry& registry = fem::ElementRegistry::instance();
    auto kernel = registry.create(name);
    fem::ElementContext ctx;
    ctx.coords = coords;
    ctx.E = E;
    ctx.nu = nu;
    ctx.thickness = thickness;
    ctx.plane_mode = plane_mode;
    Eigen::VectorXd xi_v(3);
    xi_v << xi, eta, zeta;
    return kernel->B_at(ctx, xi_v);
}

Eigen::VectorXd test_element_pressure(const std::string& name,
                                      const Eigen::MatrixXd& coords,
                                      int face_id, double p,
                                      double thickness) {
    const fem::ElementRegistry& registry = fem::ElementRegistry::instance();
    auto kernel = registry.create(name);
    fem::ElementContext ctx;
    ctx.coords = coords;
    ctx.E = 0.0;
    ctx.nu = 0.0;
    ctx.thickness = thickness;
    ctx.plane_mode = 0;
    return kernel->pressure_force(ctx, face_id, p);
}

Eigen::MatrixXd test_element_matrix(const std::string& name,
                                    const Eigen::MatrixXd& coords,
                                    double E, double nu,
                                    double thickness, int plane_mode) {
    const fem::ElementRegistry& registry = fem::ElementRegistry::instance();
    auto kernel = registry.create(name);
    fem::ElementContext ctx;
    ctx.coords = coords;
    ctx.E = E;
    ctx.nu = nu;
    ctx.thickness = thickness;
    ctx.plane_mode = plane_mode;
    return kernel->compute(ctx).Ke;
}

} // namespace

PYBIND11_MODULE(femcore, m) {
    m.doc() = "femcore: C++ numerical core for the hybrid Python-C++ FE solver";

    py::class_<fem::SolverInput>(m, "SolverInput")
        .def(py::init<>())
        .def_readwrite("coords", &fem::SolverInput::coords)
        .def_readwrite("elem_type_name", &fem::SolverInput::elem_type_name)
        .def_readwrite("elem_conn_flat", &fem::SolverInput::elem_conn_flat)
        .def_readwrite("elem_conn_offsets", &fem::SolverInput::elem_conn_offsets)
        .def_readwrite("elem_mat", &fem::SolverInput::elem_mat)
        .def_readwrite("mat_E", &fem::SolverInput::mat_E)
        .def_readwrite("mat_nu", &fem::SolverInput::mat_nu)
        .def_readwrite("mat_thickness", &fem::SolverInput::mat_thickness)
        .def_readwrite("mat_plane_mode", &fem::SolverInput::mat_plane_mode)
        .def_readwrite("dof_map", &fem::SolverInput::dof_map)
        .def_readwrite("n_dofs_total", &fem::SolverInput::n_dofs_total)
        .def_readwrite("n_dofs_free", &fem::SolverInput::n_dofs_free)
        .def_readwrite("prescribed_dofs", &fem::SolverInput::prescribed_dofs)
        .def_readwrite("prescribed_vals", &fem::SolverInput::prescribed_vals)
        .def_readwrite("point_load_dofs", &fem::SolverInput::point_load_dofs)
        .def_readwrite("point_load_vals", &fem::SolverInput::point_load_vals)
        .def_readwrite("pressure_elem", &fem::SolverInput::pressure_elem)
        .def_readwrite("pressure_face", &fem::SolverInput::pressure_face)
        .def_readwrite("pressure_val", &fem::SolverInput::pressure_val)
        .def_readwrite("n_nodes", &fem::SolverInput::n_nodes)
        .def_readwrite("n_elem", &fem::SolverInput::n_elem)
        .def_readwrite("n_mat", &fem::SolverInput::n_mat)
        .def_readwrite("dimension", &fem::SolverInput::dimension);

    py::class_<fem::SolverOutput>(m, "SolverOutput")
        .def(py::init<>())
        .def_readwrite("displacement", &fem::SolverOutput::displacement)
        .def_readwrite("reaction", &fem::SolverOutput::reaction)
        .def_readwrite("elem_stress", &fem::SolverOutput::elem_stress)
        .def_readwrite("elem_strain", &fem::SolverOutput::elem_strain)
        .def_readwrite("elem_von_mises", &fem::SolverOutput::elem_von_mises)
        .def_readwrite("elem_principal", &fem::SolverOutput::elem_principal)
        .def_readwrite("converged", &fem::SolverOutput::converged)
        .def_readwrite("residual_norm", &fem::SolverOutput::residual_norm);

    m.def("fem_solve", &fem::fem_solve,
          "Solve the linear elastostatic system described by a SolverInput.",
          py::arg("solver_input"));

    m.def("element_meta", &element_meta,
          "Return a dict of element metadata for an ABAQUS element name.",
          py::arg("name"));

    m.def("registered_elements", &registered_elements,
          "Return the sorted list of registered ABAQUS element names.");

    m.def("registered_metadata", &registered_metadata,
          "Return {ABAQUS name: element_meta} for every registered element.");

    m.def("test_element_matrix", &test_element_matrix,
          "Compute a raw element stiffness matrix for a single element (test hook).",
          py::arg("name"), py::arg("coords"), py::arg("E"), py::arg("nu"),
          py::arg("thickness"), py::arg("plane_mode") = 0);

    m.def("test_element_b_at", &test_element_b_at,
          "Evaluate the strain-displacement matrix B at natural coords (test hook).",
          py::arg("name"), py::arg("coords"), py::arg("xi"), py::arg("eta"),
          py::arg("zeta") = 0.0, py::arg("E") = 0.0, py::arg("nu") = 0.0,
          py::arg("thickness") = 1.0, py::arg("plane_mode") = 0);

    m.def("test_element_pressure", &test_element_pressure,
          "Compute the equivalent nodal pressure load vector for one face (test hook).",
          py::arg("name"), py::arg("coords"), py::arg("face_id"), py::arg("p"),
          py::arg("thickness") = 1.0);
}