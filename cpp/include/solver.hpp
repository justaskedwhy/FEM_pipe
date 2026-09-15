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