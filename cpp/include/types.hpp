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
// ReducedSystem — M7 output after applying prescribed displacement BCs.
// K_r * U_r = F_r on the free-free block.
// =============================================================
struct ReducedSystem {
    Eigen::SparseMatrix<double> Kr;               // Free-free block of global stiffness matrix
    Eigen::VectorXd             Fr;               // Reduced force RHS vector (F_f - K_fc * U_c)
    std::vector<int32_t>        free_dofs;        // Global DOF indices of free DOFs
    std::vector<int32_t>        constrained_dofs; // Global DOF indices of constrained DOFs
    Eigen::VectorXd             constrained_vals; // Prescribed displacement values
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