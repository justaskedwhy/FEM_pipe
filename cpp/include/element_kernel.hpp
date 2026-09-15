// cpp/include/element_kernel.hpp
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <stdexcept>

namespace fem {

// =============================================================
// ElementContext — Uniform element-local driver input.
// coords is (npe x 3); 2D elements store z = 0.
// =============================================================
struct ElementContext {
    Eigen::Matrix<double, Eigen::Dynamic, 3> coords; // (npe x 3) node physical coordinates
    double E;             // Young's modulus (Pa)
    double nu;            // Poisson's ratio
    double thickness;     // Out-of-plane thickness (2D only; 1.0 for 3D)
    int    plane_mode;    // 0 = plane stress, 1 = plane strain, 2 = 3D
};

// =============================================================
// ElementMatrices — Uniform element-local driver output.
// B is stacked row-wise across Gauss points:
// rows [g*n_strain, (g+1)*n_strain - 1] hold Gauss point g's B.
// =============================================================
struct ElementMatrices {
    Eigen::MatrixXd Ke;              // (ndof_e x ndof_e) dense element stiffness matrix
    Eigen::MatrixXd B;               // ((n_gauss * n_strain) x ndof_e) stacked strain-displacement matrix
    Eigen::MatrixXd D;               // (n_strain x n_strain) elastic constitutive matrix
    Eigen::MatrixXd gauss_points;    // (n_gauss x 3) natural coordinates of Gauss points
    Eigen::VectorXd gauss_weights;   // (n_gauss) Gauss integration weights
    Eigen::VectorXd detJ;            // (n_gauss) Jacobian determinants at Gauss points
};

// =============================================================
// ElementKernel — Pluggable element formulation interface.
// =============================================================
class ElementKernel {
public:
    virtual ~ElementKernel() = default;

    // ---------------------------------------------------------------
    // 1. Metadata Methods (Cheap, Constant Queries)
    // ---------------------------------------------------------------
    // ABAQUS canonical identifier of this element (e.g. "CPS4", "CPE4").
    // Injected by the registry at construction time (see REGISTER_ELEMENT).
    // T3 and Q4 kernels are each registered under two ABAQUS names; the
    // overloads share one formulation class but report distinct name().
    virtual std::string name() const { return registry_name_; }
    virtual int dim()                          const = 0; // Spatial dimension: 2 or 3
    virtual int nodes_per_element()            const = 0; // npe: 3, 4, or 8
    virtual int dofs_per_element()             const = 0; // ndof_e = npe * dim
    virtual int n_strain_components()          const = 0; // 3 (2D) or 6 (3D)
    virtual int n_gauss_points()               const = 0; // Number of integration points
    virtual int vtk_cell_type()                const = 0; // VTK XML cell type code

    void set_registry_name(std::string n) { registry_name_ = std::move(n); }

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

private:
    std::string registry_name_;
};

} // namespace fem