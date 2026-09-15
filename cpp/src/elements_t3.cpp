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