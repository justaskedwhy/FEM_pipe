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