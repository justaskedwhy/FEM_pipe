#include "element_kernel.hpp"
#include "element_registry.hpp"
#include <stdexcept>
#include <cmath>

namespace fem {

class H8Kernel : public ElementKernel {
public:
    int dim()                        const override { return 3; }
    int nodes_per_element()          const override { return 8; }
    int dofs_per_element()           const override { return 24; }
    int n_strain_components()        const override { return 6; }
    int n_gauss_points()             const override { return 8; }
    int vtk_cell_type()              const override { return 12; /* VTK_HEXAHEDRON */ }

    ElementMatrices compute(const ElementContext& ctx) const override {
        static const double gp = 1.0 / std::sqrt(3.0);
        static const double gauss[8][3] = {
            {-gp, -gp, -gp}, { gp, -gp, -gp}, { gp,  gp, -gp}, {-gp,  gp, -gp},
            {-gp, -gp,  gp}, { gp, -gp,  gp}, { gp,  gp,  gp}, {-gp,  gp,  gp}};
        static const double node_ref[8][3] = {
            {-1, -1, -1}, { 1, -1, -1}, { 1,  1, -1}, {-1,  1, -1},
            {-1, -1,  1}, { 1, -1,  1}, { 1,  1,  1}, {-1,  1,  1}};

        Eigen::MatrixXd D = build_D_3d(ctx.E, ctx.nu);
        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(24, 24);
        Eigen::MatrixXd B_stacked(48, 24);
        Eigen::VectorXd detJ_vec(8);
        Eigen::MatrixXd gp_coords(8, 3);

        for (int g = 0; g < 8; ++g) {
            const double xi = gauss[g][0], eta = gauss[g][1], zeta = gauss[g][2];
            gp_coords(g, 0) = xi;
            gp_coords(g, 1) = eta;
            gp_coords(g, 2) = zeta;

            Eigen::Matrix<double, 8, 3> dN_dxi;
            for (int a = 0; a < 8; ++a) {
                dN_dxi(a, 0) = 0.125 * node_ref[a][0]
                    * (1.0 + node_ref[a][1] * eta) * (1.0 + node_ref[a][2] * zeta);
                dN_dxi(a, 1) = 0.125 * node_ref[a][1]
                    * (1.0 + node_ref[a][0] * xi) * (1.0 + node_ref[a][2] * zeta);
                dN_dxi(a, 2) = 0.125 * node_ref[a][2]
                    * (1.0 + node_ref[a][0] * xi) * (1.0 + node_ref[a][1] * eta);
            }

            Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
            for (int a = 0; a < 8; ++a) {
                J += dN_dxi.row(a).transpose() * ctx.coords.row(a);
            }

            const double detJ = J.determinant();
            if (detJ <= 0.0) {
                throw std::runtime_error("H8Kernel: Degenerate Jacobian at Gauss point.");
            }
            detJ_vec[g] = detJ;

            Eigen::Matrix3d Jinv = J.inverse();
            Eigen::Matrix<double, 8, 3> dN_dx = dN_dxi * Jinv;

            Eigen::MatrixXd B = build_B(dN_dx);
            B_stacked.block<6, 24>(6 * g, 0) = B;
            Ke += B.transpose() * D * B * detJ;
        }

        ElementMatrices out;
        out.Ke = Ke;
        out.B = B_stacked;
        out.D = D;
        out.gauss_points = gp_coords;
        out.gauss_weights = Eigen::VectorXd::Constant(8, 1.0);
        out.detJ = detJ_vec;
        return out;
    }

    Eigen::VectorXd pressure_force(
        const ElementContext& ctx, int face_id, double p) const override {
        if (face_id < 1 || face_id > n_faces()) {
            throw std::runtime_error("H8Kernel: face_id out of bounds");
        }
        std::vector<int> fn = face_nodes(face_id);
        if (fn.size() != 4) {
            throw std::runtime_error("H8Kernel: internal face topology error");
        }

        static const double gp = 1.0 / std::sqrt(3.0);
        static const double gauss[4][2] = {{-gp,-gp}, {gp,-gp}, {gp,gp}, {-gp,gp}};
        static const double ref[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};

        Eigen::Vector3d elem_c = Eigen::Vector3d::Zero();
        for (int a = 0; a < 8; ++a) elem_c += ctx.coords.row(a);
        elem_c /= 8.0;

        Eigen::VectorXd fe = Eigen::VectorXd::Zero(24);
        for (int g = 0; g < 4; ++g) {
            const double xi = gauss[g][0], eta = gauss[g][1];
            double N[4], dNdxi[4], dNdeta[4];
            Eigen::Vector3d r = Eigen::Vector3d::Zero();
            Eigen::Vector3d r_xi = Eigen::Vector3d::Zero();
            Eigen::Vector3d r_eta = Eigen::Vector3d::Zero();
            for (int i = 0; i < 4; ++i) {
                N[i] = 0.25 * (1.0 + xi * ref[i][0]) * (1.0 + eta * ref[i][1]);
                dNdxi[i] = 0.25 * ref[i][0] * (1.0 + eta * ref[i][1]);
                dNdeta[i] = 0.25 * ref[i][1] * (1.0 + xi * ref[i][0]);
                r     += N[i] * ctx.coords.row(fn[i]);
                r_xi  += dNdxi[i] * ctx.coords.row(fn[i]);
                r_eta += dNdeta[i] * ctx.coords.row(fn[i]);
            }

            Eigen::Vector3d n = r_xi.cross(r_eta);
            const double detJ = n.norm();
            if (detJ < 1e-14) {
                throw std::runtime_error("H8Kernel: Degenerate face.");
            }
            n /= detJ;
            if (n.dot(r - elem_c) < 0.0) n = -n;

            const Eigen::Vector3d t = -p * n;
            for (int i = 0; i < 4; ++i) {
                fe(3*fn[i] + 0) += N[i] * t[0] * detJ;
                fe(3*fn[i] + 1) += N[i] * t[1] * detJ;
                fe(3*fn[i] + 2) += N[i] * t[2] * detJ;
            }
        }
        return fe;
    }

    Eigen::MatrixXd B_at(
        const ElementContext& ctx,
        const Eigen::VectorXd& xi) const override {
        static const double node_ref[8][3] = {
            {-1, -1, -1}, { 1, -1, -1}, { 1,  1, -1}, {-1,  1, -1},
            {-1, -1,  1}, { 1, -1,  1}, { 1,  1,  1}, {-1,  1,  1}};

        const double xi_v = xi(0), eta = xi(1), zeta = xi(2);

        Eigen::Matrix<double, 8, 3> dN_dxi;
        for (int a = 0; a < 8; ++a) {
            dN_dxi(a, 0) = 0.125 * node_ref[a][0]
                * (1.0 + node_ref[a][1] * eta) * (1.0 + node_ref[a][2] * zeta);
            dN_dxi(a, 1) = 0.125 * node_ref[a][1]
                * (1.0 + node_ref[a][0] * xi_v) * (1.0 + node_ref[a][2] * zeta);
            dN_dxi(a, 2) = 0.125 * node_ref[a][2]
                * (1.0 + node_ref[a][0] * xi_v) * (1.0 + node_ref[a][1] * eta);
        }

        Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 8; ++a) {
            J += dN_dxi.row(a).transpose() * ctx.coords.row(a);
        }

        const double detJ = J.determinant();
        if (detJ <= 0.0) {
            throw std::runtime_error("H8Kernel: Degenerate Jacobian at evaluation point.");
        }

        Eigen::Matrix3d Jinv = J.inverse();
        Eigen::Matrix<double, 8, 3> dN_dx = dN_dxi * Jinv;
        return build_B(dN_dx);
    }

    int n_faces() const override { return 6; }

    std::vector<int> face_nodes(int face_id) const override {
        static const std::vector<std::vector<int>> f = {
            {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
            {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
        if (face_id < 1 || face_id > n_faces()) {
            throw std::runtime_error("H8Kernel: Invalid face_id");
        }
        return f[face_id - 1];
    }

    std::vector<int> vtk_permutation() const override { return {0, 1, 2, 3, 4, 5, 6, 7}; }

private:
    static Eigen::MatrixXd build_B(const Eigen::Matrix<double, 8, 3>& dN_dx) {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 24);
        for (int a = 0; a < 8; ++a) {
            const double dNx = dN_dx(a, 0);
            const double dNy = dN_dx(a, 1);
            const double dNz = dN_dx(a, 2);
            B(0, 3*a + 0) = dNx;
            B(1, 3*a + 1) = dNy;
            B(2, 3*a + 2) = dNz;
            B(3, 3*a + 0) = dNy;
            B(3, 3*a + 1) = dNx;
            B(4, 3*a + 1) = dNz;
            B(4, 3*a + 2) = dNy;
            B(5, 3*a + 0) = dNz;
            B(5, 3*a + 2) = dNx;
        }
        return B;
    }

    static Eigen::MatrixXd build_D_3d(double E, double nu) {
        Eigen::MatrixXd D(6, 6);
        D.setZero();
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double a = 1.0 - nu;
        const double b = (1.0 - 2.0 * nu) / 2.0;
        D << a, nu, nu, 0.0, 0.0, 0.0,
             nu, a, nu, 0.0, 0.0, 0.0,
             nu, nu, a, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, b, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, b, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, b;
        D *= c;
        return D;
    }
};

} // namespace fem

REGISTER_ELEMENT("C3D8", fem::H8Kernel);
REGISTER_ELEMENT("C3D8R", fem::H8Kernel);