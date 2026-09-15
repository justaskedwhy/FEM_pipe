#include "element_kernel.hpp"
#include "element_registry.hpp"
#include <stdexcept>
#include <cmath>

namespace fem {

class T4Kernel : public ElementKernel {
public:
    int dim()                        const override { return 3; }
    int nodes_per_element()          const override { return 4; }
    int dofs_per_element()           const override { return 12; }
    int n_strain_components()        const override { return 6; }
    int n_gauss_points()             const override { return 1; }
    int vtk_cell_type()              const override { return 10; /* VTK_TETRA */ }

    ElementMatrices compute(const ElementContext& ctx) const override {
        Eigen::Matrix<double, 4, 3> dN_dxi;
        dN_dxi << -1.0, -1.0, -1.0,
                   1.0,  0.0,  0.0,
                   0.0,  1.0,  0.0,
                   0.0,  0.0,  1.0;

        Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 4; ++a) {
            J += dN_dxi.row(a).transpose() * ctx.coords.row(a);
        }

        const double detJ = J.determinant();
        if (detJ <= 0.0) {
            throw std::runtime_error("T4Kernel: Degenerate or inverted tetrahedron element.");
        }
        const double V = detJ / 6.0;

        Eigen::Matrix3d Jinv = J.inverse();
        Eigen::Matrix<double, 4, 3> dN_dx = dN_dxi * Jinv;

        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 12);
        for (int a = 0; a < 4; ++a) {
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

        Eigen::MatrixXd D = build_D_3d(ctx.E, ctx.nu);
        Eigen::MatrixXd Ke = V * B.transpose() * D * B;

        ElementMatrices out;
        out.Ke = Ke;
        out.B = B;
        out.D = D;
        out.gauss_points = Eigen::MatrixXd::Zero(1, 3);
        out.gauss_points(0, 0) = 0.25;
        out.gauss_points(0, 1) = 0.25;
        out.gauss_points(0, 2) = 0.25;
        out.gauss_weights = Eigen::VectorXd::Constant(1, 1.0);
        out.detJ = Eigen::VectorXd::Constant(1, detJ);
        return out;
    }

    Eigen::VectorXd pressure_force(
        const ElementContext& ctx, int face_id, double p) const override {
        if (face_id < 1 || face_id > n_faces()) {
            throw std::runtime_error("T4Kernel: face_id out of bounds");
        }
        std::vector<int> fn = face_nodes(face_id);
        if (fn.size() != 3) {
            throw std::runtime_error("T4Kernel: internal face topology error");
        }

        const Eigen::Vector3d a = ctx.coords.row(fn[0]);
        const Eigen::Vector3d b = ctx.coords.row(fn[1]);
        const Eigen::Vector3d c = ctx.coords.row(fn[2]);

        Eigen::Vector3d n = (b - a).cross(c - a);
        const double twoA = n.norm();
        if (twoA < 1e-14) {
            throw std::runtime_error("T4Kernel: Degenerate face.");
        }
        n /= twoA;

        // Orient the normal outward from the element centroid.
        Eigen::Vector3d elem_c = Eigen::Vector3d::Zero();
        for (int i = 0; i < 4; ++i) elem_c += ctx.coords.row(i);
        elem_c /= 4.0;
        const Eigen::Vector3d face_c = (a + b + c) / 3.0;
        if (n.dot(face_c - elem_c) < 0.0) n = -n;

        const Eigen::Vector3d t = -p * n;
        const double A = 0.5 * twoA;
        const double share = A / 3.0;

        Eigen::VectorXd fe = Eigen::VectorXd::Zero(12);
        for (int i = 0; i < 3; ++i) {
            fe(3*fn[i] + 0) = t[0] * share;
            fe(3*fn[i] + 1) = t[1] * share;
            fe(3*fn[i] + 2) = t[2] * share;
        }
        return fe;
    }

    Eigen::MatrixXd B_at(
        const ElementContext& ctx,
        const Eigen::VectorXd& /*xi*/) const override {
        return compute(ctx).B;
    }

    int n_faces() const override { return 4; }

    std::vector<int> face_nodes(int face_id) const override {
        static const std::vector<std::vector<int>> f = {
            {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
        if (face_id < 1 || face_id > n_faces()) {
            throw std::runtime_error("T4Kernel: Invalid face_id");
        }
        return f[face_id - 1];
    }

    std::vector<int> vtk_permutation() const override { return {0, 1, 2, 3}; }

private:
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

REGISTER_ELEMENT("C3D4", fem::T4Kernel);