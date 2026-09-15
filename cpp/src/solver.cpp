#include "solver.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace fem {

namespace {

void validate_input(const SolverInput& in) {
    if (in.coords.rows() != in.n_nodes || in.coords.cols() != 3) {
        throw FemError("coords must have shape (n_nodes, 3)");
    }
    if (in.elem_conn_offsets.size() != static_cast<size_t>(in.n_elem) + 1) {
        throw FemError("elem_conn_offsets must have size n_elem + 1");
    }
    if (in.elem_conn_offsets.empty() || in.elem_conn_offsets.front() != 0 ||
        in.elem_conn_offsets.back() != static_cast<int32_t>(in.elem_conn_flat.size())) {
        throw FemError("elem_conn_offsets[0] must be 0 and offsets[n_elem] == flat.size()");
    }
    if (in.elem_type_name.size() != static_cast<size_t>(in.n_elem)) {
        throw FemError("elem_type_name must have size n_elem");
    }
    if (in.elem_mat.size() != static_cast<size_t>(in.n_elem)) {
        throw FemError("elem_mat must have size n_elem");
    }
    if (in.mat_E.size() != in.n_mat || in.mat_nu.size() != in.n_mat ||
        in.mat_thickness.size() != in.n_mat ||
        in.mat_plane_mode.size() != static_cast<size_t>(in.n_mat)) {
        throw FemError("material arrays must have size n_mat");
    }
    if (in.dof_map.rows() != in.n_nodes || in.dof_map.cols() != 3) {
        throw FemError("dof_map must have shape (n_nodes, 3)");
    }
    if (in.prescribed_dofs.size() != in.prescribed_vals.size()) {
        throw FemError("prescribed_dofs/values length mismatch");
    }
    if (in.point_load_dofs.size() != in.point_load_vals.size()) {
        throw FemError("point_load_dofs/values length mismatch");
    }
    if (in.pressure_elem.size() != in.pressure_val.size() ||
        in.pressure_face.size() != in.pressure_val.size()) {
        throw FemError("pressure arrays length mismatch");
    }

    const ElementRegistry& registry = ElementRegistry::instance();
    for (int32_t e = 0; e < in.n_elem; ++e) {
        if (!registry.has(in.elem_type_name[e])) {
            throw FemError("Unregistered element type: " + in.elem_type_name[e]);
        }
    }
    for (int32_t e = 0; e < in.n_elem; ++e) {
        if (in.elem_mat[e] < 0 || in.elem_mat[e] >= in.n_mat) {
            throw FemError("elem_mat index out of bounds for element " + std::to_string(e));
        }
    }
    for (int32_t m = 0; m < in.n_mat; ++m) {
        if (!(in.mat_E[m] > 0.0)) throw FemError("Young's modulus must be positive");
        if (!(in.mat_nu[m] > -1.0 && in.mat_nu[m] < 0.5)) {
            throw FemError("Poisson's ratio must satisfy -1 < nu < 0.5");
        }
    }
}

int32_t dim_mismatch_guard(const ElementKernel& kernel, const SolverInput& in) {
    if (kernel.dim() != in.dimension) {
        throw FemError("Element dimension mismatch: " + kernel.name());
    }
    return kernel.dim();
}

} // namespace

ElementContext build_context(const SolverInput& in, int32_t e) {
    const auto& registry = ElementRegistry::instance();
    auto kernel = registry.create(in.elem_type_name[e]);
    dim_mismatch_guard(*kernel, in);

    const int32_t npe = kernel->nodes_per_element();
    const int32_t start = in.elem_conn_offsets[e];
    ElementContext ctx;
    ctx.coords.resize(npe, 3);
    for (int32_t a = 0; a < npe; ++a) {
        const int32_t n = in.elem_conn_flat[start + a];
        if (n < 0 || n >= in.n_nodes) {
            throw FemError("Element " + std::to_string(e) + " references invalid node " + std::to_string(n));
        }
        ctx.coords.row(a) = in.coords.row(n);
    }
    const int32_t m = in.elem_mat[e];
    ctx.E = in.mat_E[m];
    ctx.nu = in.mat_nu[m];
    ctx.thickness = in.mat_thickness[m];
    ctx.plane_mode = in.mat_plane_mode[m];
    return ctx;
}

ElementDofMap build_element_dof_map(const SolverInput& in, int32_t e) {
    const auto& registry = ElementRegistry::instance();
    auto kernel = registry.create(in.elem_type_name[e]);
    const int32_t dim = dim_mismatch_guard(*kernel, in);

    const int32_t npe = kernel->nodes_per_element();
    const int32_t start = in.elem_conn_offsets[e];
    ElementDofMap map;
    map.global_dofs.reserve(static_cast<size_t>(npe * dim));
    for (int32_t a = 0; a < npe; ++a) {
        const int32_t n = in.elem_conn_flat[start + a];
        for (int32_t d = 0; d < dim; ++d) {
            const int32_t dof = in.dof_map(n, d);
            if (dof < 0) {
                throw FemError("Inactive DOF in element " + std::to_string(e));
            }
            map.global_dofs.push_back(dof);
        }
    }
    return map;
}

Eigen::SparseMatrix<double> assemble_stiffness(const SolverInput& in, const ElementRegistry& registry) {
    Eigen::SparseMatrix<double> K(in.n_dofs_total, in.n_dofs_total);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(in.n_elem) * 64);

    for (int32_t e = 0; e < in.n_elem; ++e) {
        auto kernel = registry.create(in.elem_type_name[e]);
        const int32_t ndof = kernel->dofs_per_element();
        const ElementContext ctx = build_context(in, e);
        const ElementMatrices emat = kernel->compute(ctx);
        const ElementDofMap dofmap = build_element_dof_map(in, e);
        for (int32_t i = 0; i < ndof; ++i) {
            for (int32_t j = 0; j < ndof; ++j) {
                triplets.emplace_back(dofmap.global_dofs[i], dofmap.global_dofs[j], emat.Ke(i, j));
            }
        }
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();
    return K;
}

Eigen::VectorXd assemble_forces(const SolverInput& in, const ElementRegistry& registry) {
    Eigen::VectorXd F = Eigen::VectorXd::Zero(in.n_dofs_total);

    for (size_t i = 0; i < in.point_load_dofs.size(); ++i) {
        const int32_t dof = in.point_load_dofs[i];
        if (dof < 0 || dof >= in.n_dofs_total) {
            throw FemError("Point load DOF out of bounds");
        }
        F[dof] += in.point_load_vals[i];
    }

    for (size_t i = 0; i < in.pressure_elem.size(); ++i) {
        const int32_t e = in.pressure_elem[i];
        if (e < 0 || e >= in.n_elem) {
            throw FemError("Pressure load references invalid element");
        }
        auto kernel = registry.create(in.elem_type_name[e]);
        const ElementContext ctx = build_context(in, e);
        const Eigen::VectorXd fe =
            kernel->pressure_force(ctx, in.pressure_face[i], in.pressure_val[i]);
        const ElementDofMap dofmap = build_element_dof_map(in, e);
        for (int32_t j = 0; j < fe.size(); ++j) {
            F[dofmap.global_dofs[j]] += fe[j];
        }
    }

    return F;
}

ReducedSystem apply_boundary_conditions(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd& F,
    const SolverInput& in) {

    std::vector<int32_t> is_constrained(in.n_dofs_total, 0);
    for (size_t i = 0; i < in.prescribed_dofs.size(); ++i) {
        const int32_t dof = in.prescribed_dofs[i];
        if (dof < 0 || dof >= in.n_dofs_total) {
            throw FemError("Prescribed DOF out of bounds");
        }
        is_constrained[dof] = 1;
    }

    ReducedSystem red;
    red.free_dofs.reserve(in.n_dofs_free);
    for (int32_t d = 0; d < in.n_dofs_total; ++d) {
        if (!is_constrained[d]) red.free_dofs.push_back(d);
    }
    for (int32_t d = 0; d < in.n_dofs_total; ++d) {
        if (is_constrained[d]) red.constrained_dofs.push_back(d);
    }

    red.constrained_vals.resize(red.constrained_dofs.size());
    std::vector<int32_t> pos(in.n_dofs_total, -1);
    for (size_t i = 0; i < in.prescribed_dofs.size(); ++i) {
        const int32_t dof = in.prescribed_dofs[i];
        if (dof < 0 || dof >= in.n_dofs_total) {
            throw FemError("Prescribed DOF out of bounds");
        }
        pos[dof] = static_cast<int32_t>(i);
    }
    for (size_t i = 0; i < red.constrained_dofs.size(); ++i) {
        const int32_t dof = red.constrained_dofs[i];
        red.constrained_vals[i] = in.prescribed_vals[pos[dof]];
    }

    const int32_t nf = static_cast<int32_t>(red.free_dofs.size());
    const int32_t nc = static_cast<int32_t>(red.constrained_dofs.size());

    Eigen::SparseMatrix<double> Kff(nf, nf);
    std::vector<Eigen::Triplet<double>> ff;
    ff.reserve(static_cast<size_t>(nf) * 16);
    const int32_t* outer = K.outerIndexPtr();
    const int32_t* inner = K.innerIndexPtr();
    const double* vals = K.valuePtr();
    for (int32_t i = 0; i < nf; ++i) {
        const int32_t row = red.free_dofs[i];
        for (Eigen::Index k = outer[row]; k < outer[row + 1]; ++k) {
            const int32_t col = inner[k];
            for (int32_t j = 0; j < nf; ++j) {
                if (red.free_dofs[j] == col) ff.emplace_back(i, j, vals[k]);
            }
        }
    }
    Kff.setFromTriplets(ff.begin(), ff.end());
    Kff.makeCompressed();
    red.Kr = Kff;

    Eigen::VectorXd Uc = red.constrained_vals;
    Eigen::VectorXd Kfc_Uc = Eigen::VectorXd::Zero(nf);
    for (int32_t i = 0; i < nf; ++i) {
        const int32_t row = red.free_dofs[i];
        for (Eigen::Index k = outer[row]; k < outer[row + 1]; ++k) {
            const int32_t col = inner[k];
            for (int32_t j = 0; j < nc; ++j) {
                if (red.constrained_dofs[j] == col) Kfc_Uc[i] += vals[k] * Uc[j];
            }
        }
    }

    red.Fr = Eigen::VectorXd(nf);
    for (int32_t i = 0; i < nf; ++i) {
        red.Fr[i] = F[red.free_dofs[i]] - Kfc_Uc[i];
    }
    return red;
}

Eigen::VectorXd solve_linear_system(const ReducedSystem& reduced) {
    const int32_t n = static_cast<int32_t>(reduced.Kr.rows());
    if (n == 0) {
        return Eigen::VectorXd(0);
    }
    Eigen::VectorXd Ur(n);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(reduced.Kr);
    if (ldlt.info() == Eigen::Success) {
        Ur = ldlt.solve(reduced.Fr);
        if (ldlt.info() == Eigen::Success) return Ur;
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.compute(reduced.Kr);
    if (lu.info() == Eigen::Success) {
        Ur = lu.solve(reduced.Fr);
        if (lu.info() == Eigen::Success) return Ur;
    }

    throw FemError("Linear solve failed (LDLT and LU)");
}

Eigen::VectorXd reconstruct_full_displacement(
    const Eigen::VectorXd& Ur,
    const ReducedSystem& reduced,
    int32_t n_dofs_total) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero(n_dofs_total);
    for (size_t i = 0; i < reduced.free_dofs.size(); ++i) {
        U[reduced.free_dofs[i]] = Ur[static_cast<Eigen::Index>(i)];
    }
    for (size_t i = 0; i < reduced.constrained_dofs.size(); ++i) {
        U[reduced.constrained_dofs[i]] = reduced.constrained_vals[i];
    }
    return U;
}

namespace {

void expand_2d_stress_strain(
    const Eigen::VectorXd& s2, const Eigen::VectorXd& e2,
    int plane_mode, double nu,
    Eigen::VectorXd& s6, Eigen::VectorXd& e6) {
    s6.setZero();
    e6.setZero();
    s6[0] = s2[0];
    s6[1] = s2[1];
    s6[3] = s2[2];
    e6[0] = e2[0];
    e6[1] = e2[1];
    e6[3] = e2[2];
    if (plane_mode == 0) {
        e6[2] = -nu / (1.0 - nu) * (e2[0] + e2[1]);
    } else {
        s6[2] = nu * (s2[0] + s2[1]);
    }
}

double von_mises(const Eigen::VectorXd& s) {
    const double sxx = s[0], syy = s[1], szz = s[2];
    const double sxy = s[3], syz = s[4], szx = s[5];
    const double a = (sxx - syy) * (sxx - syy);
    const double b = (syy - szz) * (syy - szz);
    const double c = (szz - sxx) * (szz - sxx);
    return std::sqrt(0.5 * (a + b + c) + 3.0 * (sxy * sxy + syz * syz + szx * szx));
}

Eigen::Vector3d principal_stresses(const Eigen::VectorXd& s) {
    Eigen::Matrix3d T;
    T << s[0], s[3], s[5],
         s[3], s[1], s[4],
         s[5], s[4], s[2];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(T);
    Eigen::Vector3d ev = solver.eigenvalues();
    std::sort(ev.data(), ev.data() + 3, std::greater<double>());
    return ev;
}

} // namespace

void compute_element_results(
    const SolverInput& in,
    const Eigen::VectorXd& U,
    const ElementRegistry& registry,
    SolverOutput& out) {

    out.elem_stress.resize(in.n_elem, 6);
    out.elem_strain.resize(in.n_elem, 6);
    out.elem_von_mises.resize(in.n_elem);
    out.elem_principal.resize(in.n_elem, 3);

    for (int32_t e = 0; e < in.n_elem; ++e) {
        auto kernel = registry.create(in.elem_type_name[e]);
        const int32_t dim = kernel->dim();
        const int32_t npe = kernel->nodes_per_element();
        const int32_t ndof = kernel->dofs_per_element();
        const int32_t nstr = kernel->n_strain_components();
        const ElementContext ctx = build_context(in, e);
        const ElementMatrices emat = kernel->compute(ctx);
        const ElementDofMap dofmap = build_element_dof_map(in, e);

        Eigen::VectorXd Ue(ndof);
        for (int32_t i = 0; i < ndof; ++i) Ue[i] = U[dofmap.global_dofs[i]];

        const int32_t ng = emat.B.rows() / nstr;
        Eigen::VectorXd e_bar = Eigen::VectorXd::Zero(nstr);
        for (int32_t g = 0; g < ng; ++g) {
            const Eigen::VectorXd eg = emat.B.block(nstr * g, 0, nstr, ndof) * Ue;
            e_bar += eg;
        }
        e_bar /= static_cast<double>(ng);

        const Eigen::VectorXd s_bar = emat.D * e_bar;

        Eigen::VectorXd s6(6), e6(6);
        if (dim == 3) {
            for (int32_t i = 0; i < 6; ++i) { s6[i] = s_bar[i]; e6[i] = e_bar[i]; }
        } else {
            expand_2d_stress_strain(s_bar, e_bar, ctx.plane_mode, ctx.nu, s6, e6);
        }

        out.elem_strain.row(e) = e6;
        out.elem_stress.row(e) = s6;
        out.elem_von_mises[e] = von_mises(s6);
        out.elem_principal.row(e) = principal_stresses(s6);

        if (npe * dim != ndof) throw FemError("Internal: kernel metadata mismatch");
    }

    out.displacement.resize(in.n_nodes, 3);
    out.displacement.setZero();
    for (int32_t n = 0; n < in.n_nodes; ++n) {
        for (int32_t d = 0; d < in.dimension; ++d) {
            const int32_t dof = in.dof_map(n, d);
            if (dof >= 0) out.displacement(n, d) = U[dof];
        }
    }
}

Eigen::VectorXd compute_reactions(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd& U,
    const Eigen::VectorXd& F) {
    return K * U - F;
}

SolverOutput fem_solve(const SolverInput& in) {
    SolverOutput out;
    validate_input(in);
    if (in.n_dofs_free < 0) {
        throw FemError("Negative free DOFs");
    }

    const ElementRegistry& registry = ElementRegistry::instance();

    const Eigen::SparseMatrix<double> K = assemble_stiffness(in, registry);
    const Eigen::VectorXd F = assemble_forces(in, registry);

    const ReducedSystem reduced = apply_boundary_conditions(K, F, in);

    try {
        const Eigen::VectorXd Ur = solve_linear_system(reduced);
        const Eigen::VectorXd U = reconstruct_full_displacement(Ur, reduced, in.n_dofs_total);

        compute_element_results(in, U, registry, out);

        const Eigen::VectorXd R = K * U - F;
        out.reaction.resize(in.n_nodes, 3);
        out.reaction.setZero();
        for (size_t i = 0; i < reduced.constrained_dofs.size(); ++i) {
            const int32_t dof = reduced.constrained_dofs[i];
            int32_t n = 0, d = 0;
            bool found = false;
            for (n = 0; n < in.n_nodes && !found; ++n) {
                for (d = 0; d < 3; ++d) {
                    if (in.dof_map(n, d) == dof) { found = true; break; }
                }
            }
            if (found) out.reaction(--n, d) = R[dof];
        }

        const int32_t nf = static_cast<int32_t>(reduced.free_dofs.size());
        double resid_sq = 0.0;
        for (int32_t i = 0; i < nf; ++i) {
            const double r = R[reduced.free_dofs[i]];
            resid_sq += r * r;
        }
        const double resid = std::sqrt(resid_sq);
        const double normF_free = [&]() {
            double s = 0.0;
            for (int32_t i = 0; i < nf; ++i) {
                const double f = F[reduced.free_dofs[i]];
                s += f * f;
            }
            return std::sqrt(s);
        }();
        out.residual_norm = (normF_free > 1e-12) ? resid / normF_free : resid;
        out.converged = out.residual_norm < 1e-8;
    } catch (const std::exception&) {
        out.converged = false;
        out.residual_norm = 0.0;
    }
    return out;
}

} // namespace fem