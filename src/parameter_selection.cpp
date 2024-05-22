#include "parameter_selection.hpp"
#include "utils.hpp"
#include "sc_mappings_disk.hpp"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>
#include <iostream>

std::vector<double> gsl_to_std_vector(const gsl_vector* gsl_vec) {
    std::vector<double> std_vec(gsl_vec->size);
    for (size_t i = 0; i < gsl_vec->size; ++i) {
        std_vec[i] = gsl_vector_get(gsl_vec, i);
    }
    return std_vec;
}

double logistic(double x) {
    return 1 / (1 + std::exp(-x));
}

struct sc_params {
    std::vector<complex<double>> points;
    std::vector<double> betas;
};

int sc_constraints(const gsl_vector* x, void* params_, gsl_vector* f) {
    sc_params* params = (sc_params*) params_;
    auto points = params->points;
    auto betas = params->betas;


    auto tmp = gsl_to_std_vector(x);
    std::sort(tmp.begin(), tmp.end());
    for (auto x : tmp) {
        std::cout << "xi " << x << std::endl;
    }
    std::vector<double> A_ks;
    A_ks.push_back(0.);
    A_ks.push_back(1.5);
    A_ks.push_back(1.);
    for (auto x : tmp) {
        A_ks.push_back(logistic(x));
    }

    for (auto a : A_ks) { std::cout << a << std::endl; }

    gsl_workspace ws = setup_ws();
    
    // TODO: check this
    std::cout << "--------------------\n";
    for (auto j = 2; j < points.size() - 2; j++) {
        /* complex<double> zj = std::polar(1., A_ks[j] * M_PI); */
        /* complex<double> zj1 = std::polar(1., A_ks[j+1] * M_PI); */
        complex<double> zero (0.);
        auto int_0_zj = schwarz_christoffel_singularity(A_ks, betas, zero, j, ws);
        auto int_0_zj1 = schwarz_christoffel_singularity(A_ks, betas, zero, j+1, ws);

        auto int_0_z0 = schwarz_christoffel_singularity(A_ks, betas, zero, 0, ws);
        auto int_0_z1 = schwarz_christoffel_singularity(A_ks, betas, zero, 1, ws);

        std::cout << int_0_zj1 << " " <<  int_0_zj << " " << std::abs(int_0_zj1 - int_0_zj) << " " << std::abs(int_0_z1 - int_0_z0) << std::endl;
        double lhs = std::abs(int_0_zj1 - int_0_zj) / std::abs(int_0_z1 - int_0_z0);
        double rhs = std::abs(points[j+1] - points[j]) / std::abs(points[1] - points[0]);
        std::cout << "rhs " << rhs << std::endl;
        std::cout << "lhs " << rhs << std::endl;

        gsl_vector_set(f, j-1, lhs - rhs);
    }

    return GSL_SUCCESS;
}

std::vector<double> angles(std::vector<complex<double>> points) {
    auto n = points.size();
    std::vector<double> betas;
    for (int i = 0; i < n; i++) {
        auto w1 = points[i];
        auto w2 = points[(i + 1) % n];
        auto w3 = points[(i + 2) % n];
        /* std::cout << "-----------------" << std::endl; */
        /* std::cout << w1 << " " << w2 << " " << w3 << std::endl; */

        w1 = w1 - w2;
        w3 = w3 - w2;
        /* std::cout << w1 << " " << 0 << " " << w3 << std::endl; */
        auto alpha = std::polar(1., -std::arg(w1));
        w3 = w3 * alpha;
        /* std::cout << w3 << " " << w1 * alpha << std::endl; */
        auto theta = std::atan2(w3.imag(), w3.real()) / M_PI;
        if (theta < 0) theta += 2;
        betas.push_back(1 - theta);
    }

    return betas;
}

std::tuple<std::vector<double>, std::vector<double>, std::complex<double>, std::complex<double>> solve_constraints(std::vector<complex<double>> points) {
    auto betas = angles(points);
    for (auto b : betas) {
        std::cout << "beta " << b << std::endl;
    }
    const size_t n = points.size() - 3;  // Number of equations and variables
    gsl_multiroot_function f = {&sc_constraints, n, nullptr};

    sc_params params {points, betas};

    f.params = (void*) &params;

    gsl_vector* x = gsl_vector_alloc(n);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i, ((double) i + 1) / ((double) n + 2));  // Initial guess for x1
    }

    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrid;
    // TODO: see this https://stackoverflow.com/questions/24028760/gsl-multiroot-iteration-trying-nan
    /* const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_dnewton; */
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_fsolver_set(solver, &f, x);

    int status;
    size_t iter = 0;
    do {
        iter++;
        std::cout << "iter " << iter << std::endl;
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status)
            break;

        status = gsl_multiroot_test_residual(solver->f, 1e-5);
    } while (status == GSL_CONTINUE && iter < 1000);

    auto A_ks_ = gsl_to_std_vector(solver->x);
    std::vector<double> A_ks;
    gsl_multiroot_fsolver_free(solver);
    A_ks.push_back(0);
    A_ks.push_back(1.5);
    A_ks.push_back(1.);
    for (auto x : A_ks_) {
        A_ks.push_back(x);
    }
    gsl_vector_free(x);

    complex<double> zero (0.);
    auto p0 = points[0];
    auto p1 = points[1];
    auto q0 = schwarz_christoffel_singularity(A_ks, betas, zero, 0, {});
    auto q1 = schwarz_christoffel_singularity(A_ks, betas, zero, 1, {});

    std::vector<double> A {
        q0.real(), -q0.imag(), 1, 0,
        q0.imag(),  q0.real(), 0, 1,
        q1.real(), -q1.imag(), 1, 0,
        q1.imag(),  q1.real(), 0, 1,
    };
    std::vector<double> b_ {p0.real(), p0.imag(), p1.real(), p1.imag()} ;


    gsl_matrix_view m = gsl_matrix_view_array(A.data(), 4, 4);
    gsl_vector_view b = gsl_vector_view_array(b_.data(), 4);

    gsl_vector *v_ = gsl_vector_alloc(4);

    gsl_permutation *p = gsl_permutation_alloc(4);

    int s;
    gsl_linalg_LU_decomp(&m.matrix, p, &s);
    gsl_linalg_LU_solve(&m.matrix, p, &b.vector, v_);

    auto v = gsl_to_std_vector(v_);

    gsl_permutation_free(p);
    gsl_vector_fprintf(stdout, v_, "%g");
    gsl_vector_free(v_);

    std::cout << "p0, p1" << p0 << " " << p1 << std::endl;
    std::cout << "q0, q1" << q0 << " " << q1 << std::endl;
    std::cout << "c1 " << v[0] << " " << v[1] << std::endl;
    std::cout << "c2 " << v[2] << " " << v[3] << std::endl;


    return {A_ks, betas, {v[0], v[1]}, {v[2],v[3]} };
}

