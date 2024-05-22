#include "sc_mappings_disk.hpp"
#include <iostream>

template<typename T>
inline Polar<T> pow(Polar<T> z, T beta) {
    // the discontinuity will be at theta = 3pi/2
    if (std::abs(z.theta + 3 * M_PI_4) < M_PI_4) {
        Polar<T> offset(static_cast<T>(1), static_cast<T>(2 * M_PI));
        z = z * offset;
    }

    return Polar<T>(pow(z.r, beta), z.theta * beta);
}

template<typename T>
complex<T> sc_mapping_disk(
    std::vector<T> A_ks, 
    std::vector<T> beta_ks, 
    complex<T> z, 
    T t
) {
    complex<T> z0 = 0.;
    complex<T> zeta = z0 + t * (z - z0);
    Polar<T> out(1,0);
    for (int i = 0; i < (int) A_ks.size(); i++) {
        complex<T> B_k = std::polar(1., M_PI * A_ks[i]);
        out = out * pow(Polar<T>(1. - zeta / B_k), beta_ks[i]).inv();
    }
    out = out * z;

    return out.to_complex(); 
}

template<typename T>
complex<T> sc_mapping_disk_singularity(
    std::vector<T> A_ks, 
    std::vector<T> beta_ks, 
    complex<T> z0, 
    int end_index,
    T t
) {
    auto z = std::polar(1., M_PI * A_ks[end_index]);
    complex<T> zeta = z0 + t * (z - z0);
    Polar<T> out(1,0);
    for (int i = 0; i < (int) A_ks.size(); i++) {
        if (i == end_index) continue;
        complex<T> B_k = std::polar(1., M_PI * A_ks[i]);
        out = out * pow(Polar<T>(1. - zeta / B_k), beta_ks[i]).inv();
    }
    out = out * (z - z0);
    complex<T> B_k = std::polar(1., M_PI * A_ks[end_index]);
    out = out * pow(Polar<T>(1. - z0 / B_k), beta_ks[end_index]).inv();

    return out.to_complex(); // should be multiplied by (1 - t)^-beta_i with i = end_index
}

template<typename T>
struct sc_params {std::vector<T> A_ks; std::vector<T> beta_ks; complex<T> z; int end_index;};

template<typename T>
complex<T> sc_mapping_disk_f(T t, void *params_) {
    sc_params<T> *params = (sc_params<T>*) params_;
    return sc_mapping_disk(params->A_ks, params->beta_ks, params->z, t);
}
template<typename T> T sc_mapping_disk_real(T t, void *params_) { return sc_mapping_disk_f(t, params_).real(); }
template<typename T> T sc_mapping_disk_imag(T t, void *params_) { return sc_mapping_disk_f(t, params_).imag(); }

template<typename T>
complex<T> sc_mapping_disk_singularity_f(T t, void *params_) {
    sc_params<T> *params = (sc_params<T>*) params_;
    return sc_mapping_disk_singularity(params->A_ks, params->beta_ks, params->z, params->end_index, t);
}
template<typename T> T sc_mapping_disk_singularity_real(T t, void *params_) { return sc_mapping_disk_singularity_f(t, params_).real(); }
template<typename T> T sc_mapping_disk_singularity_imag(T t, void *params_) { return sc_mapping_disk_singularity_f(t, params_).imag(); }

template<typename T>
complex<T> schwarz_christoffel(std::vector<T> A_ks, std::vector<T> beta_ks, complex<T> z, gsl_workspace w) {
    double result_real, _error_real;
    gsl_function F_real;
    F_real.function = &sc_mapping_disk_real;
    sc_params<T> params {A_ks, beta_ks, z, 0};
    F_real.params = &params;

    gsl_integration_qags(&F_real, 0, 1, 0, 1e-5, 1000, w.ws, &result_real, &_error_real);

    double result_imag, _error_imag;
    gsl_function F_imag;
    F_imag.function = &sc_mapping_disk_imag;
    F_imag.params = &params;

    gsl_integration_qags(&F_imag, 0, 1, 0, 1e-5, 1000, w.ws, &result_imag, &_error_imag);

    return std::complex<T>(result_real, result_imag);
}

template<typename T>
complex<T> schwarz_christoffel_singularity(std::vector<T> A_ks, std::vector<T> beta_ks, complex<T> start, int end_index, gsl_workspace w) {
    w = setup_ws(-beta_ks[end_index]);

    double result_real, _error_real;
    gsl_function F_real;
    F_real.function = &sc_mapping_disk_singularity_real;
    sc_params<T> params {A_ks, beta_ks, start, end_index};
    F_real.params = &params;

    /* gsl_integration_qags(&F_real, 0, 1, 0, 1e-5, 1000, w.ws, &result_real, &_error_real); */
    gsl_integration_qaws(&F_real, 0, 1, w.table, 1e-5, 1e-5, 1000, w.ws, &result_real, &_error_real);

    double result_imag, _error_imag;
    gsl_function F_imag;
    F_imag.function = &sc_mapping_disk_singularity_imag;
    F_imag.params = &params;

    gsl_integration_qaws(&F_imag, 0, 1, w.table, 1e-5, 1e-5, 1000, w.ws, &result_imag, &_error_imag);

    return std::complex<T>(result_real, result_imag);
}

template complex<double> schwarz_christoffel(std::vector<double> A_ks, std::vector<double> beta_ks, complex<double> z, gsl_workspace w);
template complex<double> schwarz_christoffel_singularity(std::vector<double> A_ks, std::vector<double> beta_ks, complex<double> z, int end_index, gsl_workspace w);
