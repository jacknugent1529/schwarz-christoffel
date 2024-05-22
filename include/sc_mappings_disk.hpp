#pragma once

#include <vector>
#include "polar.hpp"
#include "utils.hpp"
#include <gsl/gsl_integration.h>
#include <complex>
using std::complex;

template<typename T>
complex<T> schwarz_christoffel(std::vector<T> A_ks, std::vector<T> beta_ks, complex<T> z, gsl_workspace w);

template<typename T>
complex<T> schwarz_christoffel_singularity(std::vector<T> A_ks, std::vector<T> beta_ks, complex<T> z, int end_index, gsl_workspace w);
