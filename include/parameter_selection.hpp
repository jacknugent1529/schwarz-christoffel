#pragma once
#include <vector>
#include <complex>
using std::complex;
#include <tuple>

std::tuple<std::vector<double>, std::vector<double>, std::complex<double>, std::complex<double>> solve_constraints(std::vector<complex<double>> points);
