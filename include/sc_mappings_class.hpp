#pragma once
#include "sc_mappings_disk.hpp"
#include <iostream>
using std::complex;

template<typename T>
class SchwarzChristoffel {
private:
    std::vector<T> A_ks;
    std::vector<T> beta_ks;
    std::vector<complex<T>> A_k_values;
    std::vector<gsl_workspace> workspaces;
    complex<T> c1;
    complex<T> c2;

public:
    SchwarzChristoffel(std::vector<T> A_ks_, std::vector<T> beta_ks_, complex<T> c1_ = 0., complex<T> c2_ = 0.) : A_ks(A_ks_), beta_ks(beta_ks_), c1(c1_), c2(c2_) {
        // Sort the A_ks in ascending order
        std::vector<std::pair<T, T>> pairs;
        for (int i = 0; i < A_ks.size(); i++) {
            pairs.emplace_back(A_ks[i], beta_ks[i]);
        }
        std::sort(pairs.begin(), pairs.end());
        
        A_ks.clear();
        beta_ks.clear();
        for (const auto& pair : pairs) {
            A_ks.push_back(pair.first);
            beta_ks.push_back(pair.second);
        }

        // Calculate the SC mapping values at each A_k
        for (int i = 0; i < A_ks.size(); i++) {
            complex<T> zero (0.);
            auto res = schwarz_christoffel_singularity(A_ks, beta_ks, zero, i, {});
            A_k_values.push_back(res);
        }
        
    }

    complex<T> operator()(complex<T> z) {
        // Find the index of the nearest A_k to z
        int nearest_index = 0;
        T min_distance = abs(z - std::polar(1.,M_PI * A_ks[0]));
        for (int i = 1; i < A_ks.size(); i++) {
            T distance = abs(z - std::polar(1.,M_PI * A_ks[i]));
            if (distance < min_distance) {
                min_distance = distance;
                nearest_index = i;
            }
        }

        if (min_distance < 1e-4) {
            return c1 * A_k_values[nearest_index] + c2;
        }


        // Integrate from the nearest A_k to z
        complex<T> int_z_end = -schwarz_christoffel_singularity(A_ks, beta_ks, z, nearest_index, {});

        return c1 * (int_z_end + A_k_values[nearest_index]) + c2;
    }
};
