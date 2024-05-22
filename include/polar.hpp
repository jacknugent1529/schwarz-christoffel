#pragma once
#include <cmath>
#include <complex>
using std::complex;

template<typename T>
class Polar {
public: T r;
    T theta;

    Polar(const complex<T> z) {
        r = abs(z);
        theta = arg(z);
    }

    Polar(const T r_, const T theta_) {
        r = r_;
        theta = theta_;
    }

    complex<T> to_complex() const {
        return std::polar(r, theta);
    }

    Polar<T> operator*(const Polar<T>other) const {
        T r_new = r * other.r;
        T theta_new = theta + other.theta;
        return Polar<T>(r_new, theta_new);
    }

    Polar<T> inv() const {
        return Polar<T>(1 / r, -theta);
    }
};

