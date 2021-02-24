/**
 * \file       gaussian.hpp
 * \author     Bryan Flynt
 * \date       Feb 17, 2021
 */
#pragma once

#include <cassert>
#include <cmath>

namespace tbox {
namespace interpolation {
namespace kernel {

template<typename T>
struct Gaussian {

    Gaussian() = delete;
    Gaussian(const Gaussian& other) = default;
    Gaussian(Gaussian&& other) = default;
    ~Gaussian() = default;
    Gaussian& operator=(const Gaussian& other) = default;
    Gaussian& operator=(Gaussian&& other) = default;

    Gaussian(const T epsilon = 1.0) : epsilon_(epsilon) {
    }

    T operator()(const T rad) const {
        using std::exp;
        return exp(-(epsilon_*rad * epsilon_*rad));
    }

private:
    T epsilon_;
};

}  // namespace kernel
}  // namespace interpolation
}  // namespace tbox