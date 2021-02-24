/**
 * \file       euclidean.hpp
 * \author     Bryan Flynt
 * \date       Feb 17, 2021
 */
#pragma once

#include <cassert>
#include <cmath>

namespace tbox {
namespace interpolation {
namespace norm {

template <typename T>
struct Euclidean {
    template <typename ArrayType>
    T operator()(const ArrayType& Pa, const ArrayType& Pb) const {
        assert(Pa.size() == Pb.size());
        using std::sqrt;
        T ans(0);
        for (std::size_t i = 0; i < Pa.size(); ++i) {
            ans += T((Pa[i] - Pb[i]) * (Pa[i] - Pb[i]));
        }
        return sqrt(ans);
    }
};

}  // namespace norm
}  // namespace interpolation
}  // namespace tbox
