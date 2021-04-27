/*
 * math.hpp
 *
 *  Created on: April 27, 2021
 *      Author: bflynt
 */
namespace geoflow {
namespace tbox {
namespace math {

template <typename T>
constexpr T ipow(T x, long int exp) {
    if (exp == 0) {
        return x < 0 ? -1.0 : 1.0;
    }
    int sign = 1;
    if (exp < 0) {
        sign = -1;
        exp = -exp;
    }
    double ret = x;
    while (--exp) {
        ret *= x;
    }
    return sign > 0 ? ret : 1.0 / ret;
}

}  // namespace math
}  // namespace tbox
}  // namespace geoflow
