/**
 * \file       rbf.hpp
 * \author     Bryan Flynt
 * \date       Feb 17, 2021
 */
#pragma once

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cassert>
#include <cmath>
namespace ublas = boost::numeric::ublas;


namespace geoflow {
namespace tbox {
namespace rbf {
namespace detail {


/** \brief decompose the symmetric positive definit matrix A into product L L^T.
     *
     * \param MATRIX type of input matrix 
     * \param TRIA type of lower triangular output matrix
     * \param A square symmetric positive definite input matrix (only the lower triangle is accessed)
     * \param L lower triangular output matrix 
     * \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
     */
template <class MATRIX, class TRIA>
size_t cholesky_decompose(const MATRIX& A, TRIA& L) {
    using namespace ublas;

    typedef typename MATRIX::value_type T;

    assert(A.size1() == A.size2());
    assert(A.size1() == L.size1());
    assert(A.size2() == L.size2());

    const size_t n = A.size1();

    for (size_t k = 0; k < n; k++) {
        double qL_kk = A(k, k) - inner_prod(project(row(L, k), range(0, k)),
                                            project(row(L, k), range(0, k)));

        if (qL_kk <= 0) {
            return 1 + k;
        } else {
            double L_kk = sqrt(qL_kk);
            L(k, k) = L_kk;

            matrix_column<TRIA> cLk(L, k);
            project(cLk, range(k + 1, n)) = (project(column(A, k), range(k + 1, n)) - prod(project(L, range(k + 1, n), range(0, k)),
                                                                                           project(row(L, k), range(0, k)))) /
                                            L_kk;
        }
    }
    return 0;
}

/** \brief decompose the symmetric positive definit matrix A into product L L^T.
     *
     * \param MATRIX type of matrix A
     * \param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
     * \param A output: the lower triangle of A is replaced by the cholesky factor
     * \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
     */
template <class MATRIX>
size_t cholesky_decompose(MATRIX& A) {
    using namespace ublas;

    typedef typename MATRIX::value_type T;

    const MATRIX& A_c(A);

    const size_t n = A.size1();

    for (size_t k = 0; k < n; k++) {
        double qL_kk = A_c(k, k) - inner_prod(project(row(A_c, k), range(0, k)),
                                              project(row(A_c, k), range(0, k)));

        if (qL_kk <= 0) {
            return 1 + k;
        } else {
            double L_kk = sqrt(qL_kk);

            matrix_column<MATRIX> cLk(A, k);
            project(cLk, range(k + 1, n)) = (project(column(A_c, k), range(k + 1, n)) - prod(project(A_c, range(k + 1, n), range(0, k)),
                                                                                             project(row(A_c, k), range(0, k)))) /
                                            L_kk;
            A(k, k) = L_kk;
        }
    }
    return 0;
}

/** \brief decompose the symmetric positive definit matrix A into product L L^T.
     *
     * \param MATRIX type of matrix A
     * \param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
     * \param A output: the lower triangle of A is replaced by the cholesky factor
     * \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
     */
template <class MATRIX>
size_t incomplete_cholesky_decompose(MATRIX& A) {
    using namespace ublas;

    typedef typename MATRIX::value_type T;

    // read access to a const matrix is faster
    const MATRIX& A_c(A);

    const size_t n = A.size1();

    for (size_t k = 0; k < n; k++) {
        double qL_kk = A_c(k, k) - inner_prod(project(row(A_c, k), range(0, k)),
                                              project(row(A_c, k), range(0, k)));

        if (qL_kk <= 0) {
            return 1 + k;
        } else {
            double L_kk = sqrt(qL_kk);

            // aktualisieren
            for (size_t i = k + 1; i < A.size1(); ++i) {
                T* Aik = A.find_element(i, k);

                if (Aik != 0) {
                    *Aik = (*Aik - inner_prod(project(row(A_c, k), range(0, k)),
                                              project(row(A_c, i), range(0, k)))) /
                           L_kk;
                }
            }

            A(k, k) = L_kk;
        }
    }

    return 0;
}

/** \brief solve system L L^T x = b inplace
     *
     * \param L a triangular matrix
     * \param x input: right hand side b; output: solution x
     */
template <class TRIA, class VEC>
void cholesky_solve(const TRIA& L, VEC& x, ublas::lower) {
    using namespace ublas;
    //   ::inplace_solve(L, x, lower_tag(), typename TRIA::orientation_category () );
    inplace_solve(L, x, lower_tag());
    inplace_solve(trans(L), x, upper_tag());
}

}  // namespace detail

namespace kernel {


template<typename T>
T ThinPlateSpline(const T rad) {
    using std::log;
    const T ans = rad*rad*log(rad);
    return std::isnan(ans) ? 0: ans;
}

template<typename T>
T InverseMultiquadratic(const T rad) {
    using std::sqrt;
    return static_cast<T>(1) / sqrt(1 + rad*rad);
}

template<typename T>
T InverseQuadratic(const T rad) {
    return static_cast<T>(1) / (static_cast<T>(1) + rad*rad);
}

template<typename T>
T Gaussian(const T rad) {
    using std::exp;
    return exp(-rad*rad);
}

} // namespace kernel




template <typename T>
class Interpolator {

public:

    Interpolator(const std::function<T(const T)>& rbf_kernel) 
        : kernel_(rbf_kernel) {
    }

    template <typename Coordinates, typename Solutions>
    void setData(const Coordinates& points, const Solutions& values) {
        const auto num_sources = points.size();
        const auto num_dimensions = points[0].size();
        const auto num_variables = 1;

        source_values_.resize(num_sources);
        source_points_.resize(num_sources, num_dimensions);

        for (auto p = 0; p < num_sources; ++p) {
            source_values_(p) = values[p];
            for (auto d = 0; d < num_dimensions; ++d) {
                source_points_(p, d) = points[p][d];
            }
        }
    }

    void computeWeights() {
        auto const num_points = source_points_.size1();

        css_.resize(num_points, num_points);

        // Place radius distances within css to start
        epsilon_ = 0;
        for (auto i = 0; i < num_points; ++i) {
            for (auto j = i; j < num_points; ++j) {
                css_(i, j) = ublas::norm_2(ublas::row(source_points_, i) - ublas::row(source_points_, j));
                epsilon_ += css_(i, j);
            }
        }
        epsilon_ = (static_cast<T>(num_points * (num_points + 1)) / 2) / epsilon_;

        // Calc RBF using scaled radius
        for (auto i = 0; i < num_points; ++i) {
            for (auto j = i; j < num_points; ++j) {
                css_(i, j) = kernel_(css_(i, j) * epsilon_);
                css_(j, i) = css_(i, j);
            }
        }

        std::cout << "Css = " << std::endl;
        for(auto r = 0; r < css_.size1(); ++r){
            for(auto c = 0; c < css_.size2(); ++c){
                std::cout << css_(r,c) << "  ";
            }
            std::cout << std::endl;
        }

        weights_ = source_values_;
        // ublas::permutation_matrix<std::size_t> pm(css_.size1());
        // ublas::lu_factorize(css_, pm);
        // ublas::lu_substitute(css_, pm, weights_);

        // Factor Matrix Css into Lower triangular
        ublas::matrix<T> L(num_points, num_points);
        detail::cholesky_decompose(css_, L);
        detail::cholesky_solve(L, weights_, ublas::lower());     


        std::cout << "Weights = " << std::endl;
        for(auto i = 0; i < weights_.size(); ++i){
            std::cout << weights_(i) << std::endl;
        }
    }

    template <typename Coordinate>
    T calcValue(const Coordinate& point) const {

        const auto num_sources    = source_points_.size1();
        const auto num_dimensions = point.size();
        const auto num_variables  = 1;

        // Copy interpolation point over
        ublas::vector<T> ipoint(num_dimensions);
        for (auto d = 0; d < num_dimensions; ++d) {
            ipoint(d) = point[d];
        }

        // Calculate the Afs matrix (vector)
        ublas::vector<T> Afs(num_sources);
        for (auto i = 0; i < num_sources; ++i) {
            auto rad = ublas::norm_2(ipoint - ublas::row(source_points_, i));
            Afs(i)   = kernel_(rad * epsilon_);
        }

        // Dot Product to get answer
        return ublas::inner_prod(Afs, weights_);
    }

   private:
    std::function<T(const T)> kernel_;
    ublas::matrix<T> source_points_;
    ublas::vector<T> source_values_;
    ublas::matrix<T> css_;
    ublas::vector<T> weights_;
    T epsilon_;
};

}  // namespace rbf
}  // namespace tbox
}  // namespace geoflow
