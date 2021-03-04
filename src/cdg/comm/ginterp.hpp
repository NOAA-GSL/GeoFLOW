//==================================================================================
// Module       : ginterp.hpp
// Date         : 1/6/21 (BTF)
// Description  : Encapsulates the methods and data associated with
//                a geometry--free global interpolation operator
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <array>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "boost/mpi.hpp"
#include "gcomm.hpp"
#include "gtmatrix.hpp"
#include "gtvector.hpp"
#include "tbox/assert.hpp"
#include "tbox/pio.hpp"
#include "tbox/spatial.hpp"
#include "tbox/tracer.hpp"

namespace detail_extractor {
template <typename PairType>
struct pair_extractor {
    using bound_type = typename PairType::first_type;
    using key_type = typename PairType::second_type;

    const bound_type& operator()(PairType const& pair) const {
        return pair.first;
    }
};
}  // namespace detail_extractor

template <typename ValueType>
class GInterp {
   public:
    GInterp() = default;
    GInterp(const GInterp& other) = default;
    GInterp(GInterp&& other) = default;
    ~GInterp() = default;
    GInterp& operator=(const GInterp& other) = default;
    GInterp& operator=(GInterp&& other) = default;

    // template <typename Coordinates>
    // GBOOL init(const ValueType tolerance, Coordinates& source_coords, Coordinates& target_coords);

    // template <typename Solutions>
    // Solutions interpolate(const Solutions& source_soln);

    template <typename Coordinates, typename Solutions>
    static void interpolate(const ValueType tolerance,
                            const Coordinates& source_coords, const Solutions& source_soln,
                            const Coordinates& target_coords, Solutions& target_soln);

    void display() const;

   private:
    using rank_type = int;
    using size_type = std::size_t;
    using value_type = ValueType;

    template <typename Coordinates, typename Solutions>
    static void serial_interp_(const Coordinates& search_location, const Solutions& search_value, const Coordinates& target_coords, Solutions& target_value);

};

    template <typename T>
    template <typename Coordinates, typename Solutions>
    void GInterp<T>::interpolate(const T tolerance,
                                 const Coordinates& source_coords, const Solutions& source_soln,
                                 const Coordinates& target_coords, Solutions& target_soln) {
        GEOFLOW_TRACE();
        using namespace geoflow::tbox;

        // Types
        using index_bound_type = tbox::spatial::bound::Box<value_type, GDIM>;
        using coord_value_type = typename Coordinates::value_type;
        using soln_value_type = typename Solutions::value_type;
        using index_value_type = std::pair<index_bound_type, soln_value_type>;  // (bound, values)
        using index_extractor_type = detail_extractor::pair_extractor<index_value_type>;
        using shared_index_type = tbox::spatial::shared::index::RTree<index_value_type, index_extractor_type>;

        // Get MPI communicator, etc.
        namespace mpi = boost::mpi;
        mpi::communicator world;
        auto my_rank = world.rank();
        auto num_ranks = world.size();

        // ----------------------------------------------------------
        //      Build Indexer for each Local Source Location
        // ----------------------------------------------------------

        // Build "index_value_type" for each local coordinate and place into a spatial index
        GEOFLOW_TRACE_START("Build Local Source Index");
        shared_index_type local_source_indexer;
        for (auto i = 0; i < source_coords.size(); ++i) {
            // Build a slightly larger (tolerance) bounding box
            auto min_xyz = source_coords[i];
            auto max_xyz = source_coords[i];
            for (auto d = 0; d < min_xyz.size(); ++d) {
                min_xyz[d] -= tolerance;
                max_xyz[d] += tolerance;
            }

            // Place into coordinate index tree
            auto source_bound = index_bound_type(min_xyz, max_xyz);
            auto index_value = index_value_type(source_bound, source_soln[i]);
            local_source_indexer.insert(index_value);
        }
        GEOFLOW_TRACE_STOP();

        // ----------------------------------------------------------
        //      Build Bound for each Local Target Location
        // ----------------------------------------------------------

        // Build "index_value_type" for each local coordinate and place into a spatial index
        GEOFLOW_TRACE_START("Build Local Target Bounds");
        index_bound_type target_region;
        target_region.reset();
        std::vector<index_bound_type> target_bounds(target_coords.size());
        for (auto i = 0; i < target_coords.size(); ++i) {
            target_bounds[i] = index_bound_type(target_coords[i], target_coords[i]);
            target_region.stretch(target_bounds[i]);
        }
        GEOFLOW_TRACE_STOP();

        // ----------------------------------------------------------
        //  Gather/Scatter my Local Bounds to All Ranks
        // ----------------------------------------------------------

        // Perform MPI_Allgather so everyone gets a bounds region for each processor
        GEOFLOW_TRACE_START("AllGather Target Bounds");
        std::vector<index_bound_type> target_bounds_by_rank;
        mpi::all_gather(world, target_region, target_bounds_by_rank);
        GEOFLOW_TRACE_STOP();

        // pio::perr << "Source Bounds = " << local_source_indexer.bounds() << std::endl;
        // for(auto i = 0; i < target_bounds_by_rank.size(); ++i){
        //     pio::perr << "Target Bounds for Rank " << i << " = " << target_bounds_by_rank[i] << std::endl;
        // }

        // ----------------------------------------------------------
        //  Build and send buffer of intersecting source values
        // ----------------------------------------------------------

        // For each Ranks Bounding Region
        // - Build list of local source values to send to each overlapping target rank
        // - Submit a receive request to get indexed values from the rank
        // - Submit a send request to send indexed values from this rank
        GEOFLOW_TRACE_START("Gather/Scatter Source Values");
        std::map<rank_type, mpi::request> send_requests;
        std::map<rank_type, mpi::request> recv_requests;
        std::map<rank_type, std::vector<index_value_type>> send_to_ranks;
        std::map<rank_type, std::vector<index_value_type>> recv_from_ranks;
        for (rank_type rank = 0; rank < num_ranks; ++rank) {
            // Construct box slightly larger than target region
            auto target_bound = target_bounds_by_rank[rank];
            target_bound.scale(2.0);

            // if(my_rank == 0){
            //     pio::perr << "Stretched Target Bound = " << rank << " " << target_bound << std::endl;
            // }

            // Get all sources which intersect with target region
            std::vector<index_value_type> search_results;
            local_source_indexer.query(tbox::spatial::shared::predicate::Intersects(target_bound), std::back_inserter(search_results));

            // Form send/recv buffers
            // - Just directly copy if already on rank
            if (search_results.size() > 0) {
                if (my_rank == rank) {
                    recv_from_ranks.emplace(my_rank, search_results);
                } else {
                    int tag = 0;
                    send_to_ranks.emplace(rank, search_results);
                    recv_requests[rank] = world.irecv(rank, 0, recv_from_ranks[rank]);
                    send_requests[rank] = world.isend(rank, 0, send_to_ranks[rank]);
                }
            }
        }
        GEOFLOW_TRACE_STOP();

        // ----------------------------------------------------------
        //  Build "Semi-Global" Index of bounds+values
        // ----------------------------------------------------------

        GEOFLOW_TRACE_START("Build Semi-Global Source Index");
        shared_index_type global_source_indexer;
        global_source_indexer.insert(recv_from_ranks[my_rank].begin(), recv_from_ranks[my_rank].end());
        for (auto& [rank, req] : recv_requests) {
            req.wait();  // Wait to receive data
            auto& pairs_from_rank = recv_from_ranks[rank];
            global_source_indexer.insert(pairs_from_rank.begin(), pairs_from_rank.end());
        }
        for (auto& [rank, req] : send_requests) {
            req.wait();  // Clear out the send requests
        }
        GEOFLOW_TRACE_STOP();

        // pio::perr << "Global Bounds = " << global_source_indexer.bounds() << std::endl;

        // ----------------------------------------------------------
        //  Perform Interpolation using K-Nearest Neighbors
        // ----------------------------------------------------------

        GEOFLOW_TRACE_START("Nearest Neighbor Interpolation");
        const int k_nearest = 30;  // Number of nearest neighbors
        for (auto i = 0; i < target_bounds.size(); ++i) {
            // Get Nearest for each Target Location
            std::vector<index_value_type> search_results;
            local_source_indexer.query(tbox::spatial::shared::predicate::Nearest(target_bounds[i], k_nearest), std::back_inserter(search_results));

            // Eliminate Duplicate Locations and get location + value
            auto last = std::unique(search_results.begin(), search_results.end(), [](auto& a, auto& b) {
                return Intersects(a.first, b.first);
            });
            search_results.erase(last, search_results.end());

            // Extract Coordinates and Values from pair
            std::vector<coord_value_type> search_locations;
            std::vector<soln_value_type> search_values;
            for (auto& [bound, value] : search_results) {
                search_values.emplace_back(value);

                auto location = coord_value_type();
                for (auto d = 0; d < location.size(); ++d) {
                    location[d] = bound.center(d);
                }
                search_locations.emplace_back(location);
            }

            std::vector<coord_value_type> target_coords_vec = {target_coords[i]};
            std::vector<soln_value_type>  target_soln_vec   = {target_soln[i]};
            GInterp<T>::serial_interp_(search_locations, search_values, target_coords_vec, target_soln_vec);
            target_soln[i] = target_soln_vec[0];
        }
        GEOFLOW_TRACE_STOP();
    }

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cassert>
    namespace ublas = boost::numeric::ublas;

    namespace detail {

    template <typename Point>
    typename Point::value_type
    Radius(const Point& a, const Point& b) {
        using value_type = typename Point::value_type;
        value_type rad = 0;
        for (auto i = 0; i < a.size(); ++i) {
            rad += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return std::sqrt(rad);
    }

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
    void
    cholesky_solve(const TRIA& L, VEC& x, ublas::lower) {
        using namespace ublas;
        //   ::inplace_solve(L, x, lower_tag(), typename TRIA::orientation_category () );
        inplace_solve(L, x, lower_tag());
        inplace_solve(trans(L), x, upper_tag());
    }

    }  // namespace detail




    template <typename T>
    template <typename Coordinates, typename Solutions>
    void
    GInterp<T>::serial_interp_(const Coordinates& search_location, const Solutions& search_value, const Coordinates& target_coords, Solutions& target_value) {
        const auto num_sources = search_location.size();
        const auto num_dimensions = search_location[0].size();
        const auto num_variables = 1;

        // Get epsilon scaling parameter
        value_type avg_distance = 0;
        for (auto i = 0; i < num_sources; ++i) {
            avg_distance += detail::Radius(target_coords[0], search_location[i]);
        }
        avg_distance /= num_sources;
        value_type epsilon = 1.0 / avg_distance;

        ublas::matrix<value_type> weights(num_sources, num_variables);
        ublas::matrix<value_type> A(num_sources, num_sources);
        ublas::matrix<value_type> L(num_sources, num_sources);
        ublas::vector<value_type> Afs(num_sources);
        ublas::vector<value_type> soln(num_variables);

        // Copy solution into Matrix
        for (auto i = 0; i < num_sources; ++i) {
            for (auto j = 0; j < num_variables; ++j) {
                weights(i, j) = search_value[i];
            }
        }

        // Build RBF
        for (auto i = 0; i < num_sources; ++i) {
            for (auto j = 0; j < num_sources; ++j) {
                auto rad = detail::Radius(search_location[i], search_location[j]);
                A(i, j) = std::exp(-epsilon * rad * epsilon * rad);
            }
        }

        // Factor Matrix A
        size_t res = detail::cholesky_decompose(A, L);

        // Solve weight for each variable
        for (auto i = 0; i < num_variables; ++i) {
            ublas::matrix_column<ublas::matrix<value_type>> mc(weights, i);
            detail::cholesky_solve(L, mc, ublas::lower());
        }

        // Build Afs Matrix
        for (auto i = 0; i < num_sources; ++i) {
            auto rad = detail::Radius(target_coords[0], search_location[i]);
            Afs(i) = std::exp(-epsilon * rad * epsilon * rad);
        }

        // Copy solution out
        soln = ublas::prod(ublas::trans(Afs), weights);
        for (auto i = 0; i < num_variables; ++i) {
            target_value[i] = soln(i);
        }
    }
