//==================================================================================
// Module       : cdg_ggfx.cpp
// Date         : 1/19/21 (BTF)
// Description  : GeoFLOW test of new GGFX operations
// Copyright    : Copyright 2018. Colorado State University. All rights reserved
// Derived From : none.
//==================================================================================

#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>
#include "ggfx.hpp"
#include "tbox/pio.hpp"
#include "boost/mpi.hpp"



template<std::size_t N>
class RowMajorIndex final {
	static_assert(N > 0, "Cannot have 0 Rank Index");

public:

	// ====================================================
	// Types
	// ====================================================

	using size_type  = std::size_t;

	// ====================================================
	// Constructors
	// ====================================================

	RowMajorIndex()                             = delete;
	RowMajorIndex(const RowMajorIndex& other)   = default;
	RowMajorIndex(RowMajorIndex&& other)        = default;
	~RowMajorIndex()                            = default;

	RowMajorIndex(const std::array<size_type,N> shape) :
		shapes_(shape),
		linear_index_(0){
		this->calc_stride_();
		this->calc_index_();
	}

	// ====================================================
	// Operators
	// ====================================================

	RowMajorIndex& operator=(const RowMajorIndex& other) = default;
	RowMajorIndex& operator=(RowMajorIndex&& other)      = default;


	RowMajorIndex& operator=(const size_type& index){
		ASSERT(index < size());
		linear_index_ = index;
		this->calc_index_();
		return *this;
	}

	RowMajorIndex& operator++(){
		ASSERT(linear_index_ < size());
		linear_index_++;
		this->calc_index_();
		return *this;
	}

	RowMajorIndex operator++(int){
		ASSERT(linear_index_ < size());
		auto tmp = *this;
		linear_index_++;
		this->calc_index_();
		return tmp;
	}

	template<typename... Dims>
	size_type operator()(const Dims... args) {
		static_assert(sizeof...(args) == N);
		indexes_ = {static_cast<size_type>(args)...};
		linear_index_ = 0;
		for(size_type i = 0; i < N; ++i) {
			ASSERT(indexes_[i] < shapes_[i]);
			linear_index_ += indexes_[i] * strides_[i];
		}
		return linear_index_;
	}

	size_type operator[](const size_type rank) const {
		ASSERT(rank < N);
		return indexes_[rank];
	}

	// ====================================================
	// Conversion
	// ====================================================

	/** Implicit conversion to single index
	 */
	operator size_type() const {
		return linear_index_;
	}


	// ====================================================
	// Query
	// ====================================================

	size_type size() const {
		size_type sz = 1;
		for(size_type i = 0; i < N; ++i) {
			sz *= shapes_[i];
		}
		return sz;
	}

	size_type shape(const size_type i) const {
		ASSERT(i < N);
		return shapes_[i];
	}

	size_type stride(const size_type i) const {
		ASSERT(i < N);
		return strides_[i];
	}

	size_type index(const size_type i) const {
		ASSERT(i < N);
		return indexes_[i];
	}


	// ====================================================
	// PRIVATE
	// ====================================================

private:
	std::array<size_type,N> shapes_;
	std::array<size_type,N> strides_;
	std::array<size_type,N> indexes_;
	size_type               linear_index_;


	void calc_stride_() {
		strides_[N-1] = 1;
		for(size_type i = (N-1); i--> 0;) {
			strides_[i] = strides_[i+1] * shapes_[i+1];
		}
	}

	void calc_index_() {
		size_type fac = 1;
		for(size_type i = N; i--> 1;) {
			indexes_[i] = (linear_index_ / fac) % shapes_[i];
			fac *= shapes_[i];
		}
		indexes_[0] = linear_index_ / fac;
	}

};












int main(int argc, char **argv){
	namespace mpi  = boost::mpi;
	using namespace geoflow::tbox;
	constexpr int NDIM = GDIM;
	using size_type  = std::size_t;
	using value_type = double;
	using point_type = std::array<value_type, NDIM>;

	// Start Up MPI and find my place in it
    mpi::environment  env(argc,argv);
    mpi::communicator world;
	auto my_rank   = world.rank();
	auto num_ranks = world.size();
	pio::initialize(my_rank);
//	pio::logAllNodes("log");

	// Constants that control the region for each rank
	std::array<size_type,NDIM>  grid_dimensions;
	std::array<value_type,NDIM> grid_distance;
	grid_dimensions.fill(5);
	grid_distance.fill(1.0);

	// Build our Ranks grid region & spacing
	std::array<value_type,NDIM> sxyz; // grid start
	std::array<value_type,NDIM> dxyz; // grid spacing
	for(size_type i = 0; i < NDIM; ++i){
		dxyz[i] = grid_distance[i] / (grid_dimensions[i]-1);
		sxyz[i] = my_rank * dxyz[i] * (grid_dimensions[i]/2);
	}

	// Build our Cartesian Grid
	auto index = RowMajorIndex<NDIM>(grid_dimensions);
	auto num_points = index.size();
	std::vector<point_type> xyz(num_points);
	for(index = 0; index < num_points; ++index) {
		for(size_type d = 0; d < NDIM; ++d){
			xyz[index][d] = sxyz[d] + index[NDIM-1-d] * dxyz[d];
		}
	}

	// Init GGFX for each ranks region
	auto min_dist = std::accumulate(dxyz.begin(),dxyz.end(),9999.,[](auto a, auto b){
		return std::min(a,b);
	});
	value_type tolerance = 0.25 * min_dist;
	GGFX<value_type> ggfx;
	ggfx.init(tolerance,xyz);

	// Build Known Solution
	std::vector<value_type> soln(num_points, 0.0);
	for(size_type i = 0; i < num_points; ++i) {
		for(size_type d = 0; d < NDIM; ++d){
			soln[i] += d * xyz[i][d];
		}
	}
	const std::vector<value_type> original(soln);

	// Perform our Reduction
	ggfx.doOp(soln,GGFX<value_type>::Smooth());


	// Compare Solutions
	for(size_type i = 0; i < num_points; ++i) {
		if( original[i] != soln[i] ){
//			pio::pout << "ERROR: At index = " << i << std::endl;
			return 1;
		}
	}

//	ggfx.display();
//	std::vector<value_type> imult(num_points);
//	ggfx.get_imult(imult);
//	pio::plog << std::scientific;
//	for(index = 0; index < num_points; ++index) {
//		pio::plog << xyz[index][0] << " " << xyz[index][1] << " " << original[index] << "  " << soln[index] << "  " << imult[index] << std::endl;
//	}

	pio::finalize();
	return 0;
}
