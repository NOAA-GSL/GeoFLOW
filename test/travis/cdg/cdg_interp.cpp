//==================================================================================
// Module       : cdg_interp.cpp
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
#include <random>
#include <vector>
#include "ginterp.hpp"
#include "tbox/pio.hpp"
#include "boost/mpi.hpp"


template<typename T, std::size_t D>
T true_solution(const std::array<T,D>& a){
    T ans = 3 * a[0] * a[0];
    if constexpr (D >= 2){
        ans += 2 * a[1];
    }
    if constexpr (D >= 3){
        ans -= a[2];
    }
    return ans;
};

template <typename T, std::size_t D>
std::string array_to_string(const std::array<T, D>& a) {
    std::string ans;
    for (auto& val : a) {
        ans += std::to_string(val) + " ";
    }
    return ans;
};


int main(int argc, char **argv){
	namespace mpi  = boost::mpi;
	using namespace geoflow::tbox;
	constexpr int NDIM = GDIM;
	using size_type  = std::size_t;
	using value_type = double;
	using point_type = std::array<value_type, NDIM>;
    constexpr int num_targets = 10;
    constexpr int num_sources = 1000;
    constexpr int num_repeat  = num_sources;

	// Start Up MPI and find my place in it
    mpi::environment  env(argc,argv);
    mpi::communicator world;
	auto my_rank   = world.rank();
	auto num_ranks = world.size();
	pio::initialize(my_rank);
//	pio::logAllNodes("log");

    // Global Bounds
    std::array<value_type,NDIM> global_min_bound;
    std::array<value_type,NDIM> global_max_bound;
    global_min_bound.fill(-10.0);
    global_max_bound.fill(+10.0);

    // Rank Bounds
    // Make it a simple partition of X coordinates
    std::array<value_type,NDIM> rank_min_bound;
    std::array<value_type,NDIM> rank_max_bound;
    rank_min_bound = global_min_bound;
    rank_max_bound = global_max_bound;
    auto dx_per_rank = (global_max_bound[0] - global_min_bound[0]) / num_ranks;
    rank_min_bound[0] = global_min_bound[0] + my_rank * dx_per_rank;
    rank_max_bound[0] = global_min_bound[0] + (my_rank+1) * dx_per_rank;

    // Every rank gets random device
    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate Source Locations
    // Random within the full global domain
    std::vector<point_type> source_xyz(num_sources);
    for(size_type d = 0; d < NDIM; ++d){    
        std::uniform_real_distribution<value_type> dis(global_min_bound[d], global_max_bound[d]);
        for(size_type id = 0; id < num_sources; ++id){
            source_xyz[id][d] = dis(gen);
        }
    }

    // Ensure target region corners are included as sources for testing
    if constexpr ( NDIM == 1 ){
        source_xyz.push_back( {rank_min_bound[0]} );
        source_xyz.push_back( {rank_max_bound[0]} );
    }
    if constexpr ( NDIM == 2 ){
        source_xyz.push_back( {rank_min_bound[0],rank_min_bound[1]} );
        source_xyz.push_back( {rank_max_bound[0],rank_min_bound[1]} );
        source_xyz.push_back( {rank_max_bound[0],rank_max_bound[1]} );
        source_xyz.push_back( {rank_min_bound[0],rank_max_bound[1]} );
    }
    // if constexpr ( NDIM == 3 ){
    //     source_xyz.push_back( {rank_min_bound[0],rank_min_bound[1],rank_min_bound[2]} );
    //     source_xyz.push_back( {rank_max_bound[0],rank_min_bound[1],rank_min_bound[2]} );
    //     source_xyz.push_back( {rank_max_bound[0],rank_max_bound[1],rank_min_bound[2]} );
    //     source_xyz.push_back( {rank_min_bound[0],rank_max_bound[1],rank_min_bound[2]} );
    //     source_xyz.push_back( {rank_min_bound[0],rank_min_bound[1],rank_max_bound[2]} );
    //     source_xyz.push_back( {rank_max_bound[0],rank_min_bound[1],rank_max_bound[2]} );
    //     source_xyz.push_back( {rank_max_bound[0],rank_max_bound[1],rank_max_bound[2]} );
    //     source_xyz.push_back( {rank_min_bound[0],rank_max_bound[1],rank_max_bound[2]} );    
    // }

    // Create some Repeat Sources in data
    const auto sz = source_xyz.size();
    for(size_type i = 0; i < num_repeat; ++i){  
        source_xyz.push_back( source_xyz[ std::size_t(i * (sz-1) / num_repeat) ] );
    }

    // Generate Target Locations
    // Random within this ranks domain
    std::vector<point_type> target_xyz(num_targets);
    for(size_type d = 0; d < NDIM; ++d){    
        std::uniform_real_distribution<value_type> dis(rank_min_bound[d], rank_max_bound[d]);
        for(size_type id = 0; id < num_targets; ++id){
            target_xyz[id][d] = dis(gen);
        }
    }

    // Add a Source at 1st Target Location for Checking Exact Match
    source_xyz.push_back(target_xyz[0]);

    // Build Known Solution at Sources
    std::vector<value_type> source_soln(source_xyz.size(), 0);
    for(size_type id = 0; id < source_xyz.size(); ++id){
        source_soln[id] = true_solution(source_xyz[id]);
    }

    // Build Known Solution at Targets (for Checking Only)
    std::vector<value_type> target_exact_soln(target_xyz.size(), 0);
    for(size_type id = 0; id < target_xyz.size(); ++id){
        target_exact_soln[id] = true_solution(target_xyz[id]);
    }

    // Determine a tolerance to assume points are the same
    auto tolerance = dx_per_rank / 1000;

    // Init GInterp
    std::vector<value_type> target_soln(target_xyz.size());
	GInterp<value_type>::interpolate(tolerance,source_xyz,source_soln,target_xyz,target_soln);

	// Check Answer against known solution
	value_type error_tolerance = 0.01;
    auto percent_error = std::abs(target_soln[0] - target_exact_soln[0])/target_exact_soln[0];
	if( percent_error > error_tolerance ){
		pio::perr << "At Index = " << 0 << " Percent Error = " << percent_error << std::endl;
        pio::perr << "Returned Target Solution = " << target_soln[0] << std::endl;
        pio::perr << "Exact Target Solution    = " << target_exact_soln[0] << std::endl;
		return 1;
	}

	pio::finalize();
	return 0;
}
