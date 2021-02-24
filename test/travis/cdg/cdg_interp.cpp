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
// #include "ginterp.hpp"
#include "tbox/pio.hpp"
#include "boost/mpi.hpp"


int main(int argc, char **argv){
	namespace mpi  = boost::mpi;
	using namespace geoflow::tbox;
	constexpr int NDIM = GDIM;
	using size_type  = std::size_t;
	using value_type = double;
	using point_type = std::array<value_type, NDIM>;
    constexpr int num_targets = 1000;
    constexpr int num_sources = 1000;
    constexpr int num_repeat  = 100;

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
    std::vector<point_type> source_xyz(num_sources+num_repeat);
    for(size_type d = 0; d < NDIM; ++d){    
        std::uniform_real_distribution<value_type> dis(global_min_bound[d], global_max_bound[d]);
        for(size_type id = 0; id < num_sources; ++id){
            source_xyz[id][d] = dis(gen);
        }
    }

    // Create some Repeat Sources in data
    for(size_type i = 0; i < num_repeat; ++i){  
        source_xyz[num_sources + i] = source_xyz[ std::size_t(i * (num_sources-1) / num_repeat) ];
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

    // Build Known Solution at Sources
    std::vector<value_type> source_soln(source_xyz.size(), 0);
    for(size_type id = 0; id < source_xyz.size(); ++id){
        for(size_type d = 0; d < NDIM; ++d){    
            source_soln[id] += d * source_xyz[id][d];
        }
    }

    // Build Known Solution at Targets (for Checking Only)
    std::vector<value_type> target_exact_soln(target_xyz.size(), 0);
    for(size_type id = 0; id < target_xyz.size(); ++id){
        for(size_type d = 0; d < NDIM; ++d){    
            target_exact_soln[id] += d * target_xyz[id][d];
        }
    }

	// Determine a tolerance to assume points are the same
	auto tolerance = dx_per_rank / 1000;

	// Init GInterp
	// GInterp interp;
	// interp.init(tolerance,source_xyz,target_xyz);

	// // Perform Interpolations using solutions at sources
	// auto target_soln = interp.interpolate(source_soln);

	// // Check Answer against known solution
	// value_type error_tolerance = 0.01;
	// for(size_type i = 0; i < num_targets; ++i) {
	// 	auto percent_error = std::abs(target_soln[i] - target_exact_soln[i])/target_exact_soln[i];
	// 	if( percent_error > error_tolerance ){
	// 		pio::pout << "Index = " << i << " %Error = " << percent_error << std::endl;
	// 		return 1;
	// 	}
	// }

	pio::finalize();
	return 0;
}
