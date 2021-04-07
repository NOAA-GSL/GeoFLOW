/*
 * rbf_interpolation.cpp
 *
 *  Created on: April 7, 2021
 *      Author: bflynt
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

#include "tbox/interpolation/rbf.hpp"

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
    using namespace geoflow::tbox;
	constexpr int NDIM = 3;
	using size_type  = std::size_t;
	using value_type = double;
	using point_type = std::array<value_type, NDIM>;
    constexpr int num_targets = 1;
    constexpr int num_sources = 50;

    // Global Bounds
    std::array<value_type,NDIM> min_bound;
    std::array<value_type,NDIM> max_bound;
    min_bound.fill(-10.0);
    max_bound.fill(+10.0);

    // Every rank gets random device
    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate Source Locations
    // Random within the full global domain
    std::vector<point_type> source_xyz(num_sources);
    for(size_type d = 0; d < NDIM; ++d){    
        std::uniform_real_distribution<value_type> dis(min_bound[d], max_bound[d]);
        for(size_type id = 0; id < num_sources; ++id){
            source_xyz[id][d] = dis(gen);
        }
    }

    // Ensure target region corners are included as sources for testing
    value_type eps = 0.001;
    if constexpr ( NDIM == 1 ){
        source_xyz.push_back( {min_bound[0]-eps} );
        source_xyz.push_back( {max_bound[0]+eps} );
    }
    if constexpr ( NDIM == 2 ){
        source_xyz.push_back( {min_bound[0]-eps,min_bound[1]-eps} );
        source_xyz.push_back( {max_bound[0]+eps,min_bound[1]-eps} );
        source_xyz.push_back( {max_bound[0]+eps,max_bound[1]+eps} );
        source_xyz.push_back( {min_bound[0]-eps,max_bound[1]+eps} );
    }
    if constexpr ( NDIM == 3 ){
        source_xyz.push_back( {min_bound[0]-eps,min_bound[1]-eps,min_bound[2]-eps} );
        source_xyz.push_back( {max_bound[0]+eps,min_bound[1]-eps,min_bound[2]-eps} );
        source_xyz.push_back( {max_bound[0]+eps,max_bound[1]+eps,min_bound[2]-eps} );
        source_xyz.push_back( {min_bound[0]-eps,max_bound[1]+eps,min_bound[2]-eps} );
        source_xyz.push_back( {min_bound[0]-eps,min_bound[1]-eps,max_bound[2]+eps} );
        source_xyz.push_back( {max_bound[0]+eps,min_bound[1]-eps,max_bound[2]+eps} );
        source_xyz.push_back( {max_bound[0]+eps,max_bound[1]+eps,max_bound[2]+eps} );
        source_xyz.push_back( {min_bound[0]-eps,max_bound[1]+eps,max_bound[2]+eps} );    
    }

    // Generate Target Locations
    std::vector<point_type> target_xyz(num_targets);
    for(size_type d = 0; d < NDIM; ++d){    
        std::uniform_real_distribution<value_type> dis(min_bound[d], max_bound[d]);
        for(size_type id = 0; id < num_targets; ++id){
            target_xyz[id][d] = dis(gen);
        }
    }

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

    // Init GInterp
    rbf::Interpolator<value_type> interpolator(rbf::kernel::Gaussian<value_type>);
    interpolator.setData(source_xyz, source_soln);
    interpolator.computeWeights();

    std::vector<value_type> target_soln(target_xyz.size(), 0);
    for(auto i = 0; i < target_soln.size(); ++i){
        target_soln[i] = interpolator.calcValue(target_xyz[i]);
    }

	// Check Answer against known solution
	value_type error_tolerance = 0.01;
    auto percent_error = std::abs(target_soln[0] - target_exact_soln[0])/target_exact_soln[0];
	//if( percent_error > error_tolerance ){
		std::cout << "At Index = " << 0 << " Percent Error = " << percent_error << std::endl;
        std::cout << "Returned Target Solution = " << target_soln[0] << std::endl;
        std::cout << "Exact Target Solution    = " << target_exact_soln[0] << std::endl;
		return 1;
	//}

	return 0;
}



