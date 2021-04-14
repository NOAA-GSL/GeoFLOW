/*
 * morton.cpp
 *
 *  Created on: Nov 7, 2019
 *      Author: bflynt
 */

#include "tbox/morton.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>

using namespace geoflow::tbox;

int main() {

	const int ndims = 3;
	const int nbits = 3 * ndims;

	{ // Array Test

		// Create MortonIndexer with limits
		using Coordinate = std::array<double,3>;
		Coordinate min_box = {-10.0, -10.0, -10.0};
		Coordinate max_box = {+10.0, +10.0, +10.0};
		MortonIndexer<Coordinate> morton(min_box,max_box);

		// Create location and index
		Coordinate a = {-1.0, -1.0, -1.0};
		MortonIndex<nbits> ia;
		morton.generate(a,ia);

		// Create location and index
		Coordinate b = {+1.0, +1.0, +1.0};
		MortonIndex<nbits> ib;
		morton.generate(b,ib);

		std::vector<MortonIndex<nbits>> indexes = {ib,ia};

		std::sort(indexes.begin(), indexes.end());
		for(const auto& val: indexes){
			std::cout << val.to_string() << std::endl;
		}
	}

	{ // Vector Test

		// Create MortonIndexer
		using Coordinate = std::vector<double>;
		Coordinate min_box = {-10.0, -10.0, -10.0};
		Coordinate max_box = {+10.0, +10.0, +10.0};
		MortonIndexer<Coordinate> morton(min_box,max_box);

		// Create location and index
		Coordinate a = {-1.0, -1.0, -1.0};
		MortonIndex<nbits> ia;
		morton.generate(a,ia);

		// Create location and index
		Coordinate b = {+1.0, +1.0, +1.0};
		MortonIndex<nbits> ib;
		morton.generate(b,ib);

		std::vector<MortonIndex<nbits>> indexes = {ib,ia};

		std::sort(indexes.begin(), indexes.end());
		for(const auto& val: indexes){
			std::cout << val.to_string() << std::endl;
		}
	}


	return 0;
}
