/*
 * indices_tests.cpp
 *
 *  Created on: Jul 9, 2018
 *      Author: bryan.flynt
 */


#include "tbox/indices.hpp"

#include <iostream>
#include <vector>

using namespace geoflow::tbox;


int main() {

	const int N = 10;
	int ans;

	std::vector<double> a(N);
	double b[N];

	// Indices 0 -> 10
	ans = 0;
	for(auto i : indices(a)){
		if( ans != i ) return ans+1;
		++ans;
	}
	if( ans != N ) return 1234;

	// Indices 0 -> 10
	ans = 0;
	for(auto i : indices(b)){
		if( ans != i ) return ans+10;
		++ans;
	}
	if( ans != N ) return 2345;

	// Range 0 -> Infinity
	ans = 0;
	for(auto i : indices({12,13,15,10})){
		if( ans != i ) return ans+100;
		++ans;
	}
	if( ans != 4 ) return 3456;

	return 0;
}

