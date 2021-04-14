/*
 * zip_test.cpp
 *
 *  Created on: Jul 9, 2018
 *      Author: bryan.flynt
 */

#include "tbox/zip.hpp"

#include <iostream>
#include <list>
#include <vector>


template<typename T>
T calc(const T val){
	return val + 0.1*val;
}

using namespace geoflow::tbox;

int main() {

	const int N = 10;
	int ans;

	std::vector<double> a;
	std::list<double>   b;
	std::vector<double> c;
	double d[N];

	// Assign values to a, b & c standard way
	for(int i = 0; i < N; ++i){
		a.push_back( calc(i) );
		b.push_back( calc(i) );
		c.push_back( calc(i) );
		d[i] = calc(i);
	}

	ans = 0;
	for(auto i : zip(a,b)){
		if( std::get<0>(i) != std::get<1>(i) ) return 100 + ans;
		++ans;
	}
	if( ans != N ) return 1000;

	ans = 0;
	for(auto i : zip(a,b,c)){
		if( std::get<0>(i) != std::get<1>(i) ) return 200 + ans;
		if( std::get<1>(i) != std::get<2>(i) ) return 300 + ans;
		++ans;
	}
	if( ans != N ) return 2000;


	return 0;
}
