/*
 * enumerate_test.cpp
 *
 *  Created on: Jul 9, 2018
 *      Author: bryan.flynt
 */

#include "tbox/enumerate.hpp"

#include <iostream>
#include <vector>

template<typename T>
T calc(const T val){
	return val + 0.1*val;
}

using namespace geoflow::tbox;


int main() {

	const int N = 5;
	const int S = 2;
	int ans;

	std::vector<double> a(N);
	double b[N];

	// Assign values to a & b standard way
	for(int i = 0; i < N; ++i){
		a[i] = calc(i);
		b[i] = calc(i);
	}

	// Enumerate 0 -> 10
	ans = 0;
	for(auto i : enumerate(a)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 10;
		if( calc(indx) != std::get<1>(i) ) return ans + 20;
		++ans;
	}
	if( ans != N ) return 1000;

	// Enumerate 0 -> 10
	ans = 0;
	for(auto i : enumerate(b)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 30;
		if( calc(indx) != std::get<1>(i) ) return ans + 40;
		++ans;
	}
	if( ans != N ) return 2000;

	// Enumerate 0 -> 10
	ans = 1;
	for(auto i : enumerate(a.begin()+1, a.end(), 1)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 50;
		if( calc(indx) != std::get<1>(i) ) return ans + 60;
		++ans;
	}
	if( ans != N ) return 3000;



	// Enumerate 0 -> 10
	ans = 0;
	for(auto i : enumerate(a).step(S)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 70;
		if( calc(indx) != std::get<1>(i) ) return ans + 80;
		ans+=S;
	}
	if( ans < N ) return 4000;

	// Enumerate 0 -> 10
	ans = 0;
	for(auto i : enumerate(b).step(S)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 90;
		if( calc(indx) != std::get<1>(i) ) return ans + 100;
		ans+=S;
	}
	if( ans < N ) return 5000;

	// Enumerate 0 -> 10
	ans = 1;
	for(auto i : enumerate(a.begin()+1, a.end(), 1).step(S)){
		const int indx = std::get<0>(i);
		if( ans        != std::get<0>(i) ) return ans + 110;
		if( calc(indx) != std::get<1>(i) ) return ans + 120;
		ans+=S;
	}
	if( ans < N ) return 6000;

	return 0;
}
