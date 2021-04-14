/*
 * range_test.cpp
 *
 *  Created on: Jul 9, 2018
 *      Author: bryan.flynt
 */

#include "tbox/range.hpp"

#include <iostream>
#include <vector>


using namespace geoflow::tbox;


int main() {
	int ans;

	// Range 0 -> 10
	ans = 0;
	for(auto i : range(0,10)){
		if( ans != i ) return ans+1;
		++ans;
	}
	if( ans != 10 ) return 1234;

	// Range 2 -> 10
	ans = 2;
	for(auto i : range(2,10)){
		if( ans != i ) return ans;
		++ans;
	}
	if( ans != 10 ) return 1235;

	// Range 0 -> Infinity
	ans = 0;
	for(auto i : range(0)){
		if( 31 == i ) break;
		if( ans != i ) return ans+1;
		++ans;
	}

	// Range 3 -> Infinity
	ans = 3;
	for(auto i : range(3)){
		if( 31 == i ) break;
		if( ans != i ) return ans;
		++ans;
	}

	// Range 0 -> 10 by 2
	ans = 0;
	for(auto i : range(0,10).step(2)){
		if( ans != i ) return ans+1;
		ans+=2;
	}
	if( ans < 10 ) return 1236;

	// Range 2 -> 10
	ans = 2;
	for(auto i : range(2,10).step(2)){
		if( ans != i ) return ans;
		ans+=2;
	}
	if( ans < 10 ) return 1237;

	// Range 0 -> Infinity by 2
	ans = 0;
	for(auto i : range(0).step(2)){
		if( 31 <= i ) break;
		if( ans != i ) return ans+1;
		ans+=2;
	}

	// Range 3 -> Infinity by 2
	ans = 3;
	for(auto i : range(3).step(2)){
		if( 31 <= i ) break;
		if( ans != i ) return ans;
		ans+=2;
	}

	return 0;
}
