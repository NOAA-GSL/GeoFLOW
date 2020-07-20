/*
 * property_tree.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */


#include "tbox/property_tree.hpp"

#include <cassert>
#include <iostream>


using namespace geoflow::tbox;

int main() {

	PropertyTree pt;

	pt.load_file("input.jsn");

	auto keys = pt.getKeys();
	assert( keys.size() == 10 );

	assert( pt.isValue<int>("object_1")   == false );
	assert( pt.isArray<int>("object_1")   == false );
	assert( pt.isArray2D<int>("object_1") == false );
	assert( pt.isPropertyTree("object_1") == true  );

	auto obj1 = pt.getPropertyTree("object_1");
	assert( obj1.isValue<int>("age")      == true );
	assert( obj1.isArray<int>("age")      == false );
	assert( obj1.isArray2D<int>("age")    == false );
	assert( obj1.isPropertyTree("age")    == false  );

	//
	// Bool Values
	//
	assert( pt.isValue<bool>("bool_value")   == true  );
	assert( pt.isArray<bool>("bool_value")   == false );
	assert( pt.isArray2D<bool>("bool_value") == false );
	assert( pt.isPropertyTree("bool_value")  == false );

	assert( pt.isValue<bool>("bool_array")   == false );
	assert( pt.isArray<bool>("bool_array")   == true  );
	assert( pt.isArray2D<bool>("bool_array") == false );
	assert( pt.isPropertyTree("bool_array")  == false );
	auto bool_array = pt.getArray<bool>("bool_array");
	assert( bool_array.size() == 5 );
	assert( bool_array[0] == true );
	assert( bool_array[1] == false );
	assert( bool_array[2] == true );
	assert( bool_array[3] == false );
	assert( bool_array[4] == true );

	//
	// Integer Values
	//
	assert( pt.isValue<int>("int_value")   == true  );
	assert( pt.isArray<int>("int_value")   == false );
	assert( pt.isArray2D<int>("int_value") == false );
	assert( pt.isPropertyTree("int_value") == false );

	assert( pt.isValue<int>("int_array")   == false );
	assert( pt.isArray<int>("int_array")   == true  );
	assert( pt.isArray2D<int>("int_array") == false );
	assert( pt.isPropertyTree("int_array") == false );
	auto int_array = pt.getArray<int>("int_array");
	assert( int_array.size() == 5 );
	assert( int_array[0] == 0 );
	assert( int_array[1] == 1 );
	assert( int_array[2] == 2 );
	assert( int_array[3] == 3 );
	assert( int_array[4] == 4 );

	//
	// Real Values
	//
	assert( pt.isValue<double>("real_value")   == true  );
	assert( pt.isArray<double>("real_value")   == false );
	assert( pt.isArray2D<double>("real_value") == false );
	assert( pt.isPropertyTree("real_value")    == false );

	assert( pt.isValue<double>("real_array")   == false );
	assert( pt.isArray<double>("real_array")   == true  );
	assert( pt.isArray2D<double>("real_array") == false );
	assert( pt.isPropertyTree("real_array")    == false );
	auto real_array = pt.getArray<double>("real_array");
	assert( real_array.size() == 5 );
	assert( real_array[0] == 1.1 );
	assert( real_array[1] == 2.2 );
	assert( real_array[2] == 3.3 );
	assert( real_array[3] == 4.4 );
	assert( real_array[4] == 5.5 );

	//
	// String Values
	//
	assert( pt.isValue<std::string>("string_value")   == true  );
	assert( pt.isArray<std::string>("string_value")   == false );
	assert( pt.isArray2D<std::string>("string_value") == false );
	assert( pt.isPropertyTree("string_value")         == false );

	assert( pt.isValue<std::string>("string_array")   == false );
	assert( pt.isArray<std::string>("string_array")   == true  );
	assert( pt.isArray2D<std::string>("string_array") == false );
	assert( pt.isPropertyTree("string_array")         == false );
	auto string_array = pt.getArray<std::string>("string_array");
	assert( string_array.size() == 5 );
	assert( string_array[0] == "does" );
	assert( string_array[1] == "this" );
	assert( string_array[2] == "string" );
	assert( string_array[3] == "array" );
	assert( string_array[4] == "work");

	//
	// Missing Keys
	//
	assert( pt.keyExists("missing_key")      == false  );
	assert( pt.isValue<int>("missing_key")   == false  );
	assert( pt.isArray<int>("missing_key")   == false  );
	assert( pt.isArray2D<int>("missing_key") == false  );
	assert( pt.isPropertyTree("missing_key") == false  );

	//
	// 2-D Array Check
	//
	auto v2d = pt.getArray2D<int>("2D Array");
	for(int r = 0; r < 3; ++r){
		for(int c = 0; c < 3; ++c){
			assert( r+c == v2d[r][c] );
		}
	}

	return 0;
}
