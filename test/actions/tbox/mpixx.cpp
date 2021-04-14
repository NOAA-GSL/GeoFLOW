/*
 * mpixx.cpp
 *
 *  Created on: Nov 13, 2018
 *      Author: bflynt
 */

#include "tbox/mpixx.hpp"

#include <cassert>

int main() {

	mpixx::environment env;
	mpixx::communicator world;

	if (world.rank() == 0) {
		world.send(1, 0, std::string("Hello"));
	    std::string msg;
	    world.recv(1, 1, msg);

	    assert( msg == "world" );
	} else if (world.rank() == 1) {
	    std::string msg;
	    world.recv(0, 0, msg);
	    world.send(0, 1, std::string("world"));
	    assert( msg == "Hello" );
	}

	return 0;
}


