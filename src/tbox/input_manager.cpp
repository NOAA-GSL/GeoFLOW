/*
 * input_manager.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#include "tbox/command_line.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/pio.hpp"

namespace geoflow {
namespace tbox {

PropertyTree InputManager::ptree_;

//
// TODO: Remove the -h requirement
// It breaks encapsulation for InputManager to check for other
// classes command line flags.  An easy fix would be to require
// all runs to use the -i flag to specify the file.  If no -i is
// used then no file is read.
//
void
InputManager::initialize(int argc, char* argv[]){
	CommandLine cline(argc,argv);
	if( cline.tagExists("-h") ){
		pio::pout << "IGNORE INPUT FILE: '-h' found on command line\n";
		return; // "-h" tag negates any file reading
	}

	std::string file_name = "input.jsn";
	if( cline.tagExists("-i") ){
		file_name = cline.getValue("-i");
	}
	InputManager::loadInputFile(file_name);
}

PropertyTree
InputManager::getInputPropertyTree(){
	return ptree_;
}

void
InputManager::loadInputFile(const std::string& file_name){
	mpixx::communicator world;
	if( world.rank() == 0 ){
		ptree_.load_file(file_name);
	}
	if(world.size() > 1){
		mpixx::broadcast(world,ptree_,0);
	}
}

} // namespace tbox
} // namespace geoflow
