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
CommandLine  InputManager::cline_;

//
// TODO: Remove the -h requirement
// It breaks encapsulation for InputManager to check for other
// classes command line flags.  An easy fix would be to require
// all runs to use the -i flag to specify the file.  If no -i is
// used then no file is read.
//
void
InputManager::initialize(int argc, char* argv[]){

    // Process command lines
    cline_.process(argc,argv);

    // Read input file if requested
    if( not cline_.exists("h","help") ) {
        auto file_name = cline_.get("i","input","input.jsn");
        InputManager::loadInputFile(file_name);
    }
}

PropertyTree
InputManager::getInputPropertyTree(){
	return ptree_;
}

CommandLine InputManager::getInputCommandLine(){
    return cline_;
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
