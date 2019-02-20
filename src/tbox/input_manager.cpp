/*
 * input_manager.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#include "tbox/input_manager.hpp"
#include "tbox/mpixx.hpp"

namespace geoflow {
namespace tbox {

PropertyTree InputManager::ptree_;

void
InputManager::initialize(int argc, char* argv[]){
	std::string file_name = "input.jsn";
	if( argc == 2 ){
		file_name = argv[1];
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
