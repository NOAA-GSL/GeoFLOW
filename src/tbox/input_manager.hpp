/*
 * input_manager.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#ifndef SRC_GEOFLOW_TBOX_INPUT_MANAGER_HPP_
#define SRC_GEOFLOW_TBOX_INPUT_MANAGER_HPP_

#include <string>

#include "tbox/property_tree.hpp"
#include "tbox/command_line.hpp"


namespace geoflow {
namespace tbox {

class PropertTree;

/**
 * Class InputManager parses an input file and returns the associated
 * PropertyTree.  This manager class hides the complexity of opening the
 * input file, creating the parser, populating the tree with values and
 * broadcasting to all processors.
 *
 * All processors must call the parsing routine.  Any errors are reported
 * to pout and will result in termination of the program.
 *
 */
class InputManager{
public:

	/**
	 * @brief Initial setup of the InputManager.
	 *
	 * This function should be invoked ONLY ONCE at the start of a process
	 * to initialize and AFTER MPI is initialized (if used) by a
	 * call to one of the MPI init routines.
	 *
	 * This function attempts to read the provided file as
	 * program input for parsing into a PropertyTree.
	 */
	static void initialize(int argc, char* argv[]);

	/**
	 * Access method for the root input property tree held by InputManager.
	 * Inputs are read from the input file and held in this tree.
	 * This method returns the PropertyTree, allowing any class
	 * to access the information inside it using standard PropertyTree calls.
	 * For example, the following could appear in a class:
	 *
	 *       // get root database
	 *       PropertyTree root_pt = InputManager::getInputPropertyTree();
	 *
	 *       // get class's sub-tree
	 *       PropertyTree class_db = root_db.getPropertyTree("MyClass");
	 *       // get parameter(s) from sub-database
	 *       int dummy = class_db.getValue<int>("dummy");
	 *
	 * where "dummy" was supplied in "MyClass" entry of the input file:
	 *
	 *       MyClass {
	 *          dummy = ...
	 *       }
	 *
	 * This function is intended for classes in which there is no
	 * easy or efficient way to supply input parameters.
	 */
	static PropertyTree getInputPropertyTree();

	/**
	 * Access method for the parsed command line arguments held by the InputManager.
	 * Inputs are read from the command line and held within the map.
	 * This method returns the CommandLine, allowing any class
	 * to access the information inside it using the CommandLine interface.
	 * For example, the following could appear in a class:
	 *
	 *  // Get Command Line args
	 *  auto command_line = InputManager::getInputCommandLine();
	 *
	 *  // Test if "-b" or "--bench" was passed
	 *  bool bench_requested = command_line.exists("b","bench");
	 *
	 *  // Get -i or --input with default argument
	 *  std::string name = command_line.get("i","input","def_arg");
	 *
	 */
	static CommandLine getInputCommandLine();


	InputManager(const InputManager&)    = delete; // Don't allow copy
	void operator=(const InputManager&)  = delete; // Don't allow assign

private:

	static PropertyTree ptree_;
	static CommandLine  cline_;

	// Load data from the specified input file.
	static void loadInputFile(const std::string& filename);

};

} // namespace tbox
} // namespace geoflow

#endif /* SRC_GEOFLOW_TBOX_INPUT_MANAGER_HPP_ */
