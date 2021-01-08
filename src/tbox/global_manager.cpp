/*
 * global_manager.cpp
 *
 *  Created on: Nov 14, 2018
 *      Author: bflynt
 */

#include "tbox/global_manager.hpp"

#include "tbox/assert.hpp"
#include "tbox/error_handler.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/pio.hpp"

namespace geoflow {
namespace tbox {

bool GlobalManager::s_initialized = false;
bool GlobalManager::s_started     = false;

//
// Initialize the package.
// This routine performs the following tasks:
//
// (1) Sets abort handler to call MPI abort
// (2) Set the new abort to call MPI abort
// (3) Initialize the parallel I/O routines
// (4) Set new handler so that an error message is printed if new fails.
// (5) Triggers the input manager to read the program json file
// (6) Calls initialize on all routines registered with Manager
//
void GlobalManager::initialize(int argc, char* argv[]){
   ASSERT(!s_initialized);

	// Set ErrorHandler to use MPI
	EH::setAbortHandler(&mpixx_abort_handler);

	// Set new error handler
	std::set_new_handler(EH::badNew);

	// Start of PIO
	pio::initialize(mpixx::communicator().rank());
	//pio::logOnlyNodeZero("log");
	//pio::logAllNodes("log");

	// Read Input File
	InputManager::initialize(argc,argv);

	// Initialize all call backs
   //StartupShutdownManager::initialize();

   s_initialized = true;
}

//
// StartUp the package.
// This routine calls startup on all
// routines registered with the Manager
//
void GlobalManager::startup(){
   ASSERT(s_initialized);
   ASSERT(!s_started);
   //StartupShutdownManager::startup();
   s_started = true;
}

//
// Shutdown the package.
// This routine calls shutdown on all
// routines registered with the Manager
//
void GlobalManager::shutdown(){
   ASSERT(s_initialized);
   ASSERT(s_started);
   //StartupShutdownManager::shutdown();
   s_started = false;
}

//
// Finalize the package.
// This routine performs the following tasks:
//
// (1) Set the new abort to call MPI abort
// (2) Calls finalize on all routines registered with Manager
//
void GlobalManager::finalize(){
   ASSERT(s_initialized);

   // Finalize all call backs
   //StartupShutdownManager::finalize();

   // Finalize PIO Buffers
   pio::finalize();

   s_initialized = false;
}

bool GlobalManager::isInitialized(){
   return s_initialized;
}

bool GlobalManager::isStarted(){
   return s_started;
}

} // namespace tbox
} // namespace geoflow

