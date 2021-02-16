//
// Test for ASSERT functionality
//
// The test is designed to fail the assertion so it should not 
// be included into continuous integration testing
//

#include "tbox/assert.hpp"
#include "tbox/pio.hpp"
#include "tbox/property_tree.hpp"
#include "tbox/mpixx.hpp"
#include "tbox/global_manager.hpp"
#include "tbox/input_manager.hpp"
#include "tbox/tracer.hpp"

using namespace geoflow::tbox;
using namespace std;

int main(int argc, char **argv) {
  GEOFLOW_TRACE_INITIALIZE(); // Must be before MPI_Init (thanks GPTL)
  using namespace ::geoflow::tbox;

    // Initialize comm & global environment:
    mpixx::environment  env(argc,argv); // init GeoFLOW comm
    mpixx::communicator world;
    GlobalManager::initialize(argc,argv);
    GlobalManager::startup();
    GEOFLOW_TRACE();  // Must be after MPI_Init

    // Pass ASSERT
    ASSERT("I will Pass");
    ASSERT_MSG("I will Pass", "Its great to be " << "good");
    
    // Failure
    //ASSERT_MSG(false);
    ASSERT_MSG(false, "Failure is not an option");

     //***************************************************
    // Do shutdown, cleaning:
    //***************************************************
    pio::pout << "assert: do shutdown..." << std::endl;
    pio::finalize();
    GlobalManager::shutdown();
    GlobalManager::finalize();
    GEOFLOW_TRACE_STOP();      // GPTL requires popping main() off the stack
    GEOFLOW_TRACE_FINALIZE();  // Must be after MPI_Finalize (thanks GPTL)
    return(0);
}