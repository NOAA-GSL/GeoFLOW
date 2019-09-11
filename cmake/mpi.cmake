#
# Check for MPI
#

set(CMAKE_REQUIRED_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_REQUIRED_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})

check_cxx_source_compiles("
#include "mpi.h"
int main(){

	MPI_Init(NULL, NULL);
	
	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    MPI_Finalize();

  return 0;
}
"
MPI_WORKS) 
