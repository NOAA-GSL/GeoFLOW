#
# Find MPI Compiler & Features
#
# This wrapper uses the CMake provided function for
# locating the MPI compiler and options. After the 
# completion of that operation it will check for the 
# modern MPI::MPI_C, MPI::MPI_CXX, MPI::MPI_Fortran 
# targets and set them if they are not provided.
#
# Calling the builtin CMake function follows the 
# four step process described below.
#
# 1. Check if existing compiler CC, CXX, FC supports MPI.
# 2. Search for environment variables MPI_<lang>_COMPILER
# 3. Attempt to find an MPI compiler wrapper on system
# 4. Try to find an MPI that does not ship (Windows)
# Where: <lang> can be C, CXX, Fortran    
#
# If all goes as planned MPI code can be compiled 
# and linked using the CMake Properties 
#
# target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)
#

if(USE_MPI)
message(STATUS "")
message(STATUS "--------------------- MPI Libraries ------------------------")

#
# Determine if existing CC, CXX, FC compiler is capable of compiling 
# MPI code
#
include(CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES(
"
#include <mpi.h>
int main() {
    MPI_Init(NULL, NULL);
    MPI_Finalize();
    return 0;
}
"
CC_IS_MPI_WRAPPER)

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
"
#include <mpi.h>
int main() {
    MPI_Init(NULL, NULL);
    MPI_Finalize();
    return 0;
}
"
CXX_IS_MPI_WRAPPER)

include(CheckFortranSourceCompiles)
check_fortran_source_compiles(
"
integer :: ierr
include 'mpif.h'
call mpi_init(ierr)
call mpi_finalize(ierr)
end
"
FC_IS_WRAPPER SRC_EXT F90)

#
# First attempt to use the CMake tools to locate the library
#
find_package(MPI REQUIRED)


# For supporting CMake < 3.9:
if(NOT TARGET MPI::MPI_C)
	add_library(MPI::MPI_C IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_C_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES})
endif()
if(NOT TARGET MPI::MPI_CXX)
	add_library(MPI::MPI_CXX IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif()
if(NOT TARGET MPI::MPI_Fortran)
	add_library(MPI::MPI_Fortran IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_Fortran_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_Fortran_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
endif()


#
# Report what we found
#
if(CMAKE_VERBOSE_MAKEFILE)
message(VERBOSE "")
message(VERBOSE "MPI Found          = ${MPI_FOUND}")
message(VERBOSE "MPI Version        = ${MPI_VERSION}")
message(VERBOSE "")
message(VERBOSE "C")
message(VERBOSE " - Found       = ${MPI_C_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_C_COMPILER}")
message(VERBOSE " - Options     = ${MPI_C_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_C_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_C_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_C_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_C_LIBRARIES}")
message(VERBOSE "C++")
message(VERBOSE " - Found       = ${MPI_CXX_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_CXX_COMPILER}")
message(VERBOSE " - Options     = ${MPI_CXX_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_CXX_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_CXX_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_CXX_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_CXX_LIBRARIES}")
message(VERBOSE "Fortran:")
message(VERBOSE " - Found       = ${MPI_Fortran_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_Fortran_COMPILER}")
message(VERBOSE " - Options     = ${MPI_Fortran_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_Fortran_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_Fortran_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_Fortran_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_Fortran_LIBRARIES}")
endif(CMAKE_VERBOSE_MAKEFILE)

endif(USE_MPI)