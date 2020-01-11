#
# Find OpenMP Compiler Features
# 
# Find and configure for OpenMP usage
#
# Warning: 
# CMake < 3.4 has a bug in the Threads package that requires 
# you to have the C language enabled.
#

if(USE_OPENMP)
message("")
message("--------------------- OpenMP Libraries ------------------------")

#
# First attempt to use the CMake tools to locate the library
# - Disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
# - Use CMake built in script
# - ReAdd ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
#
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
find_package(OpenMP REQUIRED)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
	
#
# For CMake < 3.9, we need to make the target ourselves
#
if(NOT TARGET OpenMP::OpenMP_C)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_C IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_C
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_C_FLAGS})
    set_property(TARGET OpenMP::OpenMP_C
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_C_FLAGS} Threads::Threads)
endif()

if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

if(NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_Fortran IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_Fortran
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_Fortran_FLAGS})
    set_property(TARGET OpenMP::OpenMP_Fortran
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_Fortran_FLAGS} Threads::Threads)
endif()

#
# Report what we found
#
if(CMAKE_VERBOSE_MAKEFILE)
message("")
message("OpenMP Found          = ${OpenMP_FOUND}")
message("OpenMP Version        = ${OpenMP_VERSION}")
message("")
message("C")
message(" - Found      = ${OpenMP_C_FOUND}")
message(" - Comp Flags = ${OpenMP_C_FLAGS}")
message(" - Lib Names  = ${OpenMP_C_LIB_NAMES}")
message(" - Libraries  = ${OpenMP_C_LIBRARIES}")
message("C++")
message(" - Found      = ${OpenMP_CXX_FOUND}")
message(" - Comp Flags = ${OpenMP_CXX_FLAGS}")
message(" - Lib Names  = ${OpenMP_CXX_LIB_NAMES}")
message(" - Libraries  = ${OpenMP_CXX_LIBRARIES}")
message("Fortran:")
message(" - Found      = ${OpenMP_Fortran_FOUND}")
message(" - Comp Flags = ${OpenMP_Fortran_FLAGS}")
message(" - Lib Names  = ${OpenMP_Fortran_LIB_NAMES}")
message(" - Libraries  = ${OpenMP_Fortran_LIBRARIES}")
message(" - omp_lib.h  = ${OpenMP_Fortran_HAVE_OMPLIB_HEADER}")
message(" - omp_lib    = ${OpenMP_Fortran_HAVE_OMPLIB_MODULE}")
endif(CMAKE_VERBOSE_MAKEFILE)

endif(USE_OPENMP)

