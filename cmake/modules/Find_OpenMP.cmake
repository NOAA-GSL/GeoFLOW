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
message(STATUS "")
message(STATUS "--------------------- OpenMP Libraries ------------------------")

#
# First attempt to use the CMake tools to locate the library
#
find_package(OpenMP REQUIRED)

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
message(VERBOSE "")
message(VERBOSE "OpenMP Found          = ${OpenMP_FOUND}")
message(VERBOSE "")
message(VERBOSE "C")
message(VERBOSE " - Found      = ${OpenMP_C_FOUND}")
message(VERBOSE " - Version    = ${OpenMP_C_VERSION}")
message(VERBOSE " - Comp Flags = ${OpenMP_C_FLAGS}")
message(VERBOSE " - Lib Names  = ${OpenMP_C_LIB_NAMES}")
message(VERBOSE " - Libraries  = ${OpenMP_C_LIBRARIES}")
message(VERBOSE "C++")
message(VERBOSE " - Found      = ${OpenMP_CXX_FOUND}")
message(VERBOSE " - Version    = ${OpenMP_CXX_VERSION}")
message(VERBOSE " - Comp Flags = ${OpenMP_CXX_FLAGS}")
message(VERBOSE " - Lib Names  = ${OpenMP_CXX_LIB_NAMES}")
message(VERBOSE " - Libraries  = ${OpenMP_CXX_LIBRARIES}")
message(VERBOSE "Fortran:")
message(VERBOSE " - Found      = ${OpenMP_Fortran_FOUND}")
message(VERBOSE " - Version    = ${OpenMP_Fortran_VERSION}")
message(VERBOSE " - Comp Flags = ${OpenMP_Fortran_FLAGS}")
message(VERBOSE " - Lib Names  = ${OpenMP_Fortran_LIB_NAMES}")
message(VERBOSE " - Libraries  = ${OpenMP_Fortran_LIBRARIES}")
message(VERBOSE " - omp_lib.h  = ${OpenMP_Fortran_HAVE_OMPLIB_HEADER}")
message(VERBOSE " - omp_lib    = ${OpenMP_Fortran_HAVE_OMPLIB_MODULE}")

endif(USE_OPENMP)

