#
# Find OpenMP Compiler Flags
#
# This wrapper uses the CMake provided function for
# locating the OpenMP options. The following 
# details are an short description of how to use it.
#
# Languages <lang>
# =======================
# C       = C Compiler Support
# CXX     = C++ Compiler Support
# Fortran = Fortran Compiler Support
#
# Imported Targets
# =======================
#
# OpenMP::OpenMP_<lang> = Target for using OpenMP from <lang>
#
# Result Variables
# =======================
#
# OpenMP_FOUND            = OpenMP flags for all requested languages have been found
# OpenMP_<lang>_FOUND     = OpenMP settings for <lang> were found
# OpenMP_<lang>_FLAGS     = OpenMP compiler flags for <lang>.
# MPI_<lang>_INCLUDE_DIRS = Include path(s) for OpenMP header
# OpenMP_<lang>_LIB_NAMES = List of libraries for OpenMP programs for <lang>
# MPI_<lang>_LIBRARIES    = All libraries to link OpenMP programs against
# MPI_<lang>_VERSION      = MPI version implemented for <lang>
#
# OpenMP_Fortran_HAVE_OMPLIB_HEADER = OpenMP is accessible through "omp_lib.h"
# OpenMP_Fortran_HAVE_OMPLIB_MODULE = OpenMP is accessible through the "omp_lib" module
#
# Variables To Locate MPI
# =======================
#
# OpenMP_<lang>_INCLUDE_DIR  = Header search path to find the relevant OpenMP headers
#
# Warning
# =======================
# CMake < 3.4 has a bug in the Threads package that requires 
# you to have the C language enabled.
#

message(VERBOSE "")
message(VERBOSE "--------------------- OpenMP Libraries ------------------------")
#
# First attempt to use the CMake tools to locate the library
#
find_package(OpenMP REQUIRED)

#
# For CMake < 3.9, we need to make the target ourselves
#
if(OpenMP_C_FOUND AND NOT TARGET OpenMP::OpenMP_C)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_C IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_C
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_C_FLAGS})
    set_property(TARGET OpenMP::OpenMP_C
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_C_FLAGS} Threads::Threads)
endif()

if(OpenMP_CXX_FOUND AND NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

if(OpenMP_Fortran_FOUND AND NOT TARGET OpenMP::OpenMP_Fortran)
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

