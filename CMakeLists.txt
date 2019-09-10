# Check CMake Version 
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

#
#---------------------------------------------------------------------
# Load user requested variables if set (i.e. User Hardwired Values)
# NOTE: Must be done before "project" or "enable_language()"
#---------------------------------------------------------------------
#
set(USER_CONFIG "default" CACHE STRING "Default GeoFLOW Configuration File")
include("${CMAKE_SOURCE_DIR}/cmake/build.${USER_CONFIG}")

#
#---------------------------------------------------------------------
# Project Setup
#---------------------------------------------------------------------
#
project("GeoFLOW" VERSION 0.0.0.1 LANGUAGES C CXX Fortran)
#set(PROJECT_VERSION_MAJOR 0)
#set(PROJECT_VERSION_MINOR 0)
#set(PROJECT_VERSION_PATCH 0)
#set(PROJECT_VERSION_TWEAK 1)

#---------------------------------------------------------------------
# User Configure Build Process
#---------------------------------------------------------------------
set(USER_REQUEST_VERBOSE_CMAKE TRUE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_COLOR_MAKEFILE TRUE)
set(DEFAULT_BUILD_TYPE "Release")

#---------------------------------------------------------------------
# Set location of *.cmake modules
#---------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#
#---------------------------------------------------------------------
# Force build in seperate directory
#---------------------------------------------------------------------
#
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
	include(oos_error)   
	message(FATAL_ERROR)
endif()

#
#---------------------------------------------------------------------
# Directory Report
#---------------------------------------------------------------------
#
if(USER_REQUEST_VERBOSE_CMAKE)
	message("")
	message("--------------------- Directory Report -----------------------")
	message(" ") # CMAKE directories are from top level CMakeLists.txt
	message("Top Level Directories:")
	message("CMAKE_SOURCE_DIR  = ${CMAKE_SOURCE_DIR}")
	message("CMAKE_BINARY_DIR  = ${CMAKE_BINARY_DIR}")
	message(" ") # PROJECT directories are for recent project call
	message("Project Level Directories:")
	message("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
	message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
endif(USER_REQUEST_VERBOSE_CMAKE)

#---------------------------------------------------------------------
# Locations of Installation & Report
#---------------------------------------------------------------------

# Source locations in source tree
set(MY_PROJECT_PREFIX ${PROJECT_SOURCE_DIR})
set(MY_PROJECT_INCDIR ${MY_PROJECT_PREFIX}/src)
set(MY_PROJECT_SRCDIR ${MY_PROJECT_PREFIX}/src)
set(MY_PROJECT_TSTDIR ${MY_PROJECT_PREFIX}/test)

# Where to place all libraries
set(MY_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(MY_INSTALL_BINDIR ${MY_INSTALL_PREFIX}/bin)
set(MY_INSTALL_INCDIR ${MY_INSTALL_PREFIX}/include)
set(MY_INSTALL_SRCDIR ${MY_INSTALL_PREFIX}/src)
set(MY_INSTALL_LIBDIR ${MY_INSTALL_PREFIX}/lib)
set(MY_INSTALL_TSTDIR ${MY_INSTALL_PREFIX}/test)
set(MY_INSTALL_EXADIR ${MY_INSTALL_PREFIX}/examples)

if(USER_REQUEST_VERBOSE_CMAKE)
	message(" ")
	message("Using Directories:")
	message("--- Original Locations ---")
	message("MY_PROJECT_PREFIX = ${MY_PROJECT_PREFIX}")
	message("MY_PROJECT_INCDIR = ${MY_PROJECT_INCDIR}")
	message("MY_PROJECT_SRCDIR = ${MY_PROJECT_SRCDIR}")
	message("MY_PROJECT_TSTDIR = ${MY_PROJECT_TSTDIR}")
	message(" ")
	message("--- Installation Locations ---")
	message("MY_INSTALL_PREFIX = ${MY_INSTALL_PREFIX}")
	message("MY_INSTALL_BINDIR = ${MY_INSTALL_BINDIR}")
	message("MY_INSTALL_INCDIR = ${MY_INSTALL_INCDIR}")
	message("MY_INSTALL_SRCDIR = ${MY_INSTALL_SRCDIR}")
	message("MY_INSTALL_LIBDIR = ${MY_INSTALL_LIBDIR}")
	message("MY_INSTALL_TSTDIR = ${MY_INSTALL_TSTDIR}")
	message("MY_INSTALL_EXADIR = ${MY_INSTALL_EXADIR}")
endif(USER_REQUEST_VERBOSE_CMAKE)

#---------------------------------------------------------------------
# Detect OpenMP & MPI Flags/Includes/Libraries
#---------------------------------------------------------------------
message("")
message("------------------ Searching for Libraries ---------------------")
find_package(OpenMP)
find_package(MPI)

if(USER_REQUEST_VERBOSE_CMAKE)
	message("")
	message("--------------------- OpenMP Libraries ------------------------")
	message("OpenMP Found          = ${OpenMP_FOUND}")
	message("OpenMP F Comp Flags   = ${OpenMP_Fortran_FLAGS}")
	message("OpenMP F Inludes      = ${OpenMP_Fortran_INCLUDE_PATH}")
	message("OpenMP F Link Flags   = ${OpenMP_Fortran_LINK_FLAGS}")
	message("OpenMP F Libraries    = ${OpenMP_Fortran_LIBRARIES}")
	message("OpenMP C Comp Flags   = ${OpenMP_C_FLAGS}")
	message("OpenMP C Inludes      = ${OpenMP_C_INCLUDE_PATH}")
	message("OpenMP C Link Flags   = ${OpenMP_C_LINK_FLAGS}")
	message("OpenMP C Libraries    = ${OpenMP_C_LIBRARIES}")
	message("OpenMP C++ Comp Flags = ${OpenMP_CXX_FLAGS}")
	message("OpenMP C++ Inludes    = ${OpenMP_CXX_INCLUDE_PATH}")
	message("OpenMP C++ Link Flags = ${OpenMP_CXX_LINK_FLAGS}")
	message("OpenMP C++ Libraries  = ${OpenMP_CXX_LIBRARIES}")
	message("")
	message("------------------------ MPI Libraries --------------------------")
	message("MPI F Found        = ${MPI_Fortran_FOUND}")
	message("MPI F Comp Flags   = ${MPI_Fortran_FLAGS}")
	message("MPI F Inludes      = ${MPI_Fortran_INCLUDE_PATH}")
	message("MPI F Link Flags   = ${MPI_Fortran_LINK_FLAGS}")
	message("MPI F Libraries    = ${MPI_Fortran_LIBRARIES}")
	message("MPI C Found        = ${MPI_C_FOUND}")
	message("MPI C Comp Flags   = ${MPI_C_FLAGS}")
	message("MPI C Inludes      = ${MPI_C_INCLUDE_PATH}")
	message("MPI C Link Flags   = ${MPI_C_LINK_FLAGS}")
	message("MPI C Libraries    = ${MPI_C_LIBRARIES}")
	message("MPI C++ Found      = ${MPI_CXX_FOUND}")
	message("MPI C++ Comp Flags = ${MPI_CXX_FLAGS}")
	message("MPI C++ Inludes    = ${MPI_CXX_INCLUDE_PATH}")
	message("MPI C++ Link Flags = ${MPI_CXX_LINK_FLAGS}")
	message("MPI C++ Libraries  = ${MPI_CXX_LIBRARIES}")
endif(USER_REQUEST_VERBOSE_CMAKE)


#---------------------------------------------------------------------
# Detect Other Libraries
#---------------------------------------------------------------------
      
set(Boost_DEBUG OFF)              # Enable debug output from FIND_PACKAGE
set(Boost_NO_SYSTEM_PATHS ON)     # Do not search system paths before BOOST_ROOT
set(Boost_USE_STATIC_LIBS ON)     # Static link to Boost libraries
set(Boost_USE_STATIC_RUNTIME OFF) # Use Boost static linked to C++ runtime
set(Boost_USE_MULTITHREADED OFF)  # Use Boost multi-threaded code
FIND_PACKAGE(Boost 1.65.0 COMPONENTS mpi serialization REQUIRED)

if(USER_REQUEST_VERBOSE_CMAKE)
	message("")
	message("--------------------- Boost Libraries ------------------------")
	message("BOOST_ROOT = ${BOOST_ROOT}")
	message("Boost_NO_SYSTEM_PATHS = ${Boost_NO_SYSTEM_PATHS}")
	message("Boost Found               = ${Boost_FOUND}")
	message("Boost Version             = ${Boost_VERSION}")
	message("Boost Includes            = ${Boost_INCLUDE_DIRS}")
	message("Boost Link Libraries      = ${Boost_LIBRARY_DIRS}")
	message("Boost Component Libraries = ${Boost_LIBRARIES}")
endif(USER_REQUEST_VERBOSE_CMAKE)

#
#---------------------------------------------------------------------
# Compiler Feature Detection
#---------------------------------------------------------------------
# Hasn't proven to be very usefull. I find testing for the exact 
# feature of interest in the following section much better.
# Partialy because each compiler names features differently or 
# some not at all.
#
if(USER_REQUEST_VERBOSE_CMAKE)
	message("")
	message("------------------ Compiler Feature Detection ---------------------")
	message("")
	message("CMAKE_CXX_COMPILER         = ${CMAKE_CXX_COMPILER}")
	message("CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}")
	message("CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")

	message("")
	message("CMake Version Test can for the following C++ features:")
	get_property(known_features GLOBAL PROPERTY CMAKE_CXX_KNOWN_FEATURES)
	foreach(i ${known_features})
		message("${i}")
	endforeach()
	
	message("")
	message("C++ Compiler Supported Features:")
	foreach(i ${CMAKE_CXX_COMPILE_FEATURES})
		message("${i}")
	endforeach()
endif(USER_REQUEST_VERBOSE_CMAKE)

#---------------------------------------------------------------------
# Set Project Level Compiler Options
#---------------------------------------------------------------------
# ONLY those options which impact level of compiler errors
add_compile_options("-Werror")

#---------------------------------------------------------------------
# Check Platform Support
#---------------------------------------------------------------------

#message("")
#message("------------------------ Platform Support -------------------------")
#message("")

#include(cpp_chrono)   # Returns HAVE_CPP_CHRONO_H, C11_CHRONO_WORKS, C20_CHRONO_WORKS
#include(cpp_optional) # Returns HAVE_CPP_OPTIONAL_H, C17_OPTIONAL_WORKS
#include(cpp_variant)  # Returns HAVE_CPP_VARIANT_H, C17_VARIANT_WORKS
#include(posix_time)   # Returns POSIX_TIME_WORKS

#---------------------------------------------------------------------
# 
#---------------------------------------------------------------------

# Build Libraries & Applications
add_subdirectory(src)

# Build Executables
#add_subdirectory(src/apps)

# Build Tests for Libraries
add_subdirectory(test)

