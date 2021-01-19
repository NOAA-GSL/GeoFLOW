# Check CMake Version 
cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)

Crash CMAKE

#
#---------------------------------------------------------------------
# Load user requested variables if set (i.e. User Hardwired Values)
# NOTE: Must be done before "project" or "enable_language()"
#---------------------------------------------------------------------
# Default: 
# > cmake -D GCONF=default ..
#


set(GCONF "default" CACHE STRING "Default GeoFLOW Configuration File")
set(GCONF_FILE "${CMAKE_SOURCE_DIR}/cmake/build.${GCONF}")
message("Loading Configuration -> ${GCONF_FILE}")
include(${GCONF_FILE})

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
set(CMAKE_VERBOSE_MAKEFILE ${GEOFLOW_VERBOSE_CMAKE})
set(CMAKE_BUILD_TYPE ${GEOFLOW_BUILD_TYPE})
set(CMAKE_COLOR_MAKEFILE TRUE)

# Control message level of CMake
# Should be set on command line as follow by user
# --log-level=<ERROR|WARNING|NOTICE|STATUS|VERBOSE|DEBUG|TRACE>
# --loglevel=  Older versions of CMake
#set(CMAKE_MESSAGE_LOG_LEVEL STATUS)
#set(CMAKE_MESSAGE_LOG_LEVEL VERBOSE)

#---------------------------------------------------------------------
# Set location of *.cmake modules
#---------------------------------------------------------------------
#set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

#
#---------------------------------------------------------------------
# Force build in seperate directory
#---------------------------------------------------------------------
#
include(InSourceBuild)

#
#---------------------------------------------------------------------
# Directory Report
#---------------------------------------------------------------------
#
message(VERBOSE "")
message(VERBOSE "")
message(VERBOSE "--------------------- Directory Report -----------------------")
message(VERBOSE "") # CMAKE directories are from top level CMakeLists.txt
message(VERBOSE "Top Level Directories:")
message(VERBOSE "CMAKE_SOURCE_DIR  = ${CMAKE_SOURCE_DIR}")
message(VERBOSE "CMAKE_BINARY_DIR  = ${CMAKE_BINARY_DIR}")
message(VERBOSE "") # PROJECT directories are for recent project call
message(VERBOSE "Project Level Directories:")
message(VERBOSE "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
message(VERBOSE "PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")


#---------------------------------------------------------------------
# Locations of Installation & Report
#---------------------------------------------------------------------

# Source locations in source tree
set(MY_PROJECT_PREFIX ${PROJECT_SOURCE_DIR})
set(MY_PROJECT_INCDIR ${MY_PROJECT_PREFIX}/src)
set(MY_PROJECT_SRCDIR ${MY_PROJECT_PREFIX}/src)
set(MY_PROJECT_TSTDIR ${MY_PROJECT_PREFIX}/test)
set(MY_PROJECT_DOCDIR ${MY_PROJECT_PREFIX}/docs)

# Where to place all libraries
set(MY_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(MY_INSTALL_BINDIR ${MY_INSTALL_PREFIX}/bin)
set(MY_INSTALL_INCDIR ${MY_INSTALL_PREFIX}/include)
set(MY_INSTALL_SRCDIR ${MY_INSTALL_PREFIX}/src)
set(MY_INSTALL_LIBDIR ${MY_INSTALL_PREFIX}/lib)
set(MY_INSTALL_TSTDIR ${MY_INSTALL_PREFIX}/test)
set(MY_INSTALL_EXADIR ${MY_INSTALL_PREFIX}/examples)
set(MY_INSTALL_DOCDIR ${MY_INSTALL_PREFIX}/docs)

message(VERBOSE " ")
message(VERBOSE "Using Directories:")
message(VERBOSE "--- Original Locations ---")
message(VERBOSE "MY_PROJECT_PREFIX = ${MY_PROJECT_PREFIX}")
message(VERBOSE "MY_PROJECT_INCDIR = ${MY_PROJECT_INCDIR}")
message(VERBOSE "MY_PROJECT_SRCDIR = ${MY_PROJECT_SRCDIR}")
message(VERBOSE "MY_PROJECT_TSTDIR = ${MY_PROJECT_TSTDIR}")
message(VERBOSE " ")
message(VERBOSE "--- Installation Locations ---")
message(VERBOSE "MY_INSTALL_PREFIX = ${MY_INSTALL_PREFIX}")
message(VERBOSE "MY_INSTALL_BINDIR = ${MY_INSTALL_BINDIR}")
message(VERBOSE "MY_INSTALL_INCDIR = ${MY_INSTALL_INCDIR}")
message(VERBOSE "MY_INSTALL_SRCDIR = ${MY_INSTALL_SRCDIR}")
message(VERBOSE "MY_INSTALL_LIBDIR = ${MY_INSTALL_LIBDIR}")
message(VERBOSE "MY_INSTALL_TSTDIR = ${MY_INSTALL_TSTDIR}")
message(VERBOSE "MY_INSTALL_EXADIR = ${MY_INSTALL_EXADIR}")

#---------------------------------------------------------------------
# Detect Configuration for Record
#---------------------------------------------------------------------
include(RecordVars)

#---------------------------------------------------------------------
# Detect Library includes, flags, etc.
#---------------------------------------------------------------------
message(VERBOSE "")
message(VERBOSE "================================================================")
message(VERBOSE "                    Searching for Libraries                     ")
message(VERBOSE "================================================================")
if(${CMAKE_VERSION} VERSION_GREATER "3.11.999") 
	cmake_policy(SET CMP0074 NEW) # find_package search <name>_ROOT
endif()

if(GEOFLOW_USE_OPENMP) 
	find_package(_OpenMP REQUIRED) 
endif()

if(GEOFLOW_USE_MPI)
	find_package(_MPI REQUIRED)
endif()

find_package(_Boost REQUIRED)

if(GEOFLOW_TRACER_USE_GPTL)
	find_package(GPTL REQUIRED)
endif()

if(GEOFLOW_USE_PAPI)
	find_package(PAPI REQUIRED)
endif()

if(GEOFLOW_USE_CBLAS)
	find_package(CBLAS REQUIRED) 
endif()

if(GEOFLOW_USE_CUBLAS      OR 
   GEOFLOW_USE_OPENACC     OR
   GEOFLOW_TRACER_USE_NVTX)
	find_package(CUDAToolkit REQUIRED)
endif()

if(GEOFLOW_USE_OPENACC)
	find_package(OpenACC REQUIRED)
endif()

find_package(Doxygen)

#---------------------------------------------------------------------
# Detect Library Errors
#---------------------------------------------------------------------
# Protect against system which only install partial libraries
#

if(CUDAToolkit_FOUND)
	if(NOT TARGET CUDA::cudart)
		message( FATAL_ERROR "CUDA Runtime was Not Found within the CUDA Toolkit" )
	endif()
	if(GEOFLOW_USE_CUBLAS AND NOT TARGET CUDA::cublas)
		message( FATAL_ERROR "cuBLAS was Not Found within the CUDA Toolkit" )
	endif()
endif(CUDAToolkit_FOUND)

#
#---------------------------------------------------------------------
# Compiler Feature Detection
#---------------------------------------------------------------------
# Hasn't proven to be very usefull. I find testing for the exact 
# feature of interest in the following section much better.
# Partialy because each compiler names features differently or 
# some not at all.
#
message(VERBOSE "")
message(VERBOSE "------------------ Compiler Feature Detection ---------------------")
message(VERBOSE "")
message(VERBOSE "CMAKE_CXX_COMPILER         = ${CMAKE_CXX_COMPILER}")
message(VERBOSE "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}")
message(VERBOSE "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")

message(VERBOSE "")
message(VERBOSE "CMake Version Test can for the following C++ features:")
get_property(known_features GLOBAL PROPERTY CMAKE_CXX_KNOWN_FEATURES)
foreach(i ${known_features})
	message(VERBOSE "${i}")
endforeach()
	
message(VERBOSE "")
message(VERBOSE "C++ Compiler Supported Features:")
foreach(i ${CMAKE_CXX_COMPILE_FEATURES})
	message(VERBOSE "${i}")
endforeach()

#---------------------------------------------------------------------
# Set Project Level Compiler Options
#---------------------------------------------------------------------
# ONLY those options which impact level of compiler errors
# add_compile_options("-Werror")         # Change all Warnings to Errors
# add_compile_options("-Wfatal-errors")  # Stop at 1st Error Encountered
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

# Build analysis tools
add_subdirectory(analysis/src)

# Build Testing
if(GEOFLOW_BUILD_TESTS)
	include(CTest)
	if(BUILD_TESTING)
  		add_subdirectory(test)
	endif()
endif(GEOFLOW_BUILD_TESTS)

# Build Documentation (DOxygen)
if(DOXYGEN_FOUND)
	add_subdirectory(docs)
endif()


