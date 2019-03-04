#==================================================================================
# Module       : CMakeLists.txt
# Date         : 1/31/18 (DLR)
# Description  : GeoFIRE test cmake file
# Copyright    : Copyright 2018. Colorado State University. All rights reserved
#==================================================================================
#
cmake_minimum_required(VERSION 2.8)

# Following may be changed by invoking
#   > cmake -DCMAKE_C_COMPILER=pgcc 
# etc..  These 'sets' must be done before 'project'
# is specified:
 set(CMAKE_C_COMPILER       mpiicc)
 set(CMAKE_CXX_COMPILER     mpiicpc)
 set(CMAKE_Fortran_COMPILER mpiifort)

# Following may be used to instrument entire code
# using TAU:
#set(THOME /home/Duane.Rosenberg/lib/tau-2.27/x86_64)
#set(CMAKE_C_COMPILER       ${THOME}/bin/taucc -tau:makefile=${THOME}/lib/Makefile.tau-papi)
#set(CMAKE_CXX_COMPILER     ${THOME}/bin/taucxx -tau:makefile=${THOME}/lib/Makefile.tau-papi)
#set(CMAKE_Fortran_COMPILER ${THOME}/bin/tauf90 -tau:makefile=${THOME}/lib/Makefile.tau-papi)

project(GTest)

# The version number.
set(GTest_VERSION_MAJOR 1)
#set(GTest_VERSION_MINOR 0)
message("CMAKE_VERSION: ${CMAKE_VERSION}")

enable_language (Fortran)

# Load in project macros/functions:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gtest_macros.cmake)

# Define and set default CMake logical variables here:
option(DO_DEBUG        "DO_DEBUG"     ON)   # Compile for DEBUG
option(DO_OPENMP       "DO_OPENMP"    ON)   # Compile for OpenMP
option(DO_OPENACC      "DO_OPENACC"  OFF)   # Compile for OpenACC
option(USE_GBLAS       "USE_GBLAS"    ON)   # Use GBlas instead of C 
option(USE_MPI         "USE_MPI"      ON)   # Build with MPI?
option(HAVE_PAPI       "HAVE_PAPI"   OFF)   # Is GPTL built with PAPI support?
option(DO_AUTO_PROF "DO_AUTO_PROF"   OFF)   # Do auto-profing with GPTL?
 
# Following should be set in call to cmake as
# cmake -DCMAKE_BUILD_TYPE=DEBUG [or RELEASE]:
set(CMAKE_BUILD_TYPE DEBUG) 
if (NOT DO_DEBUG)
  set(CMAKE_BUILD_TYPE RELEASE) 
  message("CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
endif()

set (CMAKE_CXX_STANDARD 11) # good for only cmake v3+
#add_subdirectory(../src/cdg)


# Set base/common code sources. These add to 
# the list CDG_SRC:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsrc_common.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ginc_common.cmake)

# Set test target soures. The sources
# with a 'main' are in CDG_TEST_MAIN, and
# the rest are in CDG_SRC:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsrc_tests.cmake)

# Set default test id, and find desired test
# MAIN from list. Use this to set EXE target name:
#set(GTEST_NUM "1")
#list(GET CDG_TEST_MAIN ${GTEST_NUM} GTEST_NAME)
#list(FIND CDG_TEST_MAIN ${GTEST_NUM} GTEST_NUM_CHK)
#if ( ${GTEST_NUM_CHK} STREQUAL "-1" )
#  message(FATAL_ERROR "Invalid id: ${GTEST_NUM}:${GTEST_NAME} \n ")
#endif()
#string(REPLACE ".cpp" "" EXE)

# Fix compiler options here (don't let CMake do this for us):
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gcomp_options.cmake)


# Define libs/paths to link to (used based on target);
# may be changed in glocale.cmake:
set(MPI_PATH 
#   /scratch/duane.rosenberg/lib/mpich-3.2_gcc
   )

# GPTL path; may be changed in glocale.cmake:
set(GPTL_PATH 
#   /home/duane.rosenberg/lib/gptl-v5.5_gcc
#   /home/Duane.Rosenberg/lib/gptl-v5.4.4_pgi_openmpi
    /home/Duane.Rosenberg/lib/gptl-v5.4.4_intel_impi_theia
   )

# BOOST path; may be changed in glocale.cmake:
set(BOOST_ROOT 
#   /scratch/duane.rosenberg/lib/boost_1_69_0_gcc
    /home/Duane.Rosenberg/lib/boost_1_66_0_intel
   )

set(CDG_COMMON_LIBS  
    ${GPTL_PATH}/lib/libgptl.a
    ${BOOST_ROOT}/lib/libboost_chrono.a
    ${BOOST_ROOT}/lib/libboost_mpi.a
    ${BOOST_ROOT}/lib/libboost_serialization.a
    ${BOOST_ROOT}/lib/libboost_wserialization.a
    ${BOOST_ROOT}/lib/libboost_system.a
   )

if (HAVE_PAPI)
  set(CDG_COMMON_LIBS ${CDG_COMMON_LIBS} 
     /apps/papi/5.4.0/lib/libpapi.a
     )

endif()

set(CDG_COMMON_INC ${CDG_COMMON_INC}
    ${GPTL_PATH}/include
    ${BOOST_ROOT}/include
   )

if (USE_MPI)
  set(CDG_COMMON_LIBS ${CDG_COMMON_LIBS}
      ${MPI_PATH}/lib/libmpi.a
     )
set(CDG_COMMON_INC ${CDG_COMMON_INC}
    ${MPI_PATH}/include
   )
endif()


# Add executables here. If we want to build only one, XXX, 
# then invoke 'cmake --target XXX'. Derive executable name
# from source name:
foreach(srcmain ${CDG_TEST_MAIN})
  string(REPLACE ".cpp" "" GTMP ${srcmain})
  get_filename_component(GEXE ${GTMP} NAME)
  message(STATUS "Process executable ${GEXE}..." )
  add_executable(${GEXE} ${CDG_SRC} ${srcmain})
  target_include_directories(${GEXE} PUBLIC ${CDG_INC} ${CDG_COMMON_INC})
  target_link_libraries(${GEXE} ${CDG_COMMON_LIBS} )

# Set install directory for TARGETS:
  install(TARGETS ${GEXE} RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)

  message(STATUS "...................................... ${GEXE} added")
endforeach(srcmain)

