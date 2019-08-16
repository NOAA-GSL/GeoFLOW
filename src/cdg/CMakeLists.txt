#==================================================================================
# Module       : CMakeLists.txt
# Date         : 8/16/19 (DLR)
# Description  : GeoFIRE executable cmake file
# Copyright    : Copyright 2019. Colorado State University. All rights reserved.
#==================================================================================
#
cmake_minimum_required(VERSION 2.8)

# Note: call locale first, since compiler specification
#       must be done bfore setting project. This sets 
#       the compilers, and lib paths. Sets vars:
#       BOOST_INC_PATH, GPTL_LIB(INCLUDE)_PATH, GPTL_PAPI_PATH,
#       MPI_INC_PATH, MPI_LIB_PATH:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/glocale.cmake)


project(GeoFLOW_CDG)

# The version number.
set(GTest_VERSION_MAJOR 1)
#set(GTest_VERSION_MINOR 0)
message("CMAKE_VERSION: ${CMAKE_VERSION}")


set(GEOFLOW_SRC "geoflow.cpp")
set(GEOFLOW_EXE "geoflow")

enable_testing()

enable_language (Fortran)

########################################################################
# Load in project macros/functions:
########################################################################
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gmacros.cmake)

########################################################################
# Define and set default CMake GeoFLOW build options:
########################################################################
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gconfig_options.cmake)

########################################################################
# Miscellaneous CMake specifications:
########################################################################
# Following should be set in call to cmake as
# cmake -DCMAKE_BUILD_TYPE=DEBUG [or RELEASE]:
set(CMAKE_BUILD_TYPE DEBUG) 
if (NOT DO_DEBUG)
  set(CMAKE_BUILD_TYPE RELEASE) 
  message("CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
endif()

set (CMAKE_CXX_STANDARD 11) # good for only cmake v3+

########################################################################
# Add list of source files:
########################################################################
# Set base/common code sources. These add to 
# the list CDG_SRC:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsrc_common.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ginc_common.cmake)

# Set test target soures. The sources
# the rest are in tests CDG_TEST_SRC:
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gtests_src.cmake)


# Set compiler options here (don't auto-detect these):
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gcomp_options.cmake)

########################################################################
# Configure Boost environment:
########################################################################
set(Boost_DEBUG OFF)              # Enable debug output from FIND_PACKAGE
set(Boost_NO_SYSTEM_PATHS ON)     # Do not search system paths before BOOST_ROOT
set(Boost_USE_STATIC_LIBS ON)     # Static link to Boost libraries
set(Boost_USE_STATIC_RUNTIME OFF) # Use Boost static linked to C++ runtime
set(Boost_USE_MULTITHREADED OFF)  # Use Boost multi-threaded code
#FIND_PACKAGE(Boost 1.65.0 COMPONENTS mpi serialization REQUIRED)


########################################################################
# Set misc. variables....
########################################################################
set(CDG_COMMON_LIBS  
    ${GLOCALE_BOOST_ROOT}/lib/libboost_chrono.a
    ${GLOCALE_BOOST_ROOT}/lib/libboost_mpi.a
    ${GLOCALE_BOOST_ROOT}/lib/libboost_serialization.a
    ${GLOCALE_BOOST_ROOT}/lib/libboost_wserialization.a
    ${GLOCALE_BOOST_ROOT}/lib/libboost_system.a
   )
set(CDG_COMMON_INC ${CDG_COMMON_INC}
    ${BOOST_INC_PATH}
   )

if (USE_GPTL)
  set(CDG_COMMON_LIBS ${CDG_COMMON_LIBS} 
    ${GPTL_LIB_PATH}
     )
  set(CDG_COMMON_INC ${CDG_COMMON_INC}
      ${GPTL_INC_PATH}
     )
  if (HAVE_PAPI)
    set(CDG_COMMON_LIBS ${CDG_COMMON_LIBS} 
       ${GPTL_PAPI_PATH}
       )
  endif()
endif()

if (USE_MPI)
  set(CDG_COMMON_LIBS ${CDG_COMMON_LIBS}
      ${MPI_LIB_PATH}
     )
  set(CDG_COMMON_INC ${CDG_COMMON_INC}
      ${MPI_INC_PATH}
     )
endif()


########################################################################
# Add build targets:
########################################################################
# First, add the main executable, if required to:
# if required to:
if (DO_GEOFLOW)
  message(STATUS "Process executable ${GEOFLOW_EXE}..." )
  add_executable(${GEXE} ${CDG_SRC} ${GEOFLOW_SRC})
  target_include_directories(${GEXE} PUBLIC ${CDG_INC} ${CDG_COMMON_INC})
  target_link_libraries(${GEOFLOW_EXE} ${CDG_COMMON_LIBS} )

  # Set install directory for main target:
  install(TARGETS ${GEXE} RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
endif()

# Add the tests executable, if required to:
if (DO_TESTS)
  foreach(cdgtests ${CDG_TEST_SRC})
    string(REPLACE ".cpp" "" GTMP ${cdgtests})
    get_filename_component(GEXE ${GTMP} NAME)
    message(STATUS "Process executable ${GEXE}..." )
    add_executable(${GEXE} ${CDG_SRC} ${cdgtests})
    target_include_directories(${GEXE} PUBLIC ${CDG_INC} ${CDG_COMMON_INC})
    target_link_libraries(${GEXE} ${CDG_COMMON_LIBS} )
    add_test(${GEXE} ${GEXE})
  # Set install directory for TARGETS:
    install(TARGETS ${GEXE} RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
  
    message(STATUS "...................................... ${GEXE} added")
  endforeach(cdgtests)
endif()

