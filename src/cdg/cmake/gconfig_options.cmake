##############################################################
# GeoFLOW build configuration options
##############################################################

# Build configuration options may
# be changed here:

set   (GDIM "2")                           # Dimensionality

option(DO_TESTS       "DO_TESTS"    OFF)   # Build testing targets
option(DO_GEOFLOW     "DO_GEOFLOW"   ON)   # Build GeoFLOW target

option(DO_DEBUG       "DO_DEBUG"     ON)   # Compile for DEBUG
option(DO_OPENMP      "DO_OPENMP"    ON)   # Compile for OpenMP
option(DO_OPENACC     "DO_OPENACC"  OFF)   # Compile for OpenACC
option(USE_GBLAS      "USE_GBLAS"    ON)   # Use GBlas instead of C 
option(USE_MPI        "USE_MPI"      ON)   # Build with MPI?
option(USE_GPTL       "USE_GPTL"     ON)   # Build with GPTL?
option(HAVE_PAPI      "HAVE_PAPI"   OFF)   # Is GPTL built with PAPI support?
option(DO_AUTO_PROF   "DO_AUTO_PROF" OFF)  # Do auto-profing with GPTL?


if(NOT (DO_TESTS OR DO_GEOFLOW))
  message( FATAL_ERROR "Must set at least one of DO_GEOFLOW or DO_TESTS to bulid" )
endif()

