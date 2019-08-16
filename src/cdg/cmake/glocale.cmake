##############################################################
# Compilers, & paths for libs required to build GeoFLOW 
##############################################################

##############################################################
# Specify compilers:
##############################################################

# Following may be changed by invoking
#   > cmake -DCMAKE_C_COMPILER=pgcc 
# etc..  These 'sets' must be done before 'project'
# is specified:
set(CMAKE_C_COMPILER       mpicc)
set(CMAKE_CXX_COMPILER     mpicxx)
set(CMAKE_Fortran_COMPILER mpifort)

# using TAU:
#set(THOME /home/Duane.Rosenberg/lib/tau-2.27/x86_64)
#set(CMAKE_C_COMPILER       ${THOME}/bin/taucc -tau:makefile=${THOME}/lib/Makefile.tau-papi)
#set(CMAKE_CXX_COMPILER     ${THOME}/bin/taucxx -tau:makefile=${THOME}/lib/Makefile.tau-papi)
#set(CMAKE_Fortran_COMPILER ${THOME}/bin/tauf90 -tau:makefile=${THOME}/lib/Makefile.tau-papi)


##############################################################
# Specify paths:
##############################################################

# BOOST path:
set(GLOCALE_BOOST_ROOT 
#   /scratch/duane.rosenberg/lib/boost_1_69_0_gcc
#   /home/Duane.Rosenberg/lib/boost_1_69_0_intel_18.1
    /Users/duane.rosenberg/lib/boost_1_69_0_gcc
   )

# MPI path:
set(GLOCALE_MPI_PATH 
#   /scratch/duane.rosenberg/lib/mpich-3.2_gcc
   )

# GPTL path:
set(GLOCALE_GPTL_PATH 
#   /home/duane.rosenberg/lib/gptl-v5.5_gcc
#   /home/Duane.Rosenberg/lib/gptl-v5.4.4_pgi_openmpi
#   /home/Duane.Rosenberg/lib/gptl-v5.4.4_intel_impi_theia
    /Users/duane.rosenberg/lib/gptl-5.6.0_clang
   )

# PAPI path: 
set(GLOCALE_PAPI_PATH
    /apps/papi/5.4.0/lib/libpapi.a
   )

##############################################################
# Don't modify below this line....
##############################################################

set(BOOST_INC_PATH    "") 
set(GPTL_LIB_PATH     "") 
set(GPTL_PAPI_PATH    "") 
set(MPI_INC_PATH      "") 
set(MPI_LIB_PATH      "") 
if (NOT "${GLOCALE_BOOST_ROOT}" STREQUAL "")
  set(BOOST_INC_PATH ${GLOCALE_BOOST_ROOT}/include)
endif()

if (NOT "${GLOCALE_GPTL_PATH}" STREQUAL "")
  set(GPTL_LIB_PATH ${GLOCALE_GPTL_PATH}/include)
endif()

if (NOT "${GLOCALE_PAPI_PATH}" STREQUAL "")
  set(GPTL_PAPI_PATH ${GLOCALE_PAPI_PATH})
endif()

if (NOT "${GLOCALE_MPI_PATH}" STREQUAL "")
  set(MPI_LIB_PATH ${GLOCALE_MPI_PATH}/lib/libmpi.a)
  set(MPI_INC_PATH ${GLOCALE_INC_PATH}/include)
endif()

