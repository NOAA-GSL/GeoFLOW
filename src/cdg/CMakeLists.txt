###############################################
# Build the CDG Library
###############################################

add_library(cdg STATIC "")

#
# Get Name of Current Directory
#
set(GHOME          ${CMAKE_CURRENT_SOURCE_DIR})

#
# Define list of directories to include for "cdg"
#
target_include_directories(cdg
	PUBLIC  
	      ${GHOME}/include
          ${GHOME}/blas
          ${GHOME}/comm
          ${GHOME}/grid
          ${GHOME}/config
          ${GHOME}/io
          ${GHOME}/pdes
          ${GHOME}/sem
          ${GHOME}/solvers
          ${GHOME}/user/bdy
          ${GHOME}/user/force
          ${GHOME}/user/state
          ${GHOME}/user/terrain
          ${GHOME}/utils
	)
	
	
#
# Source files to compile into library
#
target_sources(cdg
	PRIVATE
	${GHOME}/blas/cff_blas.F
        ${GHOME}/blas/qmxmp.F
        ${GHOME}/blas/dmxmp.F
        ${GHOME}/blas/fmxmp.F
        ${GHOME}/blas/gcblas.cpp
        ${GHOME}/blas/gmtk.cpp
        
        ${GHOME}/comm/gcomm.cpp 
        
        ${GHOME}/sem/gelem_base.cpp
#       ${GHOME}/sem/ghelmholtz.cpp
        
        
        ${GHOME}/utils/gbitblock.cpp 
        ${GHOME}/utils/gutils.cpp
)

#
# Compile line options
# 
target_compile_options(cdg 
	PRIVATE
		"$<$<BOOL:${OpenACC_CXX_FOUND}>:${OpenACC_CXX_OPTIONS}>"
#		$<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>
	)

#
# Include directories for building library
# Note: 
# - The configure.hpp file is within ${MY_INSTALL_SRCDIR}
#   until it gets installed to ${MY_INSTALL_INCDIR} at the end.
#
target_include_directories(cdg
	PUBLIC
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
	)

#
# Links required by library
# 
target_link_libraries(cdg 
	PUBLIC
		tbox
	PRIVATE
		"$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_CXX>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
		"$<$<BOOL:${Boost_SERIALIZATION_FOUND}>:Boost::serialization>"
		"$<$<BOOL:${Boost_MPI_FOUND}>:Boost::mpi>"
		"$<$<BOOL:${PAPI_FOUND}>:PAPI::PAPI>"
		"$<$<BOOL:${GPTL_FOUND}>:GPTL::GPTL>"
		"$<$<BOOL:${OpenACC_CXX_FOUND}>:OpenACC::OpenACC_CXX>"	
		"$<$<TARGET_EXISTS:CUDA::cudart>:CUDA::cudart>"
                "$<$<TARGET_EXISTS:CUDA::cublas>:CUDA::cublas>"	
                "$<$<BOOL:${CBLAS_FOUND}>:CBLAS::CBLAS>"
	)

#
# Install public headers
#
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/cdg"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING
		PATTERN "*.hpp *.h"
)

#
# Install Library
#
install(TARGETS cdg
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
