###############################################
# Build the CDG Library
###############################################

add_library(cdg STATIC "")

#
# Load in the list of header files
#
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ginc_common.cmake)

#
# Load in the list of source files
#
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/gsrc_common.cmake)

#
# Source files to compile into library
#
target_sources(cdg
	PRIVATE
		${CDG_SRC}
)

#
# Request Fortran Linker
#
#set_property(TARGET cdg PROPERTY LINKER_LANGUAGE Fortran)

#
# Compiler standard required
#
option(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set_property(TARGET tbox PROPERTY CXX_STANDARD 11)

#
# Specific Compiler features required for library
# 
target_compile_features(cdg
	PUBLIC
		cxx_constexpr
		cxx_static_assert
)

#
# Compile line options
# 
target_compile_options(cdg 
	PRIVATE
#		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
#		"$<$<COMPILE_LANGUAGE:Fortran>:-lquadmath>"
	)

#
# Include directories for building library
# Note: 
# - The configure.hpp file is within ${MY_INSTALL_SRCDIR}
#   until it gets installed to ${MY_INSTALL_INCDIR} at the end.
#
target_include_directories(cdg
	PUBLIC  
		"${CDG_INC}"
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_INCLUDE_PATH}>"
		"$<$<BOOL:${USE_PAPI}>:${PAPI_INC_PATH}>"
		"$<$<BOOL:${USE_GPTL}>:${GPTL_INC_PATH}>"
	)

#
# Links required by library
# 
target_link_libraries(cdg 
	PUBLIC
		tbox
		"$<$<BOOL:${USE_PAPI}>:${PAPI_LIB_PATH}>"
		"$<$<BOOL:${USE_GPTL}>:${GPTL_LIB_PATH}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_LIBRARIES}>"
		"$<$<BOOL:${MPI_Fortran_FOUND}>:${MPI_Fortran_LIBRARIES}>"
	)

#
# Install public headers
#
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/cdg"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING
		PATTERN "*.hpp"
)

#
# Install Library
#
install(TARGETS cdg
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
