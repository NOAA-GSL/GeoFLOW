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
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
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
		
