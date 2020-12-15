######################################################
#   Build tbox Library
######################################################
#
# Name of library to build
#
add_library(tbox STATIC "")

#
# Source files to compile into library
#
target_sources(tbox
	PRIVATE
		clock.cpp
		command_line.cpp
		display_header.cpp
		error_handler.cpp
		global_manager.cpp
		input_manager.cpp
		io_buffer.cpp
		mpixx.cpp
		pio.cpp
		property_tree.cpp
		tracer.cpp	
)

#
# Compiler standard required
#
option(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set_property(TARGET tbox PROPERTY CXX_STANDARD 17)

#
# Specific Compiler features required for library
# 
target_compile_features(tbox
	PUBLIC
		cxx_constexpr
		cxx_static_assert
)

#
# Compile line options
# 
target_compile_options(tbox 
	PRIVATE
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
)

#
# Include directories for building library
# Note: 
# - The configure.hpp file is within ${MY_INSTALL_SRCDIR}
#   until it gets installed to ${MY_INSTALL_INCDIR} at the end.
#
target_include_directories(tbox
	PUBLIC  
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
)


#
# Links required by library
# 
target_link_libraries(tbox 
	PUBLIC
		"$<$<BOOL:${Boost_SERIALIZATION_FOUND}>:Boost::serialization>"
		"$<$<BOOL:${Boost_MPI_FOUND}>:Boost::mpi>"
		"$<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::nvToolsExt>"
		"$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_CXX>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
)

#target_link_options(tbox 
#	PUBLIC
#		"$<$<BOOL:${CUDA_FOUND}>:-lnvToolsExt>"
#)

#
# Install public headers
#
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/tbox"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING
		PATTERN "*.hpp"
)

#
# Install Library
#
install(TARGETS tbox
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
		
