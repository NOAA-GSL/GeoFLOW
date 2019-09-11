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
set_property(TARGET tbox PROPERTY CXX_STANDARD 11)

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
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
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
		"${Boost_INCLUDE_DIRS}"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_INCLUDE_PATH}>"
)

#
# Links required by library
# 
target_link_libraries(tbox 
	PRIVATE
		"${Boost_LIBRARIES}"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_LIBRARIES}>"
)

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
		
		
		
		
		