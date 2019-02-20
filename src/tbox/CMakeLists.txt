

add_library(tbox STATIC "")

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

target_compile_features(tbox
	PUBLIC
		cxx_alias_templates
	PRIVATE
		cxx_nullptr	
		cxx_static_assert
)

target_compile_options(tbox 
	PRIVATE
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>"	
)

# Note: The configure.hpp file is within ${MY_INSTALL_SRCDIR}
# until it gets installed to ${MY_INSTALL_INCDIR} at the end.
target_include_directories(tbox
	PUBLIC  
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
		"${Boost_INCLUDE_DIRS}"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_INCLUDE_PATH}>"
)

target_link_libraries(tbox 
	PRIVATE
		"${Boost_LIBRARIES}"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_LIBRARIES}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_LIBRARIES}>"
)

# Install public headers
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/tbox"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING PATTERN "*.hpp"
)

# Install Library
install(TARGETS tbox
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
		