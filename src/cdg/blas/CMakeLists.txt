

add_library(blas STATIC "")

target_sources(blas
	PRIVATE
		cff_blas.F
		dmxmp.F
		fmxmp.F
		qmxmp.F	
)

#target_compile_features(blas
#	PUBLIC
#		cxx_alias_templates
#	PRIVATE
#		cxx_nullptr	
#		cxx_static_assert
#)

target_compile_options(blas 
	PRIVATE
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
)

# Note: The configure.hpp file is within ${MY_INSTALL_SRCDIR}
# until it gets installed to ${MY_INSTALL_INCDIR} at the end.
target_include_directories(blas
	PUBLIC  
		"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
		"$<BUILD_INTERFACE:${MY_INSTALL_SRCDIR}>"
)

target_link_libraries(blas 
	PRIVATE
		"${Boost_LIBRARIES}"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_LIBRARIES}>"
)

# Install public headers
install(
	DIRECTORY "${MY_PROJECT_INCDIR}/blas"
	DESTINATION "${MY_INSTALL_INCDIR}"  
	FILES_MATCHING 
		PATTERN "*.hpp" 
		PATTERN "*.ipp"
)

# Install Library
install(TARGETS blas
		LIBRARY DESTINATION "${MY_INSTALL_LIBDIR}"
		ARCHIVE DESTINATION "${MY_INSTALL_LIBDIR}")
		
		
		
		
		