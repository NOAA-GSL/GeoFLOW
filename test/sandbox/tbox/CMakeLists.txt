#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(test_tbox_files
)

#
# Build list of compiler options for tests
#
set(compiler_options 
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
	)
		       
# Set include headers unique to testing
set(include_headers "")

# Set libraries the test need
# - Gets library include paths from here
set(linked_libraries 
		tbox 
		"$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_CXX>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
		"$<$<BOOL:${Boost_SERIALIZATION_FOUND}>:Boost::serialization>"
		"$<$<BOOL:${Boost_MPI_FOUND}>:Boost::mpi>"
		"$<$<BOOL:${PAPI_FOUND}>:PAPI::PAPI>"
		"$<$<BOOL:${GPTL_FOUND}>:GPTL::GPTL>"
	)


function(compile_test test_name)
	set(test_target test_${test_name})  # Target becomes test_xxxx
	set(test_source ${test_name}.cpp)   # Source becomes xxxx.cpp
	add_executable(${test_target} ${test_source})
	target_compile_options(${test_target} PRIVATE ${compiler_options})
	target_include_directories(${test_target} PRIVATE ${include_headers})
	target_link_libraries(${test_target} ${linked_libraries})
	set_target_properties(${test_target}
		PROPERTIES
    		CXX_STANDARD 17
    		CXX_STANDARD_REQUIRED YES
    		CXX_EXTENSIONS NO
	)
	#install(TARGETS ${test_target} DESTINATION  "${MY_INSTALL_TSTDIR}/sandbox/tbox")

	add_test(${test_name} ${test_target})
endfunction()

# Install test data
install(
	DIRECTORY "${MY_PROJECT_TSTDIR}/sandbox/tbox/data/"
	DESTINATION "${MY_INSTALL_TSTDIR}/sandbox/tbox"  
	FILES_MATCHING 
	PATTERN "*"
)