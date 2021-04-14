#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(test_tbox_files  
        enumerate.cpp
        indices.cpp
        input_manager.cpp
        morton.cpp
		mpixx.cpp
        property_tree.cpp
        range.cpp
        tbox_main.cpp
        zip.cpp
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
	)
		      
# 
# MPI launch commands
# CTest does not like them to be all within one command
# 
if(${GEOFLOW_USE_MPI})
	set(LAUNCH_COMMAND  "mpirun")  
	set(LAUNCH_ARGS     "-np")
	set(LAUNCH_PROCS    "4")
	set(LAUNCH_PREFLAGS "--oversubscribe")
else()
	set(LAUNCH_COMMAND  "")  
	set(LAUNCH_ARGS     "")
	set(LAUNCH_PROCS    "")
	set(LAUNCH_PREFLAGS "")
endif()

#
# Loop over each test building and adding to ctest
# - Build test name
#	             
FOREACH(test_src ${test_tbox_files})
	string(REPLACE ".cpp" "" test_name ${test_src})
	add_executable(${test_name} ${test_src})
	target_compile_options(${test_name} PRIVATE ${compiler_options})
	target_include_directories(${test_name} PRIVATE ${include_headers})
	target_link_libraries(${test_name} ${linked_libraries})
	
	add_test(NAME ${test_name} 
		     COMMAND ${LAUNCH_COMMAND} ${LAUNCH_ARGS} ${LAUNCH_PROCS} ${LAUNCH_PREFLAGS} $<TARGET_FILE:${test_name}>)
ENDFOREACH(test_src)


# Install test data
install(
	DIRECTORY "${MY_PROJECT_TSTDIR}/actions/tbox/data/"
	DESTINATION "${MY_INSTALL_TSTDIR}/actions/tbox"  
	FILES_MATCHING 
	PATTERN "*"
)
