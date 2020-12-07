#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(test_cdg_files 
		  cdg_blas.cpp 
		  cdg_cg.cpp    
		  cdg_ggfx.cpp
		  cdg_gmtk.cpp
		  cdg_mass.cpp
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
		cdg 
		"$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_CXX>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
		"$<$<BOOL:${Boost_SERIALIZATION_FOUND}>:Boost::serialization>"
		"$<$<BOOL:${Boost_MPI_FOUND}>:Boost::mpi>"
		"$<$<BOOL:${PAPI_FOUND}>:PAPI::PAPI>"
		"$<$<BOOL:${GPTL_FOUND}>:GPTL::GPTL>"
	)
	
# 
# MPI launch commands
# CTest does not like them to be all within one command
# 
if(${MPI_FOUND})
	set(LAUNCH_COMMAND "mpirun")  
	set(LAUNCH_ARGS    "-np")
	set(LAUNCH_PROCS   "4")
else()
	set(LAUNCH_COMMAND "")  
	set(LAUNCH_ARGS    "")
	set(LAUNCH_PROCS   "")
endif()	       
      
#
# Loop over each test building and adding to ctest
# - Build test name
#	             
FOREACH(test_src ${test_cdg_files})
	string(REPLACE ".cpp" "" test_name ${test_src})
	add_executable(${test_name} ${test_src})
	target_compile_options(${test_name} 
		PRIVATE 
			${compiler_options} 
#			$<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>
	)
	option(CMAKE_CXX_STANDARD_REQUIRED TRUE)
	set_property(TARGET ${test_name} PROPERTY CXX_STANDARD 17)
	target_include_directories(${test_name} PRIVATE ${include_headers})
	target_link_libraries(${test_name} ${linked_libraries})
	add_test(NAME ${test_name} 
		     COMMAND ${LAUNCH_COMMAND} ${LAUNCH_ARGS} ${LAUNCH_PROCS} $<TARGET_FILE:${test_name}>)
ENDFOREACH(test_src)

#
# Copy Test data into executable directory
#

# Install test data
install(
	DIRECTORY "${MY_PROJECT_TSTDIR}/travis/cdg/data/"
	DESTINATION "${MY_INSTALL_TSTDIR}/travis/cdg"  
	FILES_MATCHING 
	PATTERN "*"
)


