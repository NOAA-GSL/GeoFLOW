#
# === Unit Tests
#

# 
# Build list of all tests to create
#
set(test_cdg_files
		  gtest_blas.cpp    
		  gtest_ggfx.cpp
		  gtest_gmtk.cpp
		# gtest_mass.cpp
		# gtest_advect.cpp  # Only compiles with GPTL
		# gtest_binding.cpp
		# gtest_burgers.cpp
		# gtest_class_acc.cpp
		# gtest_contig.cpp
		# gtest_derivs.cpp
		# gtest_ggrid.cpp
		# gtest_gll.cpp
		# gtest_gmm.cpp
		# gtest_helm.cpp
		# gtest_morton.cpp
		# gtest_pingpong.cpp
		# gtest_rk.cpp
		# gtest_vec_access.cpp
		# gtest_vec.cpp
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
	target_compile_options(${test_name} PRIVATE ${compiler_options})
	target_include_directories(${test_name} PRIVATE ${include_headers})
	target_link_libraries(${test_name} ${linked_libraries})
	add_test(NAME ${test_name} 
		     COMMAND ${LAUNCH_COMMAND} ${LAUNCH_ARGS} ${LAUNCH_PROCS} $<TARGET_FILE:${test_name}>)
ENDFOREACH(test_src)

