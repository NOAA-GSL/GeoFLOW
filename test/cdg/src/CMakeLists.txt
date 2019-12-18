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
		  gtest_mass.cpp
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
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
	)
		       
# Set include headers unique to testing
set(include_headers "")

# Set libraries the test need
# - Gets library include paths from here
set(linked_libraries 
		cdg 
	)
		             
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
	add_test(${test_name} ${test_name})
ENDFOREACH(test_src)

