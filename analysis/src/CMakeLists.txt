#
# === Analysis tools to Build ===
#

# 
# List of analysis tools to build
#
set(analysis_files
        gcdg_pdf1d.cpp
)

#		       	
# Include headers unique to analysis tools
# Note:
# - Only used if "main.cpp" has a "main.hpp" that needs including
#
 set(analysis_include_headers "")

#
# Compiler options for each tool
#
set(analysis_compiler_options 
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
)

#
# Set the library names needed by tools
# Note:
# - These are only things not previously linked by a GeoFLOW library
#
set(analysis_linked_libraries 
	PRIVATE
		tbox 
		cdg
		"$<$<BOOL:${OpenMP_FOUND}>:OpenMP::OpenMP_CXX>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
		"$<$<BOOL:${Boost_SERIALIZATION_FOUND}>:Boost::serialization>"
		"$<$<BOOL:${Boost_MPI_FOUND}>:Boost::mpi>"
		"$<$<BOOL:${PAPI_FOUND}>:PAPI::PAPI>"
		"$<$<BOOL:${GPTL_FOUND}>:GPTL::GPTL>"
)	
				             
#
# Loop over each building & installing
#	             
FOREACH(analysis_src ${analysis_files})
	string(REPLACE ".cpp" "" analysis_name ${analysis_src})
	add_executable(${analysis_name} ${analysis_src})
	target_compile_options(${analysis_name} PRIVATE ${analysis_compiler_options})
	target_include_directories(${analysis_name} PRIVATE ${analysis_include_headers})
	target_link_libraries(${analysis_name} ${analysis_linked_libraries})
	install(TARGETS ${analysis_name} DESTINATION "${MY_INSTALL_BINDIR}")
ENDFOREACH(analysis_src)

