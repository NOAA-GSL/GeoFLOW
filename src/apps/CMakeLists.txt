#
# === Final Applications to Build ===
#

# 
# List of final applications to build
#
set(apps_files  
        geoflow.cpp
   )

#		       	
# Include headers unique to main applications
# Note:
# - Only used if "main.cpp" has a "main.hpp" that needs including
#
set(apps_include_headers "")

#
# Compiler options for each application
#
set(apps_compiler_options 
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
		"$<$<BOOL:${OpenMP_FOUND}>:${OpenMP_CXX_FLAGS}>"
		"$<$<BOOL:${MPI_CXX_FOUND}>:${MPI_CXX_FLAGS}>" 
	)

#
# Set the library names needed by applications
# Note:
# - These are only things not previously linked by a GeoFLOW library
#
set(apps_linked_libraries tbox)	
				             
#
# Loop over each building & installing
#	             
FOREACH(app_src ${apps_files})
	string(REPLACE ".cpp" "" app_name ${app_src})
	add_executable(${app_name} ${app_src})
	target_compile_options(${app_name} PRIVATE ${apps_compiler_options})
	target_include_directories(${app_name} PRIVATE ${apps_include_headers})
	target_link_libraries(${app_name} ${apps_linked_libraries})
	install(TARGETS ${app_name} DESTINATION "${MY_INSTALL_BINDIR}")
ENDFOREACH(app_src)

