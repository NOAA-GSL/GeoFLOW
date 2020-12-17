#
# === Hack-a-Thon Applications to Build ===
#

#
# =========================================
#            Data Installation
# =========================================
#
install(
	DIRECTORY "${MY_PROJECT_SRCDIR}/apps/"
	DESTINATION "${MY_INSTALL_BINDIR}"  
	FILES_MATCHING 
	PATTERN "*.jsn"
)

#
# =========================================
#            Executables To Build
# =========================================
#

# 
# List of final applications to build
#
set(apps_files
       hackathon_driver.cpp
	)
	
#		       	
# Include headers unique to main applications
# Note:
# - Only used if "main.cpp" has a "main.hpp" that needs including
#
set(apps_include_headers 		
	PUBLIC
)

#
# Compiler options for each application
#
set(apps_compiler_options 
		"$<$<CXX_COMPILER_ID:AppleClang>:-Wno-comment>"
)

#
# Set the library names needed by applications
# Note:
# - These are only things not previously linked by a GeoFLOW library
#
set(apps_linked_libraries 
	PUBLIC
		cdg
)		
				             
#
# Loop over each building & installing
#	             
FOREACH(app_src ${apps_files})
	string(REPLACE ".cpp" "" app_name ${app_src})
	add_executable(${app_name} ${app_src})
	target_compile_options(${app_name} 
		PRIVATE 
				${apps_compiler_options}
#				 $<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>
	)
	target_include_directories(${app_name} PRIVATE ${apps_include_headers})
	target_link_libraries(${app_name} ${apps_linked_libraries})
	install(TARGETS ${app_name} DESTINATION "${MY_INSTALL_BINDIR}")
ENDFOREACH(app_src)

