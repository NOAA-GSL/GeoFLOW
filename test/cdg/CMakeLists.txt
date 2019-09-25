#
# === Unit Tests
#
add_subdirectory(src)


# Install test data
install(
	DIRECTORY "${MY_PROJECT_TSTDIR}/cdg/data"
	DESTINATION "${MY_INSTALL_TSTDIR}/cdg"  
	FILES_MATCHING 
	PATTERN "*"
)
