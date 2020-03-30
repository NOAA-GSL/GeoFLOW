#
# === SandBox CDG Data Copy
#

#
# Install input files from data/
#
install(
	DIRECTORY "${MY_PROJECT_TSTDIR}/sandbox/cdg/data/"
	DESTINATION "${MY_INSTALL_TSTDIR}/sandbox/cdg/data"  
	FILES_MATCHING 
	PATTERN "*"
)



