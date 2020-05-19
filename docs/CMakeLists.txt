
# Only add the docs target is DOxygen was previously found
if(DOXYGEN_FOUND)

# Set configuration input and output files
set(DOXYGEN_IN ${MY_PROJECT_DOCDIR}/Doxyfile.in)
set(DOXYGEN_OUT ${MY_INSTALL_DOCDIR}/Doxyfile)
    
# Process the configure file into build directory
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
   
# Create the "docs" target and provide commands to execute
# - Copy Markdown Pages into Build Directory
# - Run DOxygen on source files
set(DOXYGEN_EXTRACT_ALL YES) 
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)  
add_custom_target(docs
		COMMAND cp -r ${MY_PROJECT_DOCDIR}/doxy ${MY_INSTALL_DOCDIR}/
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${MY_INSTALL_DOCDIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )     

# Install the Markdown Pages
# - Requires "make install" which then compiles code (yuck!!!)
#install(
#	DIRECTORY "${MY_PROJECT_DOCDIR}/doxy"
#	DESTINATION "${MY_INSTALL_DOCDIR}/"  
#	FILES_MATCHING 
#	PATTERN "*"
#)

else(DOXYGEN_FOUND)
	message(VERBOSE "Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)