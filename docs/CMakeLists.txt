
# Only add the docs target is DOxygen was previously found
if(DOXYGEN_FOUND)

# Set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    
# Request to configure the file
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    
set(DOXYGEN_EXTRACT_ALL YES) 
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)  
add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )  

else(DOXYGEN_FOUND)
	message(VERBOSE "Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)