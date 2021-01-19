#
# Create the "configure.hpp" file
# 
# src/configure.hpp.in --> build/src/configure.hpp
#
configure_file(
	"${MY_PROJECT_INCDIR}/configure.hpp.in"
	"${MY_INSTALL_SRCDIR}/configure.hpp"
)

#
# Make a header only library target
#
add_library(${PROJECT_NAME} INTERFACE)

#
# List the include directories for tests and library users
#
target_include_directories(${PROJECT_NAME} 
	INTERFACE
	    "$<BUILD_INTERFACE:${MY_PROJECT_INCDIR}>"
    	"$<INSTALL_INTERFACE:${MY_INSTALL_INCDIR}>"
)

#
# List of C++ features that are needed
#
target_compile_features(${PROJECT_NAME} 
	INTERFACE 
		cxx_std_20
)

install(
	TARGETS ${PROJECT_NAME}
	EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

#
# Copy all header files to the installation folder
#
install(DIRECTORY ${MY_PROJECT_INCDIR} DESTINATION ${MY_INSTALL_PREFIX})

#
# build/src/configure.hpp --> include/configure.hpp       
#
INSTALL(FILES "${MY_INSTALL_SRCDIR}/configure.hpp"
        DESTINATION "${MY_INSTALL_INCDIR}")     


add_subdirectory(tbox)
add_subdirectory(cdg)
add_subdirectory(apps)