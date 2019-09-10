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
# build/src/configure.hpp --> include/configure.hpp       
#
INSTALL(FILES "${MY_INSTALL_SRCDIR}/configure.hpp"
        DESTINATION "${MY_INSTALL_INCDIR}")     


add_subdirectory(tbox)
add_subdirectory(cdg)

add_subdirectory(apps)