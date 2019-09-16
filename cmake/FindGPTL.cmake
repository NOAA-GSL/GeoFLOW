#
# Try to find GPTL library
#
# Once done this will define
#  GPTL_FOUND 	      - If system found GPLT library
#  GPTL_INCLUDE_DIRS  - The GPTL include directories
#  GPTL_LIBRARIES     - The libraries needed to use GPTL
#  GPTL_FLAGS 	      - Compiler flags required for using GPTL
#

#
# Set possible hints where to find GPTL
#
set(GPTL_HINTS 
	       ${GPTL_ROOT}
	       $ENV{GPLT_ROOT}
)

#
# Set installation guesses
#
set(GPTL_PATHS
	/usr
	/usr/local
	/usr/local/opt
)

#
# Attempt to find the include directories
#
find_path(GPTL_INCLUDE_DIRS
	NAMES gptl.h
	HINTS 
	      ${GPTL_HINTS}
	PATH_SUFFIXES
		include 
		bluegene/include
		cray/include
		lahey/include
		linux/include
		macos/include
		pgi/include	
		xeon/include
		xeonphi/include
	PATHS 
	      ${GPTL_PATHS}
	DOC 
	    "GPTL header file gptl.h"
)
mark_as_advanced(GPTL_INCLUDE_DIRS)

#
# Attempt to find the GPTL Library
#
find_library(GPTL_LIBRARIES
	NAMES 
	      gptl
	      gptl_pmpi
	HINTS
		${GPTL_HINTS}
	PATH_SUFFIXES
		lib
		bluegene/include
		cray/include
                lahey/include
                linux/include
                macos/include
                pgi/include
		xeon/include
                xeonphi/include
	PATHS
		${GPTL_PATHS}
	DOC
		"GPTL library"
)
mark_as_advanced(GPTL_LIBRARIES)

# ========================================================================

# Notification messages
if(NOT GPTL_INCLUDE_DIRS)
    message(STATUS "Could NOT find 'gptl.h', install a GPTL or set GPTL_ROOT")
endif()
if(NOT GPTL_LIBRARIES)
    message(STATUS "Could NOT find a libgptl..., install it or set GPTL_ROOT")
endif()

# Determines whether or not GPTL was found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPTL DEFAULT_MSG GPTL_INCLUDE_DIRS GPTL_LIBRARIES)

# =========================================================================