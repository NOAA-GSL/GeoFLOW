#
# Try to find PAPI library
#
# Once done this will define
#  PAPI_FOUND 	      - If system found PAPI library
#  PAPI_INCLUDE_DIRS  - The PAPI include directories
#  PAPI_LIBRARIES     - The libraries needed to use PAPI
#  PAPI_FLAGS 	      - Compiler flags required for using PAPI
#

#
# Set possible hints where to find PAPI
#
set(PAPI_HINTS 
	       ${PAPI_ROOT}
	       $ENV{PAPI_ROOT}
)

#
# Set installation guesses
#
set(PAPI_PATHS
	/usr
	/usr/local
	/usr/local/opt
)

#
# Attempt to find the include directories
#
find_path(PAPI_INCLUDE_DIRS
	NAMES papi.h
	HINTS 
	      ${PAPI_HINTS}
	PATH_SUFFIXES
		include
	PATHS 
	      ${PAPI_PATHS}
	DOC 
	    "PAPI header file papi.h"
)
mark_as_advanced(PAPI_INCLUDE_DIRS)

#
# Attempt to find the PAPI Library
#
find_library(PAPI_LIBRARIES
	NAMES 
	      libpapi.so
	      libpapi.a
	      papi
	HINTS
		${PAPI_HINTS}
	PATH_SUFFIXES
		lib
	PATHS
		${PAPI_PATHS}
	DOC
		"PAPI library"
)
mark_as_advanced(PAPI_LIBRARIES)

# ========================================================================

#
# Custom Notification Messages
#
if(NOT PAPI_INCLUDE_DIRS)
    message(STATUS "Could NOT find 'papi.h', install PAPI or set PAPI_ROOT")
endif()
if(NOT PAPI_LIBRARIES)
    message(STATUS "Could NOT find a libpapi*, install PAPI or set PAPI_ROOT")
endif()

#
# Determines whether library was found
# - Set the PAPI_FOUND variable depending on checks it does
#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAPI DEFAULT_MSG PAPI_INCLUDE_DIRS PAPI_LIBRARIES)

# =========================================================================