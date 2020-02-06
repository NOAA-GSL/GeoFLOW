#
# Try to find PAPI library
#
#
#  -- Prefered Method --
#  PAPI_ROOT         Set this variable to the root installation of PAPI
#
#  -- Alternatively --
#  A user can specify the header and library directories independently
#
#  PAPI_INCLUDE_DIR  The directory containing papi.h
#
#  PAPI_LIBRARY      The directory containing the library
#
#
# Once done this will define
#  PAPI_FOUND 	      - If system found PAPI library
#  PAPI_INCLUDE_DIRS  - The PAPI include directories
#  PAPI_LIBRARIES     - The libraries needed to use PAPI
#

# Find_PAPI
# ----------------
#
# Finds the PAPI library.
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# This module provides the following imported targets, if found:
#
# "PAPI::PAPI"
# The PAPI library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This will define the following variables:
#
# "PAPI_FOUND"         - True if the system has the PAPI library.
# "PAPI_INCLUDE_DIRS"  - Include directories needed to use PAPI.
# "PAPI_LIBRARIES"     - Libraries needed to link to PAPI.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# The following cache variables may also be set:
#
# "PAPI_ROOT"        - The root directory of the PAPI library.
# "PAPI_INCLUDE_DIR" - The directory containing "papi.h".
# "PAPI_LIBRARY"     - The path to the PAPI library.
#

if(USE_PAPI)
message(VERBOSE "")
message(VERBOSE "--------------------- PAPI Libraries -----------------------")
message(VERBOSE "")
message(VERBOSE "Search Locations:")
message(VERBOSE "PAPI_ROOT            = ${PAPI_ROOT}")
message(VERBOSE "PAPI_INCLUDE_DIR     = ${PAPI_INCLUDE_DIR}")
message(VERBOSE "PAPI_INCLUDE_LIBRARY = ${PAPI_INCLUDE_LIBRARY}")

find_path(PAPI_INCLUDE_DIR
    NAMES papi.h
    PATHS ${PAPI_ROOT} $ENV{PAPI_ROOT}
    PATH_SUFFIXES 
    	include
)

find_library(PAPI_LIBRARY
    NAMES libpapi.so libpapi.a papi
    PATHS ${PAPI_ROOT} $ENV{PAPI_ROOT}
    PATH_SUFFIXES 
    	lib 
    	lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAPI 
	FOUND_VAR
		PAPI_FOUND
	REQUIRED_VARS
    	PAPI_LIBRARY
    	PAPI_INCLUDE_DIR
)

#
# Set the returned values
#
if(PAPI_FOUND)
  set(PAPI_LIBRARIES ${PAPI_LIBRARY})
  set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
endif()

#
# Allow these values to be in cache
#
mark_as_advanced(
    PAPI_ROOT
    PAPI_LIBRARY
    PAPI_INCLUDE_DIR
)

#
# Define Imported Targets
#
if(PAPI_FOUND AND NOT TARGET PAPI::PAPI)
  add_library(PAPI::PAPI UNKNOWN IMPORTED)
  set_target_properties(PAPI::PAPI PROPERTIES
    IMPORTED_LOCATION "${PAPI_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PAPI_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PAPI_INCLUDE_DIR}"
  )
endif()

#
# Display Results if Requested
#
message(VERBOSE "")
message(VERBOSE "Results:")
message(VERBOSE "PAPI Found     = ${PAPI_FOUND}")
message(VERBOSE "PAPI Includes  = ${PAPI_INCLUDE_DIRS}")
message(VERBOSE "PAPI Libraries = ${PAPI_LIBRARIES}")

endif(USE_PAPI)