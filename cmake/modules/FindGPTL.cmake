#
# Find_GPTL
# ----------------
#
# Finds the GPTL library.
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# This module provides the following imported targets, if found:
#
# ""GPTL::GPTL"
# The GPTL library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This will define the following variables:
#
# "GPTL_FOUND"         - True if the system has the GPTL library.
# "GPTL_INCLUDE_DIRS"  - Include directories needed to use GPTL.
# "GPTL_LIBRARIES"     - Libraries needed to link to GPTL.
#
# Cache Variables
# ^^^^^^^^^^^^^^^
#
# The following cache variables may also be set:
#
# "GPTL_ROOT"        - The root directory of the GPTL library.
# "GPTL_INCLUDE_DIR" - The directory containing "gptl.h".
# "GPTL_LIBRARY"     - The path to the GPTL library.
#

message(VERBOSE "")
message(VERBOSE "--------------------- GPTL Libraries -----------------------")
message(VERBOSE "")
message(VERBOSE "Search Locations:")
message(VERBOSE "GPTL_ROOT            = ${GPTL_ROOT}")
message(VERBOSE "GPTL_INCLUDE_DIR     = ${GPTL_INCLUDE_DIR}")
message(VERBOSE "GPTL_INCLUDE_LIBRARY = ${GPTL_INCLUDE_LIBRARY}")

find_path(GPTL_INCLUDE_DIR
    NAMES gptl.h
    PATHS ${GPTL_ROOT} $ENV{GPTL_ROOT}
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
)

find_library(GPTL_LIBRARY
    NAMES libgptl.so libgptl.a libgptl libgptl_pmpi.so libgptl_pmpi.a libgptl_pmpi
    PATHS ${GPTL_ROOT} $ENV{GPTL_ROOT}
    PATH_SUFFIXES 
    	lib 
    	bluegene/lib
		cray/lib
		lahey/lib
		linux/lib
		macos/lib
		pgi/lib
		xeon/lib
		xeonphi/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPTL 
	FOUND_VAR
		GPTL_FOUND
	REQUIRED_VARS
    	GPTL_LIBRARY
    	GPTL_INCLUDE_DIR
    FAIL_MESSAGE
    	"GPTL was not found. Set GPTL_ROOT to installation directory" 	
)

#
# Set the returned values
#
if(GPTL_FOUND)
  set(GPTL_LIBRARIES ${GPTL_LIBRARY})
  set(GPTL_INCLUDE_DIRS ${GPTL_INCLUDE_DIR})
endif()

#
# Allow these values to be in cache
#
mark_as_advanced(
    GPTL_ROOT
    GPTL_LIBRARY
    GPTL_INCLUDE_DIR
)

#
# Define Imported Targets
#
if(GPTL_FOUND AND NOT TARGET GPTL::GPTL)
  add_library(GPTL::GPTL UNKNOWN IMPORTED)
  set_target_properties(GPTL::GPTL PROPERTIES
    IMPORTED_LOCATION "${GPTL_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_GPTL_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${GPTL_INCLUDE_DIR}"
  )
endif()

#
# Display Results if Requested
#
message(VERBOSE "")
message(VERBOSE "Results:")
message(VERBOSE "GPTL Found     = ${GPTL_FOUND}")
message(VERBOSE "GPTL Includes  = ${GPTL_INCLUDE_DIRS}")
message(VERBOSE "GPTL Libraries = ${GPTL_LIBRARIES}")

