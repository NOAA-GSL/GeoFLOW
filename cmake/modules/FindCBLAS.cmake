
# Find CBLAS Headers & Library 
#
# This wrapper attempts to find the CBLAS
# library for Basic Linear Algebra Subroutines.
#
# To Use:
# =====================================================================
#
# Use this module by invoking find_package within CMake: 
#  find_package(CBLAS
#			[REQUIRED]  # Stop configuration if not found
#               )
#
# Dependancies:
# =====================================================================
# BLAS - Library must be found by this script
#
# Variables To Locate:
# =====================================================================
#
# - BLAS:
# BLAS_ROOT       = Root directory to search for BLAS
# BLAS_DIR        = Where to find the base directory of BLAS
# BLAS_INCDIR     = Where to find the header files
# BLAS_LIBDIR     = Where to find the library files
# BLAS_VERBOSE    = Print additional information 
#
# - CBLAS:
# CBLAS_ROOT      = Root directory to search for CBLAS
# CBLAS_DIR       = Where to find the base directory of CBLAS
# CBLAS_INCDIR    = Where to find the header files
# CBLAS_LIBDIR    = Where to find the library files
# CBLAS_VERBOSE    = Print additional information 
#
#
# Result Variables
# =====================================================================
#
# BLAS Variables:
# BLAS_FOUND          = True if Library found
# BLAS_LINKER_FLAGS   = Uncached list of required linker flags (excluding -l and -L)
# BLAS_COMPILER_FLAGS = Uncached list of required compiler flags (including -I for mkl headers).
# BLAS_LIBRARIES      = Uncached list of libraries (using full path name) to link against to use BLAS
# BLAS_VENDOR_FOUND   = Name of BLAS vendor 
#
#
#
#  CBLAS_LIBRARIES_DEP

# ---------------------------------------------------------------------
#                               Find BLAS
# ---------------------------------------------------------------------

# Attempt to Find BLAS (a Dependancy)
if (NOT BLAS_FOUND)
	set(BLAS_DIR ${BLAS_ROOT} CACHE PATH "Installation directory of BLAS library")
	if(CBLAS_FIND_REQUIRED)
		find_package(BLAS REQUIRED)
	else()
	    find_package(BLAS)
	endif()
endif() 

if(NOT CBLAS_FIND_QUIETLY)
	message(STATUS "BLAS_FOUND:             ${BLAS_FOUND}")
	message(STATUS "BLAS_LINKER_FLAGS:      ${BLAS_LINKER_FLAGS}")
	message(STATUS "BLAS_COMPILER_FLAGS:    ${BLAS_COMPILER_FLAGS}")
	message(STATUS "BLAS_COMPILER_FLAGS:    ${BLAS_COMPILER_FLAGS}")
	message(STATUS "BLAS_LIBRARIES:         ${BLAS_LIBRARIES}")
	message(STATUS "BLAS_VENDOR_FOUND:      ${BLAS_VENDOR_FOUND}")
endif() 

# ---------------------------------------------------------------------
#                               Find CBLAS
# ---------------------------------------------------------------------

if(BLAS_FOUND)
	include(CheckFunctionExists)
	
	# Check if CBLAS Functions exist in BLAS Library
	set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LINKER_FLAGS};${BLAS_LIBRARIES}")     
	set(CMAKE_REQUIRED_FLAGS "${BLAS_COMPILER_FLAGS}")
	unset(CBLAS_WORKS CACHE)
	check_function_exists(cblas_dscal CBLAS_WORKS)
	mark_as_advanced(CBLAS_WORKS)

	# If CBLAS inside BLAS then set variables
	if(CBLAS_WORKS)
		if(NOT CBLAS_FIND_QUIETLY)
			message(STATUS "CBLAS found inside BLAS")
		endif()
		
		# Assign variables for BLAS to CBLAS
		set(CBLAS_LIBRARIES "${BLAS_LIBRARIES}")
		set(CBLAS_LIBRARIES_DEP "${BLAS_LIBRARIES}")
		if (BLAS_LIBRARY_DIRS)
			set(CBLAS_LIBRARY_DIRS "${BLAS_LIBRARY_DIRS}")
		endif()       
		if(BLAS_INCLUDE_DIRS)
			set(CBLAS_INCLUDE_DIRS "${BLAS_INCLUDE_DIRS}")
			set(CBLAS_INCLUDE_DIRS_DEP "${BLAS_INCLUDE_DIRS_DEP}")
		endif()
		if (BLAS_LINKER_FLAGS)
			set(CBLAS_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
		endif()  
	endif(CBLAS_WORKS)

  	if (NOT CBLAS_WORKS)
  		if(NOT CBLAS_FIND_QUIETLY)
			message(STATUS "CBLAS was NOT found inside BLAS will search external")
		endif()
		set(CBLAS_DIR ${CBLAS_ROOT} CACHE PATH "Installation directory of CBLAS library")
  	
  		# ---------- Searching For cblas.h ----------
  		
  		# Build a List of possible locations to search
  		unset(cblas_include_hints)
  		
		# - Get Environment Vals
		set(ENV_CBLAS_ROOT "$ENV{CBLAS_ROOT}")
		set(ENV_CBLAS_DIR "$ENV{CBLAS_DIR}")
		set(ENV_CBLAS_INCDIR "$ENV{CBLAS_INCDIR}")
		if(ENV_CBLAS_ROOT)
			set(ENV_BLAS_DIR ${ENV_CBLAS_ROOT} CACHE PATH "Installation directory of CBLAS library")
		endif()
		
		# If provided an Include
		if(ENV_CBLAS_INCDIR)
			list(APPEND cblas_include_hints "${ENV_CBLAS_INCDIR}")
			
		# If Provided an Root or Dir
		elseif(ENV_CBLAS_DIR)	
		    list(APPEND cblas_include_hints "${ENV_CBLAS_DIR}")       
		    list(APPEND cblas_include_hints "${ENV_CBLAS_DIR}/include")       
		    list(APPEND cblas_include_hints "${ENV_CBLAS_DIR}/include/cblas")
		    
		# If nothing provided then try system stuff
		else()
		    string(REPLACE ":" ";" sys_path_env "$ENV{INCLUDE}")         
		    list(APPEND cblas_include_hints "${sys_path_env}")
		    string(REPLACE ":" ";" sys_path_env "$ENV{C_INCLUDE_PATH}")
		    list(APPEND cblas_include_hints "${sys_path_env}")
		    string(REPLACE ":" ";" sys_path_env "$ENV{CPATH}")
		    list(APPEND cblas_include_hints "${sys_path_env}")
		    string(REPLACE ":" ";" sys_path_env "$ENV{INCLUDE_PATH}")
		    list(APPEND cblas_include_hints "${sys_path_env}")       
		endif()
		list(APPEND cblas_include_hints "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")     
		list(APPEND cblas_include_hints "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")     
		list(REMOVE_DUPLICATES cblas_include_hints)
		
		# Try to find the cblas.h header file
		if(CBLAS_INCDIR)
			set(CBLAS_HEADER_DIR "NOTFOUND")
			find_path(CBLAS_HEADER_DIR
				NAMES cblas.h
				HINTS ${CBLAS_INCDIR})
		else()
			if(CBLAS_DIR)
				set(CBLAS_HEADER_DIR "NOTFOUND")
				find_path(CBLAS_HEADER_DIR
				    NAMES cblas.h
				    HINTS ${CBLAS_DIR}
				    PATH_SUFFIXES "include" "include/cblas")
		    else()
  				set(CBLAS_HEADER_DIR "NOTFOUND")
  			    find_path(CBLAS_HEADER_DIR
  			    	NAMES cblas.h
  			    	HINTS ${cblas_include_hints}
  			    	PATH_SUFFIXES "cblas")
  			endif()
  		endif()
  		mark_as_advanced(CBLAS_HEADER_DIR)
  		
  		# If found then add header to variable
  		if(CBLAS_HEADER_DIR)
  			set(CBLAS_INCLUDE_DIRS "${CBLAS_HEADER_DIR}")
  			if(NOT CBLAS_FIND_QUIETLY)
  				message(STATUS "Found cblas.h at ${CBLAS_INCLUDE_DIRS}")
  			endif()
  		else()
  			set(CBLAS_INCLUDE_DIRS "CBLAS_INCLUDE_DIRS-NOTFOUND")
  			if(NOT CBLAS_FIND_QUIETLY)
  				message(STATUS "Not Found -> cblas.h")
  			endif()
  		endif(CBLAS_HEADER_DIR)
  		
  		# ---------- Searching For Library ----------
  		unset(cblas_library_hints)
  		set(ENV_CBLAS_LIBDIR "$ENV{CBLAS_LIBDIR}")
  		if(ENV_CBLAS_LIBDIR)
  			list(APPEND cblas_library_hints "${ENV_CBLAS_LIBDIR}")
  		elseif(ENV_CBLAS_DIR)       
  		   	list(APPEND cblas_library_hints "${ENV_CBLAS_DIR}")
  		   	list(APPEND cblas_library_hints "${ENV_CBLAS_DIR}/lib")
  		else()
	       	if(APPLE)
	      		string(REPLACE ":" ";" cblas_library_hints "$ENV{DYLD_LIBRARY_PATH}")
	       	else()
	       		string(REPLACE ":" ";" cblas_library_hints "$ENV{LD_LIBRARY_PATH}")
	       	endif()
	       	list(APPEND cblas_library_hints "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
	       	list(APPEND cblas_library_hints "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
	    endif()
	    list(REMOVE_DUPLICATES cblas_library_hints)
  		
  		# Try to find the cblas. library file
  		if(CBLAS_LIBDIR)
  			set(CBLAS_LIBRARY_FULL_FILE "NOTFOUND")
  			find_library(CBLAS_LIBRARY_FULL_FILE
  				NAMES cblas
  				HINTS ${CBLAS_LIBDIR})
  		else()
  			if(CBLAS_DIR)
  				set(CBLAS_LIBRARY_FULL_FILE "NOTFOUND")
  				find_library(CBLAS_LIBRARY_FULL_FILE
  					NAMES cblas
  					HINTS ${CBLAS_DIR}
  					PATH_SUFFIXES lib lib32 lib64)
  			else()
  				set(CBLAS_LIBRARY_FULL_FILE "NOTFOUND")
  				find_library(CBLAS_LIBRARY_FULL_FILE
  					NAMES cblas
  					HINTS ${cblas_library_hints})
  			endif()
  		endif()
  		mark_as_advanced(CBLAS_LIBRARY_FULL_FILE)
  		
  		# If found then add library to variable
  		if(CBLAS_LIBRARY_FULL_FILE)
  			get_filename_component(CBLAS_LIBRARY_PATH "${CBLAS_LIBRARY_FULL_FILE}" PATH)
       		set(CBLAS_LIBRARIES    "${CBLAS_LIBRARY_FULL_FILE}")
       		set(CBLAS_LIBRARY_DIRS "${CBLAS_LIBRARY_PATH}")
       		if (NOT CBLAS_FIND_QUIETLY)
       			message(STATUS "lib cblas found at ${CBLAS_LIBRARIES}")
       		endif()
       	else ()
       		set(CBLAS_LIBRARIES    "CBLAS_LIBRARIES-NOTFOUND")
       		set(CBLAS_LIBRARY_DIRS "NOTFOUND")
       		if (NOT CBLAS_FIND_QUIETLY)
       			message(STATUS "lib cblas not found")
       		endif()
  		endif(CBLAS_LIBRARY_FULL_FILE)
  		
  		# ---------- Test Library ----------
  		if(CBLAS_LIBRARIES)
  			set(REQUIRED_INCDIRS)
  			set(REQUIRED_LDFLAGS)
  			set(REQUIRED_LIBDIRS)
  			set(REQUIRED_LIBS)
  		
  			# CBLAS
  			if(CBLAS_INCLUDE_DIRS)
  				set(REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
  			endif()
  			if(CBLAS_LIBRARY_DIRS)
  				set(REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
  			endif()
  			set(REQUIRED_LIBS "${CBLAS_LIBRARIES}")
  			# BLAS
  			if(BLAS_INCLUDE_DIRS)
  				list(APPEND REQUIRED_INCDIRS "${BLAS_INCLUDE_DIRS}")
  			endif()
  			if(BLAS_LIBRARY_DIRS)
  				list(APPEND REQUIRED_LIBDIRS "${BLAS_LIBRARY_DIRS}")
  			endif()
  			list(APPEND REQUIRED_LIBS "${BLAS_LIBRARIES}")
  			if(BLAS_LINKER_FLAGS)
  				list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
  			endif()
  			
  			# Set required libraries for linking
  			set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  			set(CMAKE_REQUIRED_LIBRARIES)
  			list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
  			foreach(lib_dir ${REQUIRED_LIBDIRS})
  				list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  			endforeach()
  			list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  			string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
  			
  			# Test the Link
  			unset(CBLAS_WORKS CACHE)
      		check_function_exists(cblas_dscal CBLAS_WORKS)
  			mark_as_advanced(CBLAS_WORKS)
  			
  			# If it worked then set final variables
  			if(CBLAS_WORKS)
  				set(CBLAS_LIBRARIES_DEP "${REQUIRED_LIBS}")
  				set(CBLAS_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
  				set(CBLAS_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
				set(CBLAS_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
				list(REMOVE_DUPLICATES CBLAS_LIBRARY_DIRS_DEP)
				list(REMOVE_DUPLICATES CBLAS_INCLUDE_DIRS_DEP)
				list(REMOVE_DUPLICATES CBLAS_LINKER_FLAGS)
  			else()
  				if(NOT CBLAS_FIND_QUIETLY)
  					message(STATUS "Searching for CBLAS: Tests of cblas_dscal with cblas and blas libraries failed")
  					message(STATUS "Attemped Includes:  ${CMAKE_REQUIRED_INCLUDES}")
  					message(STATUS "Attemped Libraries: ${CMAKE_REQUIRED_LIBRARIES}")   
  					message(STATUS "Check in CMakeFiles/CMakeError.log to find more details")
  				endif()
  			endif()
  			set(CMAKE_REQUIRED_INCLUDES)
  			set(CMAKE_REQUIRED_FLAGS)
  			set(CMAKE_REQUIRED_LIBRARIES)
  		endif(CBLAS_LIBRARIES)
  	
  	endif(NOT CBLAS_WORKS)
  	
else(BLAS_FOUND)
	message(STATUS "CBLAS requires BLAS but BLAS was not been found.")
	message(STATUS "You must find BLAS first.")
endif(BLAS_FOUND)

mark_as_advanced(CBLAS_DIR) 
mark_as_advanced(CBLAS_DIR_FOUND)

# Set CBLAS_FOUND
include(FindPackageHandleStandardArgs) 
find_package_handle_standard_args(CBLAS DEFAULT_MSG   
	CBLAS_LIBRARIES   
	CBLAS_WORKS
)

if(CBLAS_FOUND)
	if(NOT CBLAS_FIND_QUIETLY)
	       message(STATUS "CBLAS_LIBRARIES_DEP    = ${CBLAS_LIBRARIES_DEP}")
	       message(STATUS "CBLAS_LIBRARY_DIRS_DEP = ${CBLAS_LIBRARY_DIRS_DEP}")
	       message(STATUS "CBLAS_INCLUDE_DIRS_DEP = ${CBLAS_INCLUDE_DIRS_DEP}")
	       message(STATUS "CBLAS_LINKER_FLAGS     = ${CBLAS_LINKER_FLAGS}")
	endif()
endif()


# Export Modern Library
if(CBLAS_FOUND AND NOT TARGET CBLAS::CBLAS)
  add_library(CBLAS::CBLAS UNKNOWN IMPORTED)
  set_target_properties(CBLAS::CBLAS PROPERTIES
#    IMPORTED_LOCATION "${CBLAS_LIBRARIES_DEP}"
    INTERFACE_LINK_OPTIONS "${CBLAS_LINKER_FLAGS}"
    INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIRS_DEP}"
    INTERFACE_LINK_LIBRARIES "${CBLAS_LIBRARIES_DEP}"
  )
endif()

