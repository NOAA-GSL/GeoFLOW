
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

if (BLAS_FOUND)

  if (NOT CBLAS_STANDALONE)
    # check if a cblas function exists in the BLAS lib
    # this can be the case with libs such as MKL, ACML
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LINKER_FLAGS};${BLAS_LIBRARIES}")
    set(CMAKE_REQUIRED_FLAGS "${BLAS_COMPILER_FLAGS}")
    unset(CBLAS_WORKS CACHE)
    check_function_exists(cblas_dscal CBLAS_WORKS)
    check_function_exists(cblas_zgemm3m CBLAS_ZGEMM3M_FOUND)
    mark_as_advanced(CBLAS_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)

    if(CBLAS_WORKS)

      # Check for faster complex GEMM routine
      # (only C/Z, no S/D version)
      if ( CBLAS_ZGEMM3M_FOUND )
        add_definitions(-DCBLAS_HAS_ZGEMM3M -DCBLAS_HAS_CGEMM3M)
      endif()

      if(NOT CBLAS_FIND_QUIETLY)
        message(STATUS "Looking for cblas: test with blas succeeds")
      endif()
      # test succeeds: CBLAS is in BLAS
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
    endif()
  endif (NOT CBLAS_STANDALONE)

  if (CBLAS_STANDALONE OR NOT CBLAS_WORKS)

    if(NOT CBLAS_WORKS AND NOT CBLAS_FIND_QUIETLY)
      message(STATUS "Looking for cblas : test with blas fails")
    endif()
    # test fails: try to find CBLAS lib exterior to BLAS

    # Try to find CBLAS lib
    #######################

    # Looking for include
    # -------------------

    # Add system include paths to search include
    # ------------------------------------------
    unset(_inc_env)
    set(ENV_CBLAS_DIR "$ENV{CBLAS_DIR}")
    set(ENV_CBLAS_INCDIR "$ENV{CBLAS_INCDIR}")
    if(ENV_CBLAS_INCDIR)
      list(APPEND _inc_env "${ENV_CBLAS_INCDIR}")
    elseif(ENV_CBLAS_DIR)
      list(APPEND _inc_env "${ENV_CBLAS_DIR}")
      list(APPEND _inc_env "${ENV_CBLAS_DIR}/include")
      list(APPEND _inc_env "${ENV_CBLAS_DIR}/include/cblas")
    else()
      if(WIN32)
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        list(APPEND _inc_env "${_path_env}")
      else()
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
        list(APPEND _inc_env "${_path_env}")
        string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
        list(APPEND _inc_env "${_path_env}")
      endif()
    endif()
    list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(REMOVE_DUPLICATES _inc_env)


    # Try to find the cblas header in the given paths
    # -------------------------------------------------
    # call cmake macro to find the header path
    if(CBLAS_INCDIR)
      set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
      find_path(CBLAS_cblas.h_DIRS
        NAMES cblas.h
        HINTS ${CBLAS_INCDIR})
    else()
      if(CBLAS_DIR)
        set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
        find_path(CBLAS_cblas.h_DIRS
          NAMES cblas.h
          HINTS ${CBLAS_DIR}
          PATH_SUFFIXES "include" "include/cblas")
      else()
        set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
        find_path(CBLAS_cblas.h_DIRS
          NAMES cblas.h
          HINTS ${_inc_env}
          PATH_SUFFIXES "cblas")
      endif()
    endif()
    mark_as_advanced(CBLAS_cblas.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (CBLAS_cblas.h_DIRS)
      set(CBLAS_INCLUDE_DIRS "${CBLAS_cblas.h_DIRS}")
    else ()
      set(CBLAS_INCLUDE_DIRS "CBLAS_INCLUDE_DIRS-NOTFOUND")
      if(NOT CBLAS_FIND_QUIETLY)
        message(STATUS "Looking for cblas -- cblas.h not found")
      endif()
    endif()


    # Looking for lib
    # ---------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_CBLAS_LIBDIR "$ENV{CBLAS_LIBDIR}")
    if(ENV_CBLAS_LIBDIR)
      list(APPEND _lib_env "${ENV_CBLAS_LIBDIR}")
    elseif(ENV_CBLAS_DIR)
      list(APPEND _lib_env "${ENV_CBLAS_DIR}")
      list(APPEND _lib_env "${ENV_CBLAS_DIR}/lib")
    else()
      if(WIN32)
        string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
      else()
        if(APPLE)
          string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
        else()
          string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
        endif()
        list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
        list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
      endif()
    endif()
    list(REMOVE_DUPLICATES _lib_env)

    # Try to find the cblas lib in the given paths
    # ----------------------------------------------

    # call cmake macro to find the lib path
    if(CBLAS_LIBDIR)
      set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
      find_library(CBLAS_cblas_LIBRARY
        NAMES cblas
        HINTS ${CBLAS_LIBDIR})
    else()
      if(CBLAS_DIR)
        set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
        find_library(CBLAS_cblas_LIBRARY
          NAMES cblas
          HINTS ${CBLAS_DIR}
          PATH_SUFFIXES lib lib32 lib64)
      else()
        set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
        find_library(CBLAS_cblas_LIBRARY
          NAMES cblas
          HINTS ${_lib_env})
      endif()
    endif()
    mark_as_advanced(CBLAS_cblas_LIBRARY)

    # If found, add path to cmake variable
    # ------------------------------------
    if (CBLAS_cblas_LIBRARY)
      get_filename_component(cblas_lib_path "${CBLAS_cblas_LIBRARY}" PATH)
      # set cmake variables
      set(CBLAS_LIBRARIES    "${CBLAS_cblas_LIBRARY}")
      set(CBLAS_LIBRARY_DIRS "${cblas_lib_path}")
    else ()
      set(CBLAS_LIBRARIES    "CBLAS_LIBRARIES-NOTFOUND")
      set(CBLAS_LIBRARY_DIRS "CBLAS_LIBRARY_DIRS-NOTFOUND")
      if (NOT CBLAS_FIND_QUIETLY)
        message(STATUS "Looking for cblas -- lib cblas not found")
      endif()
    endif ()

    # check a function to validate the find
    if(CBLAS_LIBRARIES)

      set(REQUIRED_INCDIRS)
      set(REQUIRED_LDFLAGS)
      set(REQUIRED_LIBDIRS)
      set(REQUIRED_LIBS)

      # CBLAS
      if (CBLAS_INCLUDE_DIRS)
        set(REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
      endif()
      if (CBLAS_LIBRARY_DIRS)
        set(REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
      endif()
      set(REQUIRED_LIBS "${CBLAS_LIBRARIES}")
      # BLAS
      if (BLAS_INCLUDE_DIRS)
        list(APPEND REQUIRED_INCDIRS "${BLAS_INCLUDE_DIRS}")
      endif()
      if (BLAS_LIBRARY_DIRS)
        list(APPEND REQUIRED_LIBDIRS "${BLAS_LIBRARY_DIRS}")
      endif()
      list(APPEND REQUIRED_LIBS "${BLAS_LIBRARIES}")
      if (BLAS_LINKER_FLAGS)
        list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
      endif()

      # set required libraries for link
      set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
      set(CMAKE_REQUIRED_LIBRARIES)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
      foreach(lib_dir ${REQUIRED_LIBDIRS})
        list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
      endforeach()
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
      string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

      # test link
      unset(CBLAS_WORKS CACHE)
      include(CheckFunctionExists)
      check_function_exists(cblas_dscal CBLAS_WORKS)
      mark_as_advanced(CBLAS_WORKS)

      if(CBLAS_WORKS)

        # Check for faster complex GEMM routine
        # (only C/Z, no S/D version)
        check_function_exists(cblas_zgemm3m CBLAS_ZGEMM3M_FOUND)
        if ( CBLAS_ZGEMM3M_FOUND )
          add_definitions(-DCBLAS_HAS_ZGEMM3M -DCBLAS_HAS_CGEMM3M)
        endif()

        # save link with dependencies
        set(CBLAS_LIBRARIES_DEP "${REQUIRED_LIBS}")
        set(CBLAS_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
        set(CBLAS_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
        set(CBLAS_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
        list(REMOVE_DUPLICATES CBLAS_LIBRARY_DIRS_DEP)
        list(REMOVE_DUPLICATES CBLAS_INCLUDE_DIRS_DEP)
        list(REMOVE_DUPLICATES CBLAS_LINKER_FLAGS)
      else()
        if(NOT CBLAS_FIND_QUIETLY)
          message(STATUS "Looking for cblas : test of cblas_dscal with cblas and blas libraries fails")
          message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
          message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
          message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
        endif()
      endif()
      set(CMAKE_REQUIRED_INCLUDES)
      set(CMAKE_REQUIRED_FLAGS)
      set(CMAKE_REQUIRED_LIBRARIES)
    endif(CBLAS_LIBRARIES)

  endif (CBLAS_STANDALONE OR NOT CBLAS_WORKS)

else(BLAS_FOUND)

  if (NOT CBLAS_FIND_QUIETLY)
    message(STATUS "CBLAS requires BLAS but BLAS has not been found."
      "Please look for BLAS first.")
  endif()

endif(BLAS_FOUND)

if (CBLAS_LIBRARIES)
  list(GET CBLAS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(CBLAS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of CBLAS library" FORCE)
  else()
    set(CBLAS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of CBLAS library" FORCE)
  endif()
endif()
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
  add_library(CBLAS::CBLAS INTERFACE IMPORTED)
  set_target_properties(CBLAS::CBLAS PROPERTIES
    IMPORTED_LOCATION "${CBLAS_LIBRARIES_DEP}"
    INTERFACE_LINK_OPTIONS "${CBLAS_LINKER_FLAGS}"
    INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIRS_DEP}"
    INTERFACE_LINK_LIBRARIES "${CBLAS_LIBRARIES_DEP}"
  )
endif()

#if(CBLAS_FOUND AND NOT TARGET CBLAS::CBLAS)
#  add_library(CBLAS::CBLAS INTERFACE IMPORTED)
#  set_target_properties(CBLAS::CBLAS PROPERTIES
#    INTERFACE_LINK_LIBRARIES "${CBLAS_LIBRARIES_DEP}"
#  )
#endif()

