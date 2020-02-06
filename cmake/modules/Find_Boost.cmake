#
# Find Boost Headers & Library 
#
# This wrapper uses the CMake provided function for
# locating the Boost Library. The following 
# details are an short description of how to use it.
#
# Imported Targets
# =======================
#
# Boost::headers - Target for header-only dependencies
# Boost::<C>     - Target for specific component dependency <C> is lower-case
#
# Result Variables
# =======================
#
# Boost_FOUND        = True if headers and requested libraries were found
# Boost_INCLUDE_DIRS = Boost include directories
# Boost_LIBRARY_DIRS = Link directories for Boost libraries
# Boost_LIBRARIES    = Boost component libraries to be linked
# Boost_<C>_FOUND    = True if component <C> was found (<C> is upper-case)
# Boost_<C>_LIBRARY  = Libraries to link for component <C>
# Boost_VERSION      = same as Boost_VERSION_STRING
#
# Boost_INCLUDE_DIR         = Directory containing Boost headers
# Boost_LIBRARY_DIR_RELEASE = Directory containing release Boost libraries
# Boost_LIBRARY_DIR_DEBUG   = Directory containing debug Boost libraries
# Boost_<C>_LIBRARY_DEBUG   = Component <C> library debug variant
# Boost_<C>_LIBRARY_RELEASE = Component <C> library release variant
#
# Variables To Locate Boost
# =======================
#
# BOOST_ROOT       = Preferred installation prefix
# BOOST_INCLUDEDIR = Preferred include directory e.g. <prefix>/include
# BOOST_LIBRARYDIR = Preferred library directory e.g. <prefix>/lib
#
#

message(VERBOSE "")
message(VERBOSE "--------------------- Boost Libraries ------------------------")
message(VERBOSE "Search Locations:")
message(VERBOSE "BOOST_ROOT       = ${BOOST_ROOT}")
message(VERBOSE "BOOST_INCLUDEDIR = ${BOOST_INCLUDEDIR}")
message(VERBOSE "BOOST_LIBRARYDIR = ${BOOST_LIBRARYDIR}")


set(Boost_DEBUG                  OFF) # Enable debug output from FIND_PACKAGE
set(Boost_NO_SYSTEM_PATHS        ON)  # Do not search system paths before BOOST_ROOT
set(Boost_USE_DEBUG_LIBS         OFF) # Use debug libs
set(Boost_USE_RELEASE_LIBS       ON)  # Use release libs
set(Boost_USE_MULTITHREADED      ON)  # Use Boost multi-threaded code
set(Boost_USE_STATIC_LIBS        ON)  # Static link to Boost libraries
set(Boost_USE_STATIC_RUNTIME     OFF) # Linked statically to the C++ runtime
set(Boost_USE_DEBUG_RUNTIME      OFF) # Linked to the MS debug C++ runtime

find_package(Boost COMPONENTS mpi serialization REQUIRED)


message("")
message(VERBOSE "Results:")
message(VERBOSE "Boost Found               = ${Boost_FOUND}")
message(VERBOSE "Boost Version             = ${Boost_VERSION}")
message(VERBOSE "Boost Includes            = ${Boost_INCLUDE_DIRS}")
message(VERBOSE "Boost Link Libraries      = ${Boost_LIBRARY_DIRS}")
message(VERBOSE "Boost Component Libraries = ${Boost_LIBRARIES}")

