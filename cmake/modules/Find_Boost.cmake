#
# Find Boost Library
# 
# Find and configure for Boost usage
#
#

if(USE_BOOST)
message(VERBOSE "")
message(VERBOSE "--------------------- Boost Libraries ------------------------")
message(VERBOSE "Search Locations:")
message(VERBOSE "BOOST_ROOT            = ${BOOST_ROOT}")
message(VERBOSE "BOOST_INCLUDE_DIR     = ${BOOST_INCLUDE_DIR}")
message(VERBOSE "BOOST_INCLUDE_LIBRARY = ${BOOST_INCLUDE_LIBRARY}")


set(Boost_DEBUG OFF)              # Enable debug output from FIND_PACKAGE
set(Boost_NO_SYSTEM_PATHS ON)     # Do not search system paths before BOOST_ROOT
set(Boost_USE_STATIC_LIBS ON)     # Static link to Boost libraries
set(Boost_USE_STATIC_RUNTIME OFF) # Use Boost static linked to C++ runtime
#set(Boost_USE_MULTITHREADED OFF)  # Use Boost multi-threaded code
find_package(Boost 1.65.0 COMPONENTS mpi serialization REQUIRED)


message("")
message(VERBOSE "Results:")
message(VERBOSE "Boost Found               = ${Boost_FOUND}")
message(VERBOSE "Boost Version             = ${Boost_VERSION}")
message(VERBOSE "Boost Includes            = ${Boost_INCLUDE_DIRS}")
message(VERBOSE "Boost Link Libraries      = ${Boost_LIBRARY_DIRS}")
message(VERBOSE "Boost Component Libraries = ${Boost_LIBRARIES}")


endif(USE_BOOST)

