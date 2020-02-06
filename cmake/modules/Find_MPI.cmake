#
# Find MPI Compiler & Features
#
# This wrapper uses the CMake provided function for
# locating the MPI compiler and options. The following 
# details are an short description of how to use it.
#
# Languages <lang>
# =======================
# C       = MPI C API
# CXX     = MPI C API used by C++ code
# MPICXX  = MPI C++ API which was removed in MPI-3
# Fortran = MPI Fortran API
#
# Imported Targets
# =======================
#
# MPI::MPI_<lang> - Target for using MPI from <lang>
#
# Result Variables
# =======================
#
# MPI_<lang>_FOUND           = MPI settings for <lang> were found
# MPI_<lang>_COMPILER        = MPI compiler for <lang>
# MPI_<lang>_COMPILE_OPTIONS = Compilation options for MPI programs in <lang>
# MPI_<lang>_INCLUDE_DIRS    = Include path(s) for MPI header
# MPI_<lang>_LINK_FLAGS      = Linker flags for MPI programs
# MPI_<lang>_LIBRARIES       = All libraries to link MPI programs against
# MPI_<lang>_VERSION         = MPI version implemented for <lang>
#
# Variables To Locate MPI
# =======================
#
# MPIEXEC_EXECUTABLE = Manually specify the location of mpiexec
# MPI_HOME           = Specify the base directory of the MPI installation
# ENV{MPI_HOME}      = Environment variable of MPI_HOME
#
# IMPORTANT:
# =======================
# Calling the builtin CMake function follows the four step process.
#
# 1. Search for MPIEXEC_EXECUTABLE and, if found, use its base directory.
# 2. Search for environment variables MPI_<lang>_COMPILER
# 3. Attempt to find an MPI compiler wrapper on system
# 4. Try to find an MPI that does not ship (Windows)
#
# If all goes as planned MPI code can be compiled 
# and linked using the CMake Properties 
#
# target_link_libraries(MyTarget PUBLIC MPI::MPI_C)
# target_link_libraries(MyTarget PUBLIC MPI::MPI_CXX)
# target_link_libraries(MyTarget PUBLIC MPI::MPI_Fortran)
#
#
#

# If true, the module assumes that the compiler itself does not provide an MPI
set(MPI_ASSUME_NO_BUILTIN_MPI FALSE)

# If true, no compiler wrapper will be searched for
set(MPI_SKIP_COMPILER_WRAPPER FALSE)

# If true, the guessing step will be skipped.
set(MPI_SKIP_GUESSING FALSE)

# If true, disable the MPI-2 C++ bindings
set(MPI_CXX_SKIP_MPICXX FALSE)

message(VERBOSE "")
message(VERBOSE "--------------------- MPI Libraries ------------------------")
#
# First attempt to use the CMake tools to locate the library
#
find_package(MPI REQUIRED)

# For supporting CMake < 3.9:
if(NOT TARGET MPI::MPI_C)
	add_library(MPI::MPI_C IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_C_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_C
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES})
endif()
if(NOT TARGET MPI::MPI_CXX)
	add_library(MPI::MPI_CXX IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif()
if(NOT TARGET MPI::MPI_MPICXX)
	add_library(MPI::MPI_MPICXX IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_MPICXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_MPICXX_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_MPICXX
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_MPICXX_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_MPICXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_MPICXX_LINK_FLAGS} ${MPI_MPICXX_LIBRARIES})
endif()
if(NOT TARGET MPI::MPI_Fortran)
	add_library(MPI::MPI_Fortran IMPORTED INTERFACE)
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_Fortran_COMPILE_OPTIONS})
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_Fortran_INCLUDE_DIRS}")
    set_property(TARGET MPI::MPI_Fortran
                 PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES})
endif()


#
# Report what we found
#
message(VERBOSE "")
message(VERBOSE "MPI Found      = ${MPI_FOUND}")
message(VERBOSE "")
message(VERBOSE "C")
message(VERBOSE " - Found       = ${MPI_C_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_C_COMPILER}")
message(VERBOSE " - Options     = ${MPI_C_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_C_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_C_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_C_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_C_LIBRARIES}")
message(VERBOSE "C++")
message(VERBOSE " - Found       = ${MPI_CXX_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_CXX_COMPILER}")
message(VERBOSE " - Options     = ${MPI_CXX_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_CXX_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_CXX_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_CXX_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_CXX_LIBRARIES}")
message(VERBOSE "Fortran:")
message(VERBOSE " - Found       = ${MPI_Fortran_FOUND}")
message(VERBOSE " - Compiler    = ${MPI_Fortran_COMPILER}")
message(VERBOSE " - Options     = ${MPI_Fortran_COMPILE_OPTIONS}")
message(VERBOSE " - Definitions = ${MPI_Fortran_COMPILE_DEFINITIONS}")
message(VERBOSE " - Include Dir = ${MPI_Fortran_INCLUDE_DIRS}")
message(VERBOSE " - Link Flags  = ${MPI_Fortran_LINK_FLAGS}")
message(VERBOSE " - Libraries   = ${MPI_Fortran_LIBRARIES}")
