#
# Set CMake variables for code record
#
# Once done this will define
#
# GIT_SHORT_HASH  - Short hash for code being configured
# GIT_LONG_HASH   - Long hash for code being configured
# GIT_COMMIT_TIME - Long hash for code being configured
# CONFIGURE_TIME  - Time on system for configuration
#


#
# Test for git and get results
#
find_program(GIT_EXECUTABLE git)

if( GIT_EXECUTABLE )

execute_process(
	COMMAND
		git rev-parse --short HEAD
    RESULT_VARIABLE
        GIT_SHORT_HASH_RESULT
    OUTPUT_VARIABLE
        GIT_SHORT_HASH
)
string(REGEX REPLACE "\n" "" GIT_SHORT_HASH "${GIT_SHORT_HASH}")

execute_process(
	COMMAND
		git rev-parse HEAD
    RESULT_VARIABLE
        GIT_LONG_HASH_RESULT
    OUTPUT_VARIABLE
        GIT_LONG_HASH
)
string(REGEX REPLACE "\n" "" GIT_LONG_HASH "${GIT_LONG_HASH}")

execute_process(
	COMMAND
		git show -s --format=%ci HEAD
    RESULT_VARIABLE
        GIT_COMMIT_TIME_RESULT
    OUTPUT_VARIABLE
        GIT_COMMIT_TIME
)
string(REGEX REPLACE "\n" "" GIT_COMMIT_TIME "${GIT_COMMIT_TIME}")

endif(GIT_EXECUTABLE)

#
# Get the time when CMake was run
#
string(TIMESTAMP CONFIGURE_TIME "%Y-%m-%d %H:%M:%S +0000" UTC)

#
# Report if requested
#
message(VERBOSE "")
message(VERBOSE "--------------------- Version Report -----------------------")
message(VERBOSE "")
message(VERBOSE "PROJECT_VERSION = ${PROJECT_VERSION}")
message(VERBOSE "CONFIGURE_TIME  = ${CONFIGURE_TIME}")
message(VERBOSE "GIT_COMMIT_TIME = ${GIT_COMMIT_TIME}")
message(VERBOSE "GIT_SHORT_HASH  = ${GIT_SHORT_HASH}")
message(VERBOSE "GIT_LONG_HASH   = ${GIT_LONG_HASH}")


