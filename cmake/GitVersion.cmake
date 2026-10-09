# Writes OUTPUT (the definition of fhhos4::Version()) from library/Version.cpp.in: VERSION, followed by the commit
# (git describe) when SOURCE_DIR is a git clone. A release archive has no .git: VERSION alone.
# Usage: cmake -DSOURCE_DIR=<dir> -DVERSION=<x.y.z> -DGIT_EXECUTABLE=<git> -DOUTPUT=<file> -P GitVersion.cmake
set(FULL_VERSION ${VERSION})
if(GIT_EXECUTABLE AND EXISTS ${SOURCE_DIR}/.git)
	execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
		WORKING_DIRECTORY ${SOURCE_DIR}
		OUTPUT_VARIABLE describe OUTPUT_STRIP_TRAILING_WHITESPACE
		RESULT_VARIABLE result ERROR_QUIET)
	if(result EQUAL 0)
		set(FULL_VERSION "${VERSION} (${describe})")
	endif()
endif()
# configure_file() leaves OUTPUT untouched when its content doesn't change: no recompilation
configure_file(${SOURCE_DIR}/library/Version.cpp.in ${OUTPUT} @ONLY)
