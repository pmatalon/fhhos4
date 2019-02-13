cmake_minimum_required (VERSION 2.6)

project(dghho)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

# Include Eigen
include_directories(eigen-eigen-b3f3d4950030)

file(
	GLOB_RECURSE
	source_files
	src/*
)

add_executable(dghho ${source_files})