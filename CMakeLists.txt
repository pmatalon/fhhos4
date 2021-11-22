cmake_minimum_required (VERSION 3.1)

project(fhhos4)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# c++ flags
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast") # comment this line to have exact same results than in Debug

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g") # uncomment this line and comment the one above to use Valgrind

# enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# program options
option(ENABLE_GMSH "Build with the GMSH interface" ON)
option(ENABLE_AGMG "Build with AGMG" OFF)
option(ENABLE_CGAL "Build with CGAL" ON)
option(SMALL_INDEX "Use type 'int' for the indexes" OFF)

add_definitions(-DROOT_PATH="${CMAKE_SOURCE_DIR}")

if(SMALL_INDEX)
	add_definitions(-DSMALL_INDEX)
	message(STATUS "SMALL_INDEX: ON (int)")
else()
	message(STATUS "SMALL_INDEX: OFF (size_t)")
endif()

set(EXECUTABLE_OUTPUT_PATH bin/)

# show compiler commands
set(CMAKE_VERBOSE_MAKEFILE ON)

# Files to compile
file(
	GLOB_RECURSE
	source_files
	src/*
)

# Include AGMG compiled files
if(ENABLE_AGMG)
	message(STATUS "AGMG enabled")
	SET(OBJS ${AGMG_DIR}/dagmg.o ${AGMG_DIR}/dagmg_mumps.o)
	add_definitions(-DAGMG_ENABLED)
else()
	message(STATUS "AGMG disabled. Add -DENABLE_AGMG=ON and -DAGMG_DIR=<path/to/.o/files> to enable it.")
endif()

# Executable file
add_executable(fhhos4 ${OBJS} ${source_files})

# Include CGAL
if(ENABLE_CGAL)
	message(STATUS "CGAL enabled")
	set(CGAL_DIR "${CMAKE_SOURCE_DIR}/dependencies/cgal/CGAL-5.1")
	find_package(CGAL)
	target_link_libraries(fhhos4 CGAL::CGAL)
	add_definitions(-DCGAL_ENABLED)
else()
	message(STATUS "CGAL disabled. Add -DENABLE_CGAL=ON to enable it.")
endif()

# Include Eigen
include_directories(dependencies/eigen/eigen-eigen-b3f3d4950030)

# Include GMSH
if (ENABLE_GMSH)
	message(STATUS "GMSH enabled")
	if (GMSH_API)
		include_directories(${GMSH_API})
	else()
		message(WARNING "Option -DGMSH_API missing. Set this option if gmsh.h is not found during compilation")
	endif()
	if (GMSH_LIB)
		target_link_libraries(fhhos4 "${GMSH_LIB}") # -lgmsh
	else()
		find_library(GMSH_LIB NAMES gmsh)
		if(GMSH_LIB)
			target_link_libraries(fhhos4 "${GMSH_LIB}") # -lgmsh
			message(STATUS "GMSH library found: ${GMSH_LIB}")
		else()
			message(FATAL_ERROR "The GMSH dynamic library was not found. Note that GMSH must be compiled with the option -DENABLE_BUILD_DYNAMIC=1. Pass the library path with the option -DGMSH_LIB=<path>/libgmsh.so.")
		endif()
	endif()
	add_definitions(-DGMSH_ENABLED)
else()
	message(STATUS "GMSH disabled. Add -DENABLE_GMSH=ON to enable it.")
endif()

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(fhhos4 Threads::Threads)

# Include AGMG (requires Blas and Lapack)
if(ENABLE_AGMG)
	find_library(TPL_BLAS_LIB NAMES blas openblas PATHS ${TPL_BLAS_LIB_PATH})
	if(TPL_BLAS_LIB)
		message(STATUS "Found BLAS: ${TPL_BLAS_LIB}")
	else()
		message(FATAL_ERROR "BLAS not found.")
	endif()

	find_library(TPL_LAPACK_LIB NAMES lapack PATHS ${TPL_LAPACK_LIB_PATH})
	if(TPL_LAPACK_LIB)
		message(STATUS "Found LAPACK: ${TPL_LAPACK_LIB}")
	else()
		message(FATAL_ERROR "LAPACK not found.")
	endif()

	target_link_libraries(fhhos4 "-lgfortran" "-lm" ${TPL_LAPACK_LIB} ${TPL_BLAS_LIB})
endif()
