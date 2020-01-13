cmake_minimum_required (VERSION 2.8)

project(dghho)

option(ENABLE_GMSH "Build with the GMSH interface" OFF)
option(ENABLE_AGMG "Build with the AGMG" OFF)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

# show compiler commands
set(CMAKE_VERBOSE_MAKEFILE ON)

# Executable file
file(
	GLOB_RECURSE
	source_files
	src/*
)
add_executable(dghho ${source_files})
#add_executable(dghho "${PROJECT_SOURCE_DIR}/src/main.cpp")

# Include Eigen
include_directories(eigen-eigen-b3f3d4950030)

# Include GMSH
if (ENABLE_GMSH)
	include_directories(${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/include)
	find_library(GMSH_LIB NAMES gmsh PATHS "${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib")
	if(GMSH_LIB)
		#link_directories(${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib)
		#target_link_libraries(dghho "-lgmsh")
		target_link_libraries(dghho "${GMSH_LIB}")
		add_definitions(-DGMSH_ENABLED)
	else()
		message(FATAL_ERROR "GMSH not found. Must be compiled and put in ${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib")
	endif()
	message(STATUS "GMSH enabled")
else()
	message(STATUS "GMSH disabled")
endif()

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(dghho Threads::Threads)

# c++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

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

	find_library(TPL_AGMG_LIB NAMES agmg PATHS ${TPL_AGMG_LIB_PATH})
	if(TPL_AGMG_LIB)
		message(STATUS "Found AGMG: ${TPL_AGMG_LIB}")
	else()
		message(FATAL_ERROR "AGMG not found.")
	endif()

	target_link_libraries(dghho "-lgfortran" "-lm" ${TPL_AGMG_LIB} ${TPL_LAPACK_LIB} ${TPL_BLAS_LIB})

	add_definitions(-DAGMG_ENABLED)

	message(STATUS "AGMG enabled")
else()
	message(STATUS "AGMG disabled")
endif()
