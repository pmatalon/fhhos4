cmake_minimum_required (VERSION 3.1)

project(dghho)

# enable C++14 (for CGAL)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_GMSH "Build with the GMSH interface" ON)
option(ENABLE_AGMG "Build with AGMG" OFF)

set(EXECUTABLE_OUTPUT_PATH bin/)

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

add_definitions(-DROOT_PATH="${CMAKE_SOURCE_DIR}")

# Include CGAL
find_package(CGAL)
target_link_libraries(dghho CGAL::CGAL)

# Include Eigen
include_directories(eigen-eigen-b3f3d4950030)

# Include GMSH
if (ENABLE_GMSH)
	include_directories(${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/include)
	find_library(GMSH_LIB NAMES gmsh PATHS "${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib")
	if(GMSH_LIB)
		target_link_libraries(dghho "${GMSH_LIB}") # -lgmsh
		add_definitions(-DGMSH_ENABLED)
	else()
		message(FATAL_ERROR "GMSH dynamic libraries not found. GMSH must be compiled with the option -DENABLE_BUILD_DYNAMIC=1 and the libraries put in a standard directory or in ${CMAKE_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib")
	endif()
	message(STATUS "GMSH enabled")
else()
	message(STATUS "GMSH disabled. Add -DENABLE_GMSH=ON to enable it.")
endif()

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(dghho Threads::Threads)

# c++ flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g") # uncomment this line and comment the one above to use Valgrind

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
	message(STATUS "AGMG disabled. Add -DENABLE_AGMG=ON to enable it.")
endif()
