cmake_minimum_required(VERSION 3.20)

project(fhhos4 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# program options
option(ENABLE_GMSH "Build with the GMSH interface" ON)
option(ENABLE_AGMG "Build with AGMG" OFF)
option(ENABLE_CGAL "Build with CGAL" ON)
option(SMALL_INDEX "Use type 'int' for the indexes" ON)
option(ENABLE_1D "Enable 1D" OFF)
option(ENABLE_2D "Enable 2D" ON)
option(ENABLE_3D "Enable 3D" ON)

set(EXECUTABLE_OUTPUT_PATH bin/)

# Pass -DCMAKE_VERBOSE_MAKEFILE=ON (or `cmake --build build -- VERBOSE=1`) to see full compiler commands

# Files to compile
file(
	GLOB_RECURSE
	source_files
	CONFIGURE_DEPENDS
	src/*
)

# Include AGMG compiled files
if(ENABLE_AGMG)
	message(STATUS "AGMG enabled")
	set(OBJS ${AGMG_DIR}/dagmg.o ${AGMG_DIR}/dagmg_mumps.o)
else()
	message(STATUS "AGMG disabled. Add -DENABLE_AGMG=ON and -DAGMG_DIR=<path/to/.o/files> to enable it.")
endif()

# Executable file
add_executable(fhhos4 ${OBJS} ${source_files})

target_compile_definitions(fhhos4 PRIVATE ROOT_PATH="${CMAKE_SOURCE_DIR}")

if(SMALL_INDEX)
	target_compile_definitions(fhhos4 PRIVATE SMALL_INDEX)
	message(STATUS "SMALL_INDEX: ON (int)")
else()
	message(STATUS "SMALL_INDEX: OFF (size_t)")
endif()

if(ENABLE_1D)
	target_compile_definitions(fhhos4 PRIVATE ENABLE_1D)
	message(STATUS "ENABLE_1D: ON")
else()
	message(STATUS "ENABLE_1D: OFF")
endif()

if(ENABLE_2D)
	target_compile_definitions(fhhos4 PRIVATE ENABLE_2D)
	message(STATUS "ENABLE_2D: ON")
else()
	message(STATUS "ENABLE_2D: OFF")
endif()

if(ENABLE_3D)
	target_compile_definitions(fhhos4 PRIVATE ENABLE_3D)
	message(STATUS "ENABLE_3D: ON")
else()
	message(STATUS "ENABLE_3D: OFF")
endif()

if(ENABLE_AGMG)
	target_compile_definitions(fhhos4 PRIVATE AGMG_ENABLED)
endif()

# Include Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(fhhos4 PRIVATE Eigen3::Eigen)

# Include CGAL
if(ENABLE_CGAL)
	message(STATUS "CGAL enabled")
	# CGAL's own exported CMake files call cmake_minimum_required() internally, which resets
	# CMP0167 right before their own find_package(Boost) call and triggers a policy warning.
	# A plain cmake_policy(SET ...) gets clobbered by that nested call, so pin the default here
	# instead - it's rechecked on every cmake_minimum_required()/cmake_policy(VERSION) call.
	if(POLICY CMP0167)
		set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
	endif()
	find_package(CGAL REQUIRED)
	target_link_libraries(fhhos4 PRIVATE CGAL::CGAL)
	target_compile_definitions(fhhos4 PRIVATE CGAL_ENABLED)
else()
	message(STATUS "CGAL disabled. Add -DENABLE_CGAL=ON to enable it.")
endif()

# Include GMSH
if(ENABLE_GMSH)
	message(STATUS "GMSH enabled")
	if(GMSH_API)
		target_include_directories(fhhos4 PRIVATE ${GMSH_API})
	else()
		find_path(GMSH_API NAMES gmsh.h)
		if(GMSH_API)
			target_include_directories(fhhos4 PRIVATE ${GMSH_API})
			message(STATUS "GMSH header found: ${GMSH_API}/gmsh.h")
		else()
			message(WARNING "Option -DGMSH_API missing, and gmsh.h could not be found automatically. Set this option if gmsh.h is not found during compilation")
		endif()
	endif()
	if(GMSH_LIB)
		target_link_libraries(fhhos4 PRIVATE "${GMSH_LIB}") # -lgmsh
	else()
		find_library(GMSH_LIB NAMES gmsh)
		if(GMSH_LIB)
			target_link_libraries(fhhos4 PRIVATE "${GMSH_LIB}") # -lgmsh
			message(STATUS "GMSH library found: ${GMSH_LIB}")
		else()
			message(FATAL_ERROR "The GMSH dynamic library was not found. Note that GMSH must be compiled with the option -DENABLE_BUILD_DYNAMIC=1. Pass the library path with the option -DGMSH_LIB=<path>/libgmsh.so.")
		endif()
	endif()
	target_compile_definitions(fhhos4 PRIVATE GMSH_ENABLED)
else()
	message(STATUS "GMSH disabled. Add -DENABLE_GMSH=ON to enable it.")
endif()

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(fhhos4 PRIVATE Threads::Threads)

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

	target_link_libraries(fhhos4 PRIVATE "-lgfortran" "-lm" ${TPL_LAPACK_LIB} ${TPL_BLAS_LIB})
endif()
