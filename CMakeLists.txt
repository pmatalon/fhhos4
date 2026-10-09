cmake_minimum_required(VERSION 3.20)

# The version of fhhos4, written here only: fhhos4_version.h, fhhos4::Version() and the CMake package are generated from
# it. A release changes it, then gets the annotated tag v<version>.
project(fhhos4 VERSION 1.2.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# program options
option(ENABLE_PROGRAM "Build the program fhhos4 and its tests. OFF: the library fhhos4_AMG alone, which only needs Eigen" ON)
option(ENABLE_GMSH "Build with the GMSH interface" ON)
option(ENABLE_AGMG "Build with AGMG" OFF)
option(ENABLE_CGAL "Build with CGAL" ON)
option(SMALL_INDEX "Use type 'int' for the indexes" ON)
option(ENABLE_1D "Enable 1D" OFF)
option(ENABLE_2D "Enable 2D" ON)
option(ENABLE_3D "Enable 3D" ON)
option(ENABLE_DG "Build the DG programs (-discr dg)" OFF)
option(ENABLE_FEM "Build the FEM programs (-discr fem)" OFF)
option(ENABLE_TESTS "Build the test suite (fhhos4_tests). Turn off to save its compilation time during development." ON)
option(ENABLE_OPENMP "Use OpenMP for the parallel loops (and Eigen's internal parallelism)" ON)
option(ENABLE_NATIVE_ARCH "Optimize for the CPU of the build machine (-march=native: AVX2, FMA...). The binary may not run on other CPUs." OFF)

set(EXECUTABLE_OUTPUT_PATH bin/)
include(GNUInstallDirs)

# With Ninja, limit the number of parallel compilations (ninja runs nproc+2 by default): the
# compilation of Program.cpp takes ~4 GB of RAM, and the build time is set by this single file
# anyway, so more parallel jobs don't make the build faster.
set(MAX_COMPILE_JOBS 2 CACHE STRING "Maximum number of parallel compilations with the Ninja generator")
set_property(GLOBAL PROPERTY JOB_POOLS compile_pool=${MAX_COMPILE_JOBS})
set(CMAKE_JOB_POOL_COMPILE compile_pool)

# Pass -DCMAKE_VERBOSE_MAKEFILE=ON (or `cmake --build build -- VERBOSE=1`) to see full compiler commands

#----------------------------------------------------------------------------------------------------------------------
# Library fhhos4_AMG
#----------------------------------------------------------------------------------------------------------------------
# The AMG of fhhos4 for hybrid discretizations (U-AMG), for the codes that discretize the problem themselves
# (library/README.md). Shared, and compiled with -fvisibility=hidden: only the public API (library/include/fhhos4) is
# exported, fhhos4's internals and its copy of Eigen's functions stay private. Its own compile definitions and link
# libraries: the calling code only needs Eigen.
add_library(fhhos4_AMG SHARED library/Solver.cpp ${CMAKE_BINARY_DIR}/generated/Version.cpp)
add_library(fhhos4::Solver ALIAS fhhos4_AMG)
# SOVERSION: binary compatibility within a minor version only, since adding a parameter (a public member) to
# fhhos4::Solver changes its size.
set_target_properties(fhhos4_AMG PROPERTIES
	EXPORT_NAME AMG
	CXX_VISIBILITY_PRESET hidden
	VISIBILITY_INLINES_HIDDEN ON
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
target_include_directories(fhhos4_AMG
	PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/library/include> $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Version, used by the library and the program: the public header fhhos4_version.h, and the definition of
# fhhos4::Version(), which adds the commit (git describe). The commit changes without CMake running again: Version.cpp
# is regenerated at each build, but rewritten only when it changes, so a new commit recompiles this small file alone.
configure_file(library/fhhos4_version.h.in ${CMAKE_BINARY_DIR}/generated/include/fhhos4_version.h @ONLY)
find_package(Git QUIET)
set(git_version_command ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR} -DVERSION=${PROJECT_VERSION}
	-DGIT_EXECUTABLE=${GIT_EXECUTABLE} -DOUTPUT=${CMAKE_BINARY_DIR}/generated/Version.cpp
	-P ${CMAKE_SOURCE_DIR}/cmake/GitVersion.cmake)
execute_process(COMMAND ${git_version_command}) # Version.cpp exists after the configuration
add_custom_target(fhhos4_git_version
	COMMAND ${git_version_command}
	BYPRODUCTS ${CMAKE_BINARY_DIR}/generated/Version.cpp
	COMMENT "Version: git describe"
	VERBATIM)
add_dependencies(fhhos4_AMG fhhos4_git_version)

if(SMALL_INDEX)
	target_compile_definitions(fhhos4_AMG PRIVATE SMALL_INDEX)
	message(STATUS "SMALL_INDEX: ON (int)")
else()
	message(STATUS "SMALL_INDEX: OFF (size_t)")
endif()

# Include Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(fhhos4_AMG PUBLIC Eigen3::Eigen)

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(fhhos4_AMG PRIVATE Threads::Threads)

# Include OpenMP
if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
	target_link_libraries(fhhos4_AMG PRIVATE OpenMP::OpenMP_CXX)
	message(STATUS "OpenMP enabled")
else()
	message(STATUS "OpenMP disabled: sequential execution.")
endif()

if(ENABLE_NATIVE_ARCH)
	target_compile_options(fhhos4_AMG PRIVATE -march=native)
	message(STATUS "ENABLE_NATIVE_ARCH: ON (-march=native)")
endif()

# Installation (component library), and CMake package of the library: find_package(fhhos4) then
# target_link_libraries(<target> fhhos4::AMG), from the installation (cmake --install) or from this build tree
# (-Dfhhos4_DIR=<this build directory>). find_package(fhhos4 1.2) accepts the versions 1.x with x >= 2: the source
# compatibility of the API (the binary compatibility is the SOVERSION's).
include(CMakePackageConfigHelpers)
install(TARGETS fhhos4_AMG EXPORT fhhos4Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT library)
install(FILES library/include/fhhos4 library/include/fhhos4_error.h
	${CMAKE_BINARY_DIR}/generated/include/fhhos4_version.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT library)
install(EXPORT fhhos4Targets NAMESPACE fhhos4:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4 COMPONENT library)
configure_package_config_file(library/fhhos4Config.cmake.in ${CMAKE_BINARY_DIR}/fhhos4Config.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4)
write_basic_package_version_file(${CMAKE_BINARY_DIR}/fhhos4ConfigVersion.cmake COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_BINARY_DIR}/fhhos4Config.cmake ${CMAKE_BINARY_DIR}/fhhos4ConfigVersion.cmake
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4 COMPONENT library)
export(EXPORT fhhos4Targets NAMESPACE fhhos4:: FILE ${CMAKE_BINARY_DIR}/fhhos4Targets.cmake)

#----------------------------------------------------------------------------------------------------------------------
# Program fhhos4 and its tests
#----------------------------------------------------------------------------------------------------------------------
if(NOT ENABLE_PROGRAM)
	message(STATUS "ENABLE_PROGRAM: OFF (the library fhhos4_AMG alone)")
	return()
endif()

# Files to compile
file(
	GLOB_RECURSE
	source_files
	CONFIGURE_DEPENDS
	src/*
)
list(REMOVE_ITEM source_files ${CMAKE_SOURCE_DIR}/src/main.cpp)

# Include AGMG compiled files
if(ENABLE_AGMG)
	message(STATUS "AGMG enabled")
	set(OBJS ${AGMG_DIR}/dagmg.o ${AGMG_DIR}/dagmg_mumps.o)
else()
	message(STATUS "AGMG disabled. Add -DENABLE_AGMG=ON and -DAGMG_DIR=<path/to/.o/files> to enable it.")
endif()

# Core library: everything except main.cpp, shared by the fhhos4 executable and the test executable.
# Compile definitions/link libraries are set PUBLIC so both consumers see the exact same
# preprocessor configuration (critical since most of the code is header-only and gets
# recompiled per translation unit).
add_library(fhhos4_core STATIC ${OBJS} ${source_files})
target_include_directories(fhhos4_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(fhhos4_core PUBLIC fhhos4_AMG) # -s libuamg, and fhhos4::Error

# Executable file
add_executable(fhhos4 src/main.cpp)
target_link_libraries(fhhos4 PRIVATE fhhos4_core)

# Data (the meshes), found from the directory of the executable (FileSystem::DataDirectory()): share/fhhos4 in the
# installation, and in the build tree, where share/fhhos4/meshes links to data/meshes.
file(RELATIVE_PATH data_dir_from_bin ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_FULL_DATADIR}/fhhos4)
target_compile_definitions(fhhos4_core PUBLIC DATA_DIR_FROM_BIN="${data_dir_from_bin}")
cmake_path(ABSOLUTE_PATH data_dir_from_bin BASE_DIRECTORY ${CMAKE_BINARY_DIR}/bin NORMALIZE OUTPUT_VARIABLE build_data_dir)
file(MAKE_DIRECTORY ${build_data_dir})
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/data/meshes ${build_data_dir}/meshes SYMBOLIC)

# Installation (component program): the executable, which finds the library in the installation's lib directory, and
# the meshes
file(RELATIVE_PATH lib_dir_from_bin ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_FULL_LIBDIR})
set_target_properties(fhhos4 PROPERTIES INSTALL_RPATH "$ORIGIN/${lib_dir_from_bin}" INSTALL_RPATH_USE_LINK_PATH ON)
install(TARGETS fhhos4 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT program)
install(DIRECTORY data/meshes DESTINATION ${CMAKE_INSTALL_DATADIR}/fhhos4 COMPONENT program)

if(SMALL_INDEX)
	target_compile_definitions(fhhos4_core PUBLIC SMALL_INDEX)
endif()

if(ENABLE_1D)
	target_compile_definitions(fhhos4_core PUBLIC ENABLE_1D)
	message(STATUS "ENABLE_1D: ON")
else()
	message(STATUS "ENABLE_1D: OFF")
endif()

if(ENABLE_2D)
	target_compile_definitions(fhhos4_core PUBLIC ENABLE_2D)
	message(STATUS "ENABLE_2D: ON")
else()
	message(STATUS "ENABLE_2D: OFF")
endif()

if(ENABLE_3D)
	target_compile_definitions(fhhos4_core PUBLIC ENABLE_3D)
	message(STATUS "ENABLE_3D: ON")
else()
	message(STATUS "ENABLE_3D: OFF")
endif()

if(ENABLE_DG)
	target_compile_definitions(fhhos4_core PUBLIC DG_ENABLED)
	message(STATUS "ENABLE_DG: ON")
else()
	message(STATUS "ENABLE_DG: OFF")
endif()

if(ENABLE_FEM)
	target_compile_definitions(fhhos4_core PUBLIC FEM_ENABLED)
	message(STATUS "ENABLE_FEM: ON")
else()
	message(STATUS "ENABLE_FEM: OFF")
endif()

if(ENABLE_AGMG)
	target_compile_definitions(fhhos4_core PUBLIC AGMG_ENABLED)
endif()

# Include Eigen
target_link_libraries(fhhos4_core PUBLIC Eigen3::Eigen)

# Include CGAL (only used for the 2D polygonal elements)
if(ENABLE_CGAL AND NOT ENABLE_2D)
	message(STATUS "CGAL disabled because 2D is disabled (CGAL is only used for the 2D polygonal elements)")
elseif(ENABLE_CGAL)
	message(STATUS "CGAL enabled")
	# CGAL's own exported CMake files call cmake_minimum_required() internally, which resets
	# CMP0167 right before their own find_package(Boost) call and triggers a policy warning.
	# A plain cmake_policy(SET ...) gets clobbered by that nested call, so pin the default here
	# instead - it's rechecked on every cmake_minimum_required()/cmake_policy(VERSION) call.
	if(POLICY CMP0167)
		set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
	endif()
	find_package(CGAL REQUIRED)
	target_link_libraries(fhhos4_core PUBLIC CGAL::CGAL)
	target_compile_definitions(fhhos4_core PUBLIC CGAL_ENABLED)
else()
	message(STATUS "CGAL disabled. Add -DENABLE_CGAL=ON to enable it.")
endif()

# Include GMSH
if(ENABLE_GMSH)
	message(STATUS "GMSH enabled")
	if(GMSH_API)
		target_include_directories(fhhos4_core PUBLIC ${GMSH_API})
	else()
		find_path(GMSH_API NAMES gmsh.h)
		if(GMSH_API)
			target_include_directories(fhhos4_core PUBLIC ${GMSH_API})
			message(STATUS "GMSH header found: ${GMSH_API}/gmsh.h")
		else()
			message(WARNING "Option -DGMSH_API missing, and gmsh.h could not be found automatically. Set this option if gmsh.h is not found during compilation")
		endif()
	endif()
	if(GMSH_LIB)
		target_link_libraries(fhhos4_core PUBLIC "${GMSH_LIB}") # -lgmsh
	else()
		find_library(GMSH_LIB NAMES gmsh)
		if(GMSH_LIB)
			target_link_libraries(fhhos4_core PUBLIC "${GMSH_LIB}") # -lgmsh
			message(STATUS "GMSH library found: ${GMSH_LIB}")
		else()
			message(FATAL_ERROR "The GMSH dynamic library was not found. Note that GMSH must be compiled with the option -DENABLE_BUILD_DYNAMIC=1. Pass the library path with the option -DGMSH_LIB=<path>/libgmsh.so.")
		endif()
	endif()
	target_compile_definitions(fhhos4_core PUBLIC GMSH_ENABLED)
else()
	message(STATUS "GMSH disabled. Add -DENABLE_GMSH=ON to enable it.")
endif()

# Include pthread
target_link_libraries(fhhos4_core PUBLIC Threads::Threads)

# Include OpenMP
if(ENABLE_OPENMP)
	target_link_libraries(fhhos4_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(ENABLE_NATIVE_ARCH)
	target_compile_options(fhhos4_core PUBLIC -march=native)
endif()

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

	target_link_libraries(fhhos4_core PUBLIC "-lgfortran" "-lm" ${TPL_LAPACK_LIB} ${TPL_BLAS_LIB})
endif()

# Validation test suite (GoogleTest)
if(ENABLE_TESTS)
	enable_testing()
	find_package(GTest REQUIRED)
	include(GoogleTest)

	file(
		GLOB_RECURSE
		test_source_files
		CONFIGURE_DEPENDS
		tests/*.cpp
	)
	add_executable(fhhos4_tests ${test_source_files})
	target_include_directories(fhhos4_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
	target_link_libraries(fhhos4_tests PRIVATE fhhos4_core GTest::gtest_main)
	gtest_discover_tests(fhhos4_tests)
	message(STATUS "ENABLE_TESTS: ON")
else()
	message(STATUS "ENABLE_TESTS: OFF")
endif()
