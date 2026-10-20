cmake_minimum_required(VERSION 3.20)

# The version of fhhos4, written here only: fhhos4_version.h, fhhos4::Version() and the CMake package are generated from
# it. A release changes it, then gets the annotated tag v<version>.
project(fhhos4 VERSION 1.4.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# program options
option(ENABLE_PROGRAM "Build the program fhhos4 and its tests. OFF: the library alone, which only needs Eigen" ON)
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
option(ENABLE_PYTHON "Build the Python module of the library (fhhos4, with nanobind) and its tests" OFF)
option(ENABLE_MPI "MPI parallelism (distributed memory) of the library and the program, on several processes (mpiexec). OFF: a single process" OFF)

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
# Library fhhos4
#----------------------------------------------------------------------------------------------------------------------
# The AMG of fhhos4 for hybrid discretizations (U-AMG), for the codes that discretize the problem themselves
# (README_library.md). Shared, and compiled with -fvisibility=hidden: only the public API (library/include/fhhos4.h) is
# exported, fhhos4's internals and its copy of Eigen's functions stay private. Its own compile definitions and link
# libraries: the calling code only needs Eigen.
add_library(fhhos4_lib SHARED library/Solver.cpp ${CMAKE_BINARY_DIR}/generated/Version.cpp)
add_library(fhhos4::fhhos4 ALIAS fhhos4_lib) # fhhos4_lib: the name of the target in this project (fhhos4 is the program); libfhhos4.so, fhhos4::fhhos4 outside
# SOVERSION: binary compatibility within a minor version only, since adding a parameter (a public member) to
# fhhos4::Solver changes its size.
set_target_properties(fhhos4_lib PROPERTIES
	OUTPUT_NAME fhhos4
	EXPORT_NAME fhhos4
	CXX_VISIBILITY_PRESET hidden
	VISIBILITY_INLINES_HIDDEN ON
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
	VERSION ${PROJECT_VERSION}
	SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
target_include_directories(fhhos4_lib
	PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/library/include> $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/include>
		$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	PRIVATE ${CMAKE_SOURCE_DIR}/src)

# Version, used by the library and the program: the public header fhhos4_version.h, and the definition of
# fhhos4::Version(), which adds the commit (git describe). The commit changes without CMake running again: Version.cpp
# is regenerated at each build, but rewritten only when it changes, so a new commit recompiles this small file alone.
# fhhos4_version.h also gives the index type of the sparse matrices (SMALL_INDEX, src/Utils/Types.h): fhhos4::SparseMatrix.
if(SMALL_INDEX)
	set(FHHOS4_SPARSE_INDEX int)
else()
	set(FHHOS4_SPARSE_INDEX intmax_t)
endif()
# fhhos4_version.h also says whether the library is compiled with MPI (FHHOS4_MPI): its public header then includes mpi.h
if(ENABLE_MPI)
	set(FHHOS4_MPI 1)
else()
	set(FHHOS4_MPI 0)
endif()
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
add_dependencies(fhhos4_lib fhhos4_git_version)

if(SMALL_INDEX)
	target_compile_definitions(fhhos4_lib PRIVATE SMALL_INDEX)
	message(STATUS "SMALL_INDEX: ON (int)")
else()
	message(STATUS "SMALL_INDEX: OFF (intmax_t)")
endif()

# Include Eigen
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(fhhos4_lib PUBLIC Eigen3::Eigen)

# Include MPI (optional): public, since the public header includes mpi.h (fhhos4::Solver::Comm)
if(ENABLE_MPI)
	find_package(MPI REQUIRED COMPONENTS CXX)
	target_link_libraries(fhhos4_lib PUBLIC MPI::MPI_CXX)
	target_compile_definitions(fhhos4_lib PRIVATE MPI_ENABLED)
	message(STATUS "MPI enabled (MPI ${MPI_CXX_VERSION})")
else()
	message(STATUS "MPI disabled: a single process. Add -DENABLE_MPI=ON to enable it.")
endif()

# Include pthread
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(fhhos4_lib PRIVATE Threads::Threads)

# Include OpenMP
if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
	target_link_libraries(fhhos4_lib PRIVATE OpenMP::OpenMP_CXX)
	message(STATUS "OpenMP enabled")
else()
	message(STATUS "OpenMP disabled: sequential execution.")
endif()

if(ENABLE_NATIVE_ARCH)
	target_compile_options(fhhos4_lib PRIVATE -march=native)
	message(STATUS "ENABLE_NATIVE_ARCH: ON (-march=native)")
endif()

# Installation (component library), and CMake package of the library: find_package(fhhos4) then
# target_link_libraries(<target> fhhos4::fhhos4), from the installation (cmake --install) or from this build tree
# (-Dfhhos4_DIR=<this build directory>). find_package(fhhos4 1.2) accepts the versions 1.x with x >= 2: the source
# compatibility of the API (the binary compatibility is the SOVERSION's).
include(CMakePackageConfigHelpers)
install(TARGETS fhhos4_lib EXPORT fhhos4Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT library)
install(FILES library/include/fhhos4 library/include/fhhos4.h library/include/fhhos4_error.h
	${CMAKE_BINARY_DIR}/generated/include/fhhos4_version.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT library)
install(EXPORT fhhos4Targets NAMESPACE fhhos4:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4 COMPONENT library)
configure_package_config_file(library/fhhos4Config.cmake.in ${CMAKE_BINARY_DIR}/fhhos4Config.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4)
write_basic_package_version_file(${CMAKE_BINARY_DIR}/fhhos4ConfigVersion.cmake COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_BINARY_DIR}/fhhos4Config.cmake ${CMAKE_BINARY_DIR}/fhhos4ConfigVersion.cmake
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fhhos4 COMPONENT library)
export(EXPORT fhhos4Targets NAMESPACE fhhos4:: FILE ${CMAKE_BINARY_DIR}/fhhos4Targets.cmake)

# Python module of the library (library/python): the module fhhos4 in <build>/python, linked to libfhhos4.so (build
# tree only for now: not installed). Its test (ctest -R Python) runs with that directory in PYTHONPATH.
if(ENABLE_PYTHON)
	if(ENABLE_MPI)
		message(FATAL_ERROR "ENABLE_PYTHON: the Python module of the library compiled with MPI (ENABLE_MPI) is not available yet.")
	endif()
	find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
	execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
		OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
	find_package(nanobind CONFIG REQUIRED)
	nanobind_add_module(fhhos4_python NB_STATIC library/python/fhhos4_python.cpp)
	set_target_properties(fhhos4_python PROPERTIES OUTPUT_NAME fhhos4 LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python)
	target_link_libraries(fhhos4_python PRIVATE fhhos4_lib)
	message(STATUS "ENABLE_PYTHON: ON (module fhhos4 for ${Python_EXECUTABLE})")

	enable_testing()
	add_test(NAME PythonModule COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/library/python/tests/test_module.py)
	set_tests_properties(PythonModule PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python)
endif()

# The tests of the library on several processes (tests/MPI, their own main()), with the library's code only (src/
# headers and libfhhos4.so, not the program): one ctest test per number of processes, each run by mpiexec, with 2
# threads per process
if(ENABLE_MPI AND ENABLE_TESTS)
	enable_testing()
	find_package(GTest REQUIRED)
	file(GLOB mpi_test_source_files CONFIGURE_DEPENDS tests/MPI/*.cpp)
	add_executable(fhhos4_mpi_tests ${mpi_test_source_files})
	target_include_directories(fhhos4_mpi_tests PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/tests)
	target_compile_definitions(fhhos4_mpi_tests PRIVATE MPI_ENABLED)
	if(SMALL_INDEX)
		target_compile_definitions(fhhos4_mpi_tests PRIVATE SMALL_INDEX)
	endif()
	target_link_libraries(fhhos4_mpi_tests PRIVATE fhhos4_lib GTest::gtest Threads::Threads)
	if(ENABLE_OPENMP)
		target_link_libraries(fhhos4_mpi_tests PRIVATE OpenMP::OpenMP_CXX)
	endif()
	if(ENABLE_NATIVE_ARCH)
		target_compile_options(fhhos4_mpi_tests PRIVATE -march=native)
	endif()
	foreach(np 1 2 3 4)
		add_test(NAME MPI_np${np} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:fhhos4_mpi_tests> ${MPIEXEC_POSTFLAGS})
		set_tests_properties(MPI_np${np} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 600 PROCESSORS ${np})
	endforeach()
endif()

#----------------------------------------------------------------------------------------------------------------------
# Program fhhos4 and its tests
#----------------------------------------------------------------------------------------------------------------------
if(NOT ENABLE_PROGRAM)
	message(STATUS "ENABLE_PROGRAM: OFF (the library alone)")
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
target_link_libraries(fhhos4_core PUBLIC fhhos4_lib) # -s libuamg, and fhhos4::Error

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

# MPI: linked through fhhos4_lib
if(ENABLE_MPI)
	target_compile_definitions(fhhos4_core PUBLIC MPI_ENABLED)
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
	list(FILTER test_source_files EXCLUDE REGEX "/tests/MPI/")
	if(ENABLE_MPI)
		# The main() of the MPI tests: MPI initialized (one process, without mpiexec)
		list(APPEND test_source_files ${CMAKE_SOURCE_DIR}/tests/MPI/main.cpp)
	endif()
	add_executable(fhhos4_tests ${test_source_files})
	target_include_directories(fhhos4_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
	if(ENABLE_MPI)
		target_link_libraries(fhhos4_tests PRIVATE fhhos4_core GTest::gtest)
	else()
		target_link_libraries(fhhos4_tests PRIVATE fhhos4_core GTest::gtest_main)
	endif()
	gtest_discover_tests(fhhos4_tests)
	message(STATUS "ENABLE_TESTS: ON")

	# The program on several processes (MPI): GMSH meshes, U-AMG through the library (uncondensed system) and the
	# program's (its condensed system), the L2 errors of one process. -no-cache: no concurrent writes of GMSH's cache by
	# the tests run in parallel.
	if(ENABLE_MPI)
		set(mpi_program_cases
			"2D_k1_fcglibuamg|3|-geo square -mesh tri -mesher gmsh -n 32 -k 1 -s fcglibuamg|L2 Error = 4.48e-04"
			"2D_k2_fcguamg|4|-geo square -mesh tri -mesher gmsh -n 32 -k 2 -s fcguamg -smoothers hbgs,hrbgs|L2 Error = 1.23e-05"
			"3D_k0_fcguamg|2|-geo cube -mesh tetra -n 8 -k 0 -s fcguamg|L2 Error = 7.18e-01"
			"3D_k2_fcglibuamg|4|-geo cube -mesh tetra -n 4 -k 2 -s fcglibuamg|L2 Error = 9.83e-02"
			"InHouseMeshRefused|2|-geo square -mesh tri -mesher inhouse -n 8 -k 0 -s fcglibuamg|the meshes are GMSH's only")
		foreach(case ${mpi_program_cases})
			string(REPLACE "|" ";" fields "${case}")
			list(GET fields 0 name)
			list(GET fields 1 np)
			list(GET fields 2 arguments)
			list(GET fields 3 expected)
			separate_arguments(arguments)
			add_test(NAME ProgramMPI_${name}_np${np} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:fhhos4> ${arguments} -no-cache ${MPIEXEC_POSTFLAGS}
				WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
			set_tests_properties(ProgramMPI_${name}_np${np} PROPERTIES PASS_REGULAR_EXPRESSION "${expected}" ENVIRONMENT OMP_NUM_THREADS=2 TIMEOUT 600 PROCESSORS ${np})
		endforeach()
	endif()

else()
	message(STATUS "ENABLE_TESTS: OFF")
endif()
