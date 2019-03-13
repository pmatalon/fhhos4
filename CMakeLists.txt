cmake_minimum_required (VERSION 2.8)

project(dghho)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(TPL_ENABLE_AGMG FALSE CACHE BOOL "Whether AGMG is enabled")

# Include Eigen
include_directories(eigen-eigen-b3f3d4950030)

file(
	GLOB_RECURSE
	source_files
	src/*
)

add_executable(dghho ${source_files})

# Link to Third Party Libraries (TPL)
# AGMG (requires Blas and Lapack)
if(${TPL_ENABLE_AGMG})
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

	file(WRITE 
		src/Utils/agmg.h 
		"#define TPL_ENABLE_AGMG"
	)

	message(STATUS "AGMG enabled")
else()
	file(WRITE 
		src/Utils/agmg.h 
		"#undef TPL_ENABLE_AGMG"
	)
	message(STATUS "AGMG disabled")
endif()
