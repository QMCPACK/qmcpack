cmake_minimum_required(VERSION 3.18)  # 3.18 for LAPACK::LAPACK

project(
	boost-multi-adaptor-lapack
	VERSION 0.1
	LANGUAGES CXX
)

add_library(multi-lapack INTERFACE)

if(NOT DEFINED ENABLE_CIRCLE)
	return()
endif()

find_package(LAPACK REQUIRED)

foreach(lib ${LAPACK_LIBRARIES})
	# https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors
	if(${lib} MATCHES "mkl")
		message("Some LAPACK found matches MKL")
		target_compile_definitions(multi-lapack INTERFACE _MULTI_USING_LAPACK_MKL)
		# SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
		# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
		# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
		# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
		# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
	endif()
	if(${lib} MATCHES "hpc_sdk")
		message("Some LAPACK found matches hpc_sdk (nvhpc)")  # nvhpc may still use a different version of BLAS
		# add_definitions(-DBLAS_DOT_RETURNS_VOID)
		# target_compile_definitions(multi INTERFACE BLAS_DOT_RETURNS_VOID)
	endif()
	if(${lib} MATCHES "Accelerate")
		message("Some LAPACK found matches Accelerate (Apple) [beware of sdot and snrm2 bugs]")
		# message(WARNING "Apple Accelerate BLAS is known to have bugs in single precission function `sdot` and `smrm2`, be careful: https://stackoverflow.com/a/77017238/225186, https://fortran-lang.discourse.group/t/how-many-blas-libraries-have-this-error/4454/23, https://forums.developer.apple.com/forums/thread/717757")
		# add_definitions(-DBLAS_DOT_RETURNS_VOID)
		# target_compile_definitions(multi INTERFACE BLAS_DOT_RETURNS_VOID)
	endif()
endforeach()


target_link_libraries(multi-lapack INTERFACE multi LAPACK::LAPACK)

add_subdirectory(./test)
