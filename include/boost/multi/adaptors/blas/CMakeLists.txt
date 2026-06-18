cmake_minimum_required(VERSION 3.16)

#project(
# boost-multi-adaptors-blas
# VERSION 0.1
# LANGUAGES CXX
#)

find_package(BLAS)

if(BLAS_FOUND)
	message("Some BLAS found: linker flags: ${BLAS_LINKER_FLAGS}, libs: ${BLAS_LIBRARIES}, libs95: ${BLAS95_LIBRARIES}")
	foreach(lib ${BLAS_LIBRARIES})
		# https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors
		if(${lib} MATCHES "mkl")
			message("Some BLAS found matches MKL")
			add_definitions(-DBLAS_DOT_RETURNS_VOID)
			target_compile_definitions(multi INTERFACE BLAS_DOT_RETURNS_VOID)
			# add_definitions(-D_MULTI_USING_BLAS_MKL)
			# in some systems with MKL, regular BLAS headers need to be found for it to work
			SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
			SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
			SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
		endif()
		if(${lib} MATCHES "hpc_sdk")
			message("Some BLAS found matches hpc_sdk (nvhpc)")  # nvhpc may still use a different version of BLAS
			add_definitions(-DBLAS_DOT_RETURNS_VOID)
			target_compile_definitions(multi INTERFACE BLAS_DOT_RETURNS_VOID)
		endif()
		if(${lib} MATCHES "Accelerate")
			message("Some BLAS found matches Accelerate (Apple) [beware of sdot and snrm2 bugs]")
			# message(WARNING "Apple Accelerate BLAS is known to have bugs in single precission function `sdot` and `smrm2`, be careful: https://stackoverflow.com/a/77017238/225186, https://fortran-lang.discourse.group/t/how-many-blas-libraries-have-this-error/4454/23, https://forums.developer.apple.com/forums/thread/717757")
			# add_definitions(-DBLAS_DOT_RETURNS_VOID)
			target_compile_definitions(multi INTERFACE BLAS_DOT_RETURNS_VOID)
		endif()
	endforeach()
endif()

if(BLAS_FOUND)
	find_package(Boost NO_MODULE)  # COMPONENTS boost)  # headers unit_test_framework)
	if(Boost_FOUND)
		add_subdirectory(test)
	else()
		message(WARNING "Boost Unit Test Framework not found, BLAS-adaptor tests will not be compiled and run. If you want this feature, install Boost")
	endif()
else()
	message(WARNING
		"BLAS not found, BLAS-adaptor tests will not be compiled and run. If you want this feature install BLAS, for example please run:"
		"\n    sudo apt install libblas-dev # in Ubuntu"
		"\n    sudo dnf install blas-devel  # in Fedora")
	message(WARNING "BLA_VENDOR was set to ${BLA_VENDOR}\n")
endif()
