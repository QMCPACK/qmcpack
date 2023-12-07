# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.16)

#project(
# boost-multi-adaptors-blas
# VERSION 0.1
# LANGUAGES CXX
#)

message("CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")

set(BLA_VENDOR Intel10_64lp)
find_package(BLAS)
if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work
	message("Multi/BLAS: MKL environment detected")
	add_definitions(-DBLAS_DOT_RETURNS_VOID)
	# add_definitions(-DFORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID)
	SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64")
else()
	message("Multi/BLAS: MKL environment not detected, looking for other BLAS")
	set(BLA_VENDOR Apple)
	if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work
		message("Multi/BLAS: Applen BLAS (Accelerate) environment detected")
		add_definitions(-DBLAS_DOT_RETURNS_VOID)
	else()
		message("Multi/BLAS: Apple BLAS (Accelerate) environment NOT detected")
		if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
			set(BLA_VENDOR NVHPC)
			find_package(BLAS)
			if(BLAS_FOUND)
				message("Multi/BLAS: NVHPC BLAS environment detected")
				add_definitions(-DBLAS_DOT_RETURNS_VOID)
				# add_definitions(-DFORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID)
				# add_definitions(-DRETURN_BY_STACK)  # TODO(correaa) check
			endif()
		else()
			message("Multi/BLAS: NVHPC environment not detected, looking for other BLAS")
			unset(BLA_VENDOR)
			find_package(BLAS)
			find_path(
				BLAS_INCLUDE_DIRS
				cblas.h
				/usr/include
				/usr/local/include
				$ENV{BLAS_HOME}/include
			)
			if(BLAS_FOUND)
				message("Some BLAS found")
			endif()
		endif()
	endif()
endif()

if(BLAS_FOUND)
	message("Some BLAS found: linker flags ${BLAS_LINKER_FLAGS}, libs ${BLAS_LIBRARIES}, libs95 ${BLAS95_LIBRARIES}")
endif()

link_libraries(${BLAS_LIBRARIES})
# link_libraries(BLAS::BLAS) # needs cmake 3.18

# include_directories(${CMAKE_BINARY_DIR})

if(BLAS_FOUND)
  add_subdirectory(test)
else()
  message(WARNING "BLAS not found, BLAS-adaptor tests will not be compiled and run. If you want this feature install BLAS, for example please run:"
                "\n    sudo apt install libblas-dev"
                "\n    sudo dnf install blas-devel  # in Fedora")
endif()
