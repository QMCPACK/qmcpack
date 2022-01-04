cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(
	boost-multi-adaptors-fftw
	VERSION 0.1
	LANGUAGES CXX
)

find_package(PkgConfig REQUIRED)
pkg_search_module(
	FFTW
	REQUIRED
	fftw3
	IMPORTED_TARGET
)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

link_libraries(${BLAS_LIBRARIES})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(test)
