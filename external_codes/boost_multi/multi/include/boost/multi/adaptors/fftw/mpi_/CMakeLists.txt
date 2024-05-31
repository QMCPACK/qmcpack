cmake_minimum_required(VERSION 3.17)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#project(
#	boost-multi-adaptors-fftw
#	VERSION 0.1
#	LANGUAGES CXX
#)

find_package(PkgConfig)
pkg_search_module(
	FFTW
	fftw3
	IMPORTED_TARGET
)

link_libraries(-lfftw3_mpi)

add_subdirectory(test)
