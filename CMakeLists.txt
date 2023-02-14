# -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.16)

project(
	mpi3
	VERSION 0.79.0
	DESCRIPTION "B-MPI3 is a C++ library wrapper for version 3.1 of the MPI standard interface that simplifies the utilization and maintenance of MPI code."
	HOMEPAGE_URL "https://gitlab.com/correaa/boost-mpi3"
	LANGUAGES CXX
)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/mpi3/dummy>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
#set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

# this makes CM FetchContent friendly https://www.foonathan.net/2022/06/cmake-fetchcontent/
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  	return()
endif()

include(CTest)

enable_testing()

add_subdirectory(test)
