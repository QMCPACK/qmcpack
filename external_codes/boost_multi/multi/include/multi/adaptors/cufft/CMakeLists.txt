# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(
	boost-multi-adaptors-blas
	VERSION 0.1
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda")
endif()

find_package(CUDA QUIET)
include_directories(${CUDA_INCLUDE_DIRS})
link_libraries(${CUDA_CUFFT_LIBRARIES})

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(test)
