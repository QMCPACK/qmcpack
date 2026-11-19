cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(
	boost-multi-adaptors-cutensor
	VERSION 0.1
	LANGUAGES CXX
)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(cuTENSOR)

if (cuTENSOR_FOUND)
	message(STATUS "cuTENSOR found")
	add_subdirectory(test)
else()
	message(STATUS "cuTENSOR NOT found")
endif()
