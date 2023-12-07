cmake_minimum_required(VERSION 3.11)
project(
	boost-multi-adaptor-cuda
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 --extended-lambda --expt-relaxed-constexpr -Xcudafe \"--display_error_number\"")

	enable_testing()
	list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.1

	add_subdirectory(test)
endif()
