cmake_minimum_required(VERSION 3.11)
project(
	boost-multi-adaptor-cuda
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED)  # REQUIRED COMPONENTS unit_test_framework)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)

	enable_testing()

	add_subdirectory(test)
endif()
