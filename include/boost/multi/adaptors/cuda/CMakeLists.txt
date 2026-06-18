cmake_minimum_required(VERSION 3.11)
#project(
#   boost-multi-adaptor-cuda
#   VERSION 0.1
#   LANGUAGES CXX
#)

# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost REQUIRED COMPONENTS)  # unit_test_framework)

if(ENABLE_HIP)
  enable_language(HIP)  # may need `module load rocm`
  enable_testing()

  add_subdirectory(cublas/test)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 --extended-lambda --expt-relaxed-constexpr")

	enable_testing()
	list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.1

	add_subdirectory(cublas/test)
endif()
