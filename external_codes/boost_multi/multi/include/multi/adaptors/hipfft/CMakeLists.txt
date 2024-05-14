cmake_minimum_required(VERSION 3.11)

project(
	boost-multi-adaptors-hipfft
	VERSION 0.1
	LANGUAGES CXX
)

enable_language(HIP)
find_package(hipFFT REQUIRED)
# find_package(rocthrust REQUIRED)

# if(ENABLE_CUDA OR DEFINED CXXCUDA)
#   enable_language(CUDA)
#   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda")
# endif()

# find_package(CUDA QUIET)
#include_directories(${CUDA_INCLUDE_DIRS})
#link_libraries(${CUDA_CUFFT_LIBRARIES})

# include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(test)
