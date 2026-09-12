cmake_minimum_required(VERSION 3.17)  # for CUDAToolkit; the cuda CI job runs 3.22

project(boost-multi-adaptors-lapack-test-cuda LANGUAGES CXX CUDA)

# GPU (cuSOLVER) tests for the lapack adaptor; unlike the CPU tests, these need neither
# a CPU LAPACK nor ENABLE_CIRCLE, so they are wired directly from the root under ENABLE_CUDA

find_package(Boost REQUIRED COMPONENTS headers)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
	set(CMAKE_CUDA_ARCHITECTURES native)  # local-build convenience; "native" needs cmake 3.24, CI always passes -DCMAKE_CUDA_ARCHITECTURES=<n>
endif()

find_package(CUDAToolkit REQUIRED)

enable_testing()
include(CTest)

add_executable(lapack_cusolver.cu.x cusolver.cu)
target_link_libraries(lapack_cusolver.cu.x PRIVATE multi CUDA::cusolver Boost::headers)

add_test(NAME lapack_cusolver.cu.x COMMAND $<TARGET_FILE:lapack_cusolver.cu.x>)
