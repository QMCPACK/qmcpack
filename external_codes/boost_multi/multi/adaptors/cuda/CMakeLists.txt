cmake_minimum_required(VERSION 3.11)
project(boost-multi-adaptor-cuda VERSION 0.1 LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_language(CUDA)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++17 --extended-lambda --expt-relaxed-constexpr -Xcudafe \"--display_error_number --set_error_limit=2 --verbose_diagnostics --emit_warnings_as_errors --diag_suppress=implicit_return_from_non_void_function\"")

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

enable_testing()
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.1

add_subdirectory(thrust/test)
add_subdirectory(cublas/test)

