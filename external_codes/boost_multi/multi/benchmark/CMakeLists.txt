cmake_minimum_required(VERSION 3.10)

project(boost-multi-benchmark VERSION 0.1 LANGUAGES CXX CUDA)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)
find_package(benchmark REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_61 -std=c++17 --use_fast_math --expt-relaxed-constexpr --extended-lambda") #set(CMAKE_CUDA_STANDARD 17)

enable_testing()
include(CTest)

#file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_executable(thrust_assignment.cu.x         thrust_assignment.cu  )
add_test(NAME  thrust_assignment.cu.x COMMAND thrust_assignment.cu.x)
target_link_libraries(thrust_assignment.cu.x  PRIVATE benchmark::benchmark)

