cmake_minimum_required(VERSION 3.10)

project(
	boost-multi-benchmark-cuda
	VERSION 0.1
	LANGUAGES CXX CUDA
)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)
find_package(benchmark REQUIRED)

set(BLA_VENDOR Intel10_64lp)
find_package(BLAS)
if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work
	message("Multi/BLAS: MKL environment detected")
	add_definitions(-DRETURN_BY_STACK)
else()
	message("Multi/BLAS: MKL environment not detected, looking for other BLAS")
	unset(BLA_VENDOR)
	find_package(BLAS REQUIRED)
endif()

link_libraries(${BLAS_LIBRARIES})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_61 -std=c++17 --use_fast_math --expt-relaxed-constexpr --extended-lambda") # set(CMAKE_CUDA_STANDARD 17)

enable_testing()
include(CTest)

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

include_directories(../../..)

add_executable(axpy.cu.x axpy.cu)
add_test(NAME axpy COMMAND axpy.cu.x)
target_link_libraries(axpy.cu.x PRIVATE benchmark::benchmark)
target_link_libraries(axpy.cu.x PRIVATE -lcublas)
