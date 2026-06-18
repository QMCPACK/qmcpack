# CXX=g++-9 cmake .. -D CMAKE_BUILD_TYPE=Release -D CMAKE_CUDA_HOST_COMPILER=g++-9
cmake_minimum_required(VERSION 3.10)

project(
	boost-multi-benchmark
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework system iostreams)
find_package(benchmark REQUIRED)

set(BLA_VENDOR Intel10_64lp) # . /opt/intel/oneapi/mkl/latest/env/vars.sh
find_package(BLAS)
if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work
	message("Multi/BLAS: MKL environment detected")
	add_definitions(-DRETURN_BY_STACK)
else()
	message("Multi/BLAS: MKL environment not detected, looking for other BLAS")
	unset(BLA_VENDOR)
	find_package(BLAS REQUIRED)
	find_path(
		BLAS_INCLUDE_DIRS
		cblas.h
		/usr/include
		/usr/local/include
		$ENV{BLAS_HOME}/include
	)
endif()
link_libraries(${BLAS_LIBRARIES})

#find_package(BLAS REQUIRED)
#find_package(TBB REQUIRED HINTS "/usr")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_61 -std=c++17 --use_fast_math --expt-relaxed-constexpr --extended-lambda") # set(CMAKE_CUDA_STANDARD 17)

enable_testing()
include(CTest)

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#add_executable(thrust_assignment.cu.x thrust_assignment.cu)
#add_test(NAME thrust_assignment.cu.x COMMAND thrust_assignment.cu.x)
#target_link_libraries(thrust_assignment.cu.x PRIVATE benchmark::benchmark)

add_executable(matrix_multiplication.cpp.x matrix_multiplication.cpp)
target_link_libraries(matrix_multiplication.cpp.x benchmark::benchmark)
target_link_libraries(matrix_multiplication.cpp.x -ltbb)
add_test(NAME matrix_multiplication.cpp.x COMMAND matrix_multiplication.cpp.x)

add_executable(algorithms_gemm.cpp.x algorithms_gemm.cpp)
target_link_libraries(algorithms_gemm.cpp.x benchmark::benchmark Boost::unit_test_framework)
target_link_libraries(algorithms_gemm.cpp.x -ltbb)
add_test(NAME algorithms_gemm.cpp.x COMMAND algorithms_gemm.cpp.x)

add_executable(serialization serialization.cpp)
target_link_libraries(serialization benchmark::benchmark Boost::serialization Boost::iostreams)
add_test(NAME serialization COMMAND serialization)
