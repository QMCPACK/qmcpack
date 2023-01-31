cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

#project(
#	boost-multi-adaptors-cuda-cublas-test
#	VERSION 0.1
#	LANGUAGES CXX CUDA
#)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

find_package(BLAS REQUIRED)
find_path(
	BLAS_INCLUDE_DIRS
	cblas.h
	/usr/include
	/usr/local/include
	$ENV{BLAS_HOME}/include
)

link_libraries(${BLAS_LIBRARIES})
include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number\"")
endif()

find_package(CUDA)

enable_testing()
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.17
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS # herk.cu
	gemm.cu
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()
	# target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)
	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
	target_include_directories(${TEST_EXE}        PRIVATE ${CUDA_INCLUDE_DIRS})

	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})

	target_link_libraries(${TEST_EXE} PRIVATE ${CUDA_LIBRARIES})
	target_include_directories(${TEST_EXE} PRIVATE /opt/nvidia/hpc_sdk/Linux_x86_64/22.3/math_libs/include)

	# if(NOT ENABLE_CUDA) target_compile_options  (${TEST_EXE} PRIVATE $<$<CXX_COMPILER_ID:GNU>: -Werror -Wall -Wextra -fno-common -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
	# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
	# -Werror -Wall -Wextra -fno-common -Wpedantic -Wmove> $<$<CXX_COMPILER_ID:Intel>: -Werror -Wall -Wextra -fno-common -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized
	# -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
	# -Werror -diag-error:3846 > $<$<CXX_COMPILER_ID:MSVC>: /W4>) endif()
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
