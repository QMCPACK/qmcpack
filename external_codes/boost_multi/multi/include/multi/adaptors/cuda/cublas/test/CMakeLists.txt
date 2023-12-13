cmake_minimum_required(VERSION 3.17)  # for CUDAToolkit

set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

find_package(BLAS)
if(BLAS_FOUND)
else()
  message(ERROR "BLAS not found. Please try:"
                "\n    sudo apt install libblas-dev"
                "\n    sudo dnf install blas-devel  # in Fedora")
endif()
find_path(
	BLAS_INCLUDE_DIRS
	cblas.h
	/usr/include
	/usr/local/include
	$ENV{BLAS_HOME}/include
)

link_libraries(${BLAS_LIBRARIES})
include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})

if(ENABLE_HIP)
	enable_language(HIP)
	find_package(hipBLAS REQUIRED)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -Xcudafe \"--display_error_number\"")
	find_package(CUDAToolkit REQUIRED COMPONENTS cuBLAS)  # requires cmake 3.17
endif()

enable_testing()
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

if(ENABLE_HIP)
	add_executable(all.hip.x all.hip)
	#set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE HIP)
	target_link_libraries(all.hip.x PRIVATE multi roc::hipblas Boost::unit_test_framework)
	add_test(NAME all.hip.x COMMAND ./all.hip.x)
endif()

if(ENABLE_CUDA)
# find_package(CUDA REQUIRED)

find_package(CUDAToolkit REQUIRED COMPONENTS cuBLAS)

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
	all.cu
	axpy.cu
	dot.cu
	gemm.cu
	gemv.cu
	# herk.cu
	scal.cu
	trsm.cu
)

foreach(TEST_FILE ${TEST_SRCS})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set(TEST_EXE "${TEST_FILE}.x")
		add_executable(${TEST_EXE} ${TEST_FILE})

		# target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)
		target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
		target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})

		target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
		target_include_directories(${TEST_EXE}        PRIVATE ${CUDA_INCLUDE_DIRS})
	
		target_link_libraries(${TEST_EXE} PRIVATE Boost::unit_test_framework)
	
		target_link_libraries(${TEST_EXE} PRIVATE ${CUDA_LIBRARIES})
		target_link_libraries(${TEST_EXE} PRIVATE CUDA::cublas)

		# target_include_directories(${TEST_EXE} PRIVATE /opt/nvidia/hpc_sdk/Linux_x86_64/22.3/math_libs/include)
	
		# if(NOT ENABLE_CUDA) target_compile_options  (${TEST_EXE} PRIVATE $<$<CXX_COMPILER_ID:GNU>: -Werror -Wall -Wextra -fno-common -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
		# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
		# -Werror -Wall -Wextra -fno-common -Wpedantic -Wmove> $<$<CXX_COMPILER_ID:Intel>: -Werror -Wall -Wextra -fno-common -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized
		# -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
		# -Werror -diag-error:3846 > $<$<CXX_COMPILER_ID:MSVC>: /W4>) endif()
	
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		target_compile_options(${TEST_EXE} PRIVATE -extended-lambda)
		
		add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
	endif()
endforeach()
endif()