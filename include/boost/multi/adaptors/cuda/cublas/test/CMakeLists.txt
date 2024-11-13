cmake_minimum_required(VERSION 3.17)  # for CUDAToolkit

find_package(Boost REQUIRED COMPONENTS)  # unit_test_framework)

if((NOT
	CMAKE_CXX_COMPILER_ID
	STREQUAL
	"PGI"
   )
   AND (NOT
	   CMAKE_CXX_COMPILER_ID
		STREQUAL
        "NVHPC"
	)
   AND (NOT
		DART_COMPILER_NAME
		STREQUAL
		"nvcc"
	)
   AND (NOT
		DART_COMPILER_NAME
		STREQUAL
		"icpc"
	   )
)
	# find_package(Boost REQUIRED COMPONENTS unit_test_framework)
	# link_libraries("-lboost_unit_test_framework")
	find_package(BLAS REQUIRED)
	find_path(
		BLAS_INCLUDE_DIRS
		cblas.h
		/usr/include
		/usr/local/include
		$ENV{BLAS_HOME}/include
	)
#   include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})
	link_libraries(${BLAS_LIBRARIES})
else()
	link_libraries("-lblas") # cmake cannot detect BLAS with pgi/nvc++ but it ships with its own version
endif()

# include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})

if(ENABLE_HIP)
	enable_language(HIP)
	find_package(hipBLAS REQUIRED)
	find_package(rocthrust REQUIRED)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
	find_package(CUDAToolkit REQUIRED COMPONENTS cuBLAS)  # requires cmake 3.17
endif()

enable_testing()
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

if(ENABLE_HIP)
	add_executable(all.hip.x all.hip)
	#set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE HIP)
	target_link_libraries(all.hip.x PRIVATE multi roc::hipblas Boost::boost)
	# target_compile_definitions(all.hip.x PRIVATE BOOST_TEST_DYN_LINK=1)

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
	nrm2.cu
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

		# target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
		target_include_directories(${TEST_EXE}        PRIVATE ${CUDA_INCLUDE_DIRS})
	
		#target_compile_options(${TEST_EXE} PRIVATE --expt-relaxed-constexpr)

		# target_link_libraries(${TEST_EXE} PRIVATE Boost::unit_test_framework)
	
		target_link_libraries(${TEST_EXE} PRIVATE ${CUDA_LIBRARIES})
		target_link_libraries(${TEST_EXE} PRIVATE CUDA::cublas)
		target_link_libraries(${TEST_EXE} PRIVATE multi)  # Boost::unit_test_framework)

		# target_include_directories(${TEST_EXE} PRIVATE /opt/nvidia/hpc_sdk/Linux_x86_64/22.3/math_libs/include)
	
		# if(NOT ENABLE_CUDA) target_compile_options  (${TEST_EXE} PRIVATE $<$<CXX_COMPILER_ID:GNU>: -Werror -Wall -Wextra -fno-common -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
		# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
		# -Werror -Wall -Wextra -fno-common -Wpedantic -Wmove> $<$<CXX_COMPILER_ID:Intel>: -Werror -Wall -Wextra -fno-common -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized
		# -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
		# -Werror -diag-error:3846 > $<$<CXX_COMPILER_ID:MSVC>: /W4>) endif()
	
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		# target_compile_options(${TEST_EXE} PRIVATE -extended-lambda)
		
		add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
	endif()
endforeach()
endif()
