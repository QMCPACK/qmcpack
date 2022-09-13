# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
#[=[Multi Test suite can be run like this:
  mkdir -p build
  cd build
  cmake .. [-DENABLE_CUDA=1]
  make -j
  ctest -j --output-on-error [-T memcheck]
  exit
#]=]
cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

#project(
#	boost-multi-adaptors-blas-test
#	VERSION 0.1
#	LANGUAGES CXX
#)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if((NOT
	CMAKE_CXX_COMPILER_ID
	STREQUAL
	"PGI"
   )
   AND (NOT
		DART_COMPILER_NAME
		STREQUAL
		"nvcc"
	   )
)
	find_package(Boost REQUIRED COMPONENTS unit_test_framework)
	link_libraries("-lboost_unit_test_framework")

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
	find_package(Boost REQUIRED) # cmake cannot detect this component with pgi compiler
	link_libraries("-lboost_unit_test_framework")

	link_libraries("-lblas") # cmake cannot detect BLAS with pgi/nvc++ but it ships with its own version
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
	add_definitions(-DRETURN_BY_STACK)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
#	find_package(CUDA QUIET)
#	if(CUDA_FOUND)
#		message("CUDA found")
#		include_directories(${CUDA_INCLUDE_DIRS})
#	else()
#		message("CUDA not found")
#	endif()

#	include(FindCUDA/select_compute_arch)
#	cuda_detect_installed_gpus(INSTALLED_GPU_CCS_1)
#	string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
#	string(
#		REPLACE " "
#				";"
#				INSTALLED_GPU_CCS_3
#				"${INSTALLED_GPU_CCS_2}"
#	)
#	string(
#		REPLACE "."
#				""
#				CUDA_ARCH_LIST
#				"${INSTALLED_GPU_CCS_3}"
#	)
#	set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})

endif()

enable_testing()
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.17
include(CTest)

configure_file("config.hpp.in" ${CMAKE_BINARY_DIR}/config.hpp)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
	axpy.cpp
	copy.cpp
	dot.cpp
	herk.cpp
	gemv.cpp
	gemm.cpp
	numeric.cpp
	scal.cpp
	traits.cpp
	trsm.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		#set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
		target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()
	# target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)

	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

    target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})

	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
	if(NOT ENABLE_CUDA
	   AND (NOT
			DART_COMPILER_NAME
			STREQUAL
			"nvcc"
		   )
	)
		target_compile_options(
			${TEST_EXE}
			PRIVATE -Werror
					-Wall
					$<$<CXX_COMPILER_ID:GNU>:
					-Wextra
					-Wpedantic
					-Wformat-truncation
					-Wno-unknown-pragmas> # -Wconversion
					$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
					-Wextra
					-Wpedantic
					-Wmove>
					$<$<CXX_COMPILER_ID:Intel>:
					-Wextra
					-wd161
					-diag-disable=remark
					-Warray-bounds
					-Wchar-subscripts
					-Wcomment
					-Wenum-compare
					-Wformat
					-Wuninitialized
					-Wmaybe-uninitialized
					-Wmain
					-Wnarrowing
					-Wnonnull
					-Wparentheses
					-Wpointer-sign
					-Wreorder
					-Wno-return-type
					-Wsign-compare
					-Wsequence-point
					-Wtrigraphs
					-Wunused-function
					-Wunused-but-set-variable
					-Wunused-variable
					-Wwrite-strings
					-diag-error:3846>
                    $<$<OR:$<CXX_COMPILER_ID:PGI>,$<CXX_COMPILER_ID:NVHPC>>:  # EDG diagnostics list: http://www.ssl.berkeley.edu/~jimm/grizzly_docs/SSL/opt/intel/cc/9.0/lib/locale/en_US/mcpcom.msg
                    --display_error_number --diag_error=incompatible_assignment_operands,returning_ptr_to_local_variable,subscript_out_of_range,used_before_set,undefined_preproc_id,implicit_func_decl,implicit_return_from_non_void_function,missing_type_specifier
		            >
					$<$<CXX_COMPILER_ID:MSVC>:
					/W4>
		)
	endif()
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
