cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
	add_definitions(-DRETURN_BY_STACK)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
endif()

enable_testing()
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
	complex.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		#set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
		target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)

	target_link_libraries(${TEST_EXE} PRIVATE multi)

	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

	# target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE BOOST_TEST_DYN_LINK=1)
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_compile_definitions(${TEST_EXE} PRIVATE BOOST_TEST_MODULE="C++ Unit Tests for Multi complex")


	if(NOT ENABLE_CUDA
	   AND (NOT
			DART_COMPILER_NAME
			STREQUAL
			"nvcc"
		   )
		AND (NOT DEFINED ENABLE_CIRCLE)
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
					-Wmove
					-Wno-error=\#warnings
					>
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
	else()
		target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_TEST_DYN_LINK=1")
		# target_compile_options    (${TEST_EXE} PRIVATE -Werror -Wall)
	endif()
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
