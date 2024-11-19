cmake_minimum_required(VERSION 3.18)  # required by BLAS::BLAS

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
	find_package(Boost REQUIRED)  # COMPONENTS boost)  # headers unit_test_framework)
#   link_libraries("-lboost_unit_test_framework")

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
	# link_libraries("-lboost_unit_test_framework")

	link_libraries("-lblas") # cmake cannot detect BLAS with pgi/nvc++ but it ships with its own version
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
#   add_definitions(-DRETURN_BY_STACK)
	add_definitions(-DBLAS_DOT_RETURNS_VOID)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")  # nvhpc will find its own blas version
#   add_definitions(-DRETURN_BY_STACK)
	add_definitions(-DBLAS_DOT_RETURNS_VOID)
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
endif()

enable_testing()
include(CTest)

# configure_file("config.hpp.in" ${CMAKE_BINARY_DIR}/config.hpp)

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
	nrm2.cpp
	scal.cpp
	syrk.cpp
	traits.cpp
	trsm.cpp
)

find_package(TBB)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		#set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
		set_property(TARGET ${TEST_EXE} PROPERTY CUDA_STANDARD 17)
	endif()

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)

	# target_link_libraries     (${TEST_EXE} PRIVATE Boost::boost)  # headers)   # unit_test_framework   )
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS} )

	if(TBB_FOUND)
		target_link_libraries(${TEST_EXE} PRIVATE TBB::tbb)
	endif()

	# target_compile_options(
	#   ${TEST_EXE} PRIVATE
	#   $<$<CXX_COMPILER_ID:IntelLLVM>:
	#   -Wno-unknown-warning-option  # for -Wno-enum-constexpr-conversion in older versions
	#   -Wno-enum-constexpr-conversion  # for Boost.Test
	#   >
	# )

	target_link_libraries(${TEST_EXE} PRIVATE ${BLAS_LIBRARIES})  # there is no BLAS_INCLUDE_DIRS

    target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_libraries(${TEST_EXE} PRIVATE BLAS::BLAS)

	if(NOT ENABLE_CIRCLE)
		if(ENABLE_CUDA)
			target_compile_options(
				${TEST_EXE}
				PRIVATE
				$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
				-Wno-unknown-warning-option
				#-Wextra
				#-Wpedantic
				#-Wmove
				-Wno-enum-constexpr-conversion  # for Boost.Test
				#-Wno-error=\#warnings
				>
			)
		endif()

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
					-Wno-unknown-pragmas
					-Wparentheses
					> # -Wconversion
					$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
					# -Weverything
					-Wextra
					-Wpedantic
					-Wmove
					-Wno-error=\#warnings
					> # -Wno-c++98-compat  # use of constexpr
					$<$<CXX_COMPILER_ID:IntelLLVM>:
					-Wno-unknown-warning-option  # for -Wno-enum-constexpr-conversion in older versions
					-Wno-enum-constexpr-conversion  # for Boost.Test
					>
					$<$<CXX_COMPILER_ID:Intel>:
					-Wextra
					-wd161
					-diag-disable=remark
					-Warray-bounds
					-Wchar-subscripts
					-Wcomment
					-Wenum-compare
					-Wno-enum-constexpr-conversion  # for Boost.Test
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
	endif()
	add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
endforeach()
