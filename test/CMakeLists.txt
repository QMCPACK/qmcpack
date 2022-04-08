#[=[Multi Test suite can be run like this: -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
  mkdir -p build && cd build
  cmake .. [-DENABLE_CUDA=1]
  make -j
  ctest -j --output-on-error [-T memcheck]
  exit
#]=]
cmake_minimum_required(VERSION 3.11)

if(EXISTS "CMakeLists.txt")
	message(FATAL_ERROR "You probably don't want to run CMake in a directory with CMakeLists.txt")
endif()

project(
	boost-multi-test
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost REQUIRED COMPONENTS unit_test_framework serialization)

if((ENABLE_CUDA OR DEFINED CXXCUDA) AND NOT DEFINED FORCED_ARCH)
	enable_language(CUDA)

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
include(CTest)

file(
	GLOB TEST_SRCS
	RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
	*.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})

	if((ENABLE_CUDA OR DEFINED CXXCUDA) AND NOT DEFINED FORCED_ARCH)
#		set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
	endif()

	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
	endif()

	target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_libraries(${TEST_EXE} PRIVATE Boost::unit_test_framework Boost::serialization)

	target_compile_definitions(${TEST_EXE} PRIVATE "$<$<CUDA_COMPILER_ID:NVIDIA>:BOOST_PP_VARIADICS>")

	target_compile_options(
		${TEST_EXE}
		PRIVATE
			$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
				-Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wdouble-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
				--expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe=--display_error_number -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
			>
			$<$<AND:$<CXX_COMPILER_ID:GNU>,$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>,$<NOT:$<CUDA_COMPILER_ID:Clang>>>:
				-Werror
				-Wall
				-Wextra
				-Wcast-align
				-Wcast-qual
				-Wconversion
				-Wdouble-promotion
				-Wduplicated-branches
				-Wduplicated-cond
				-Wformat-truncation
				-Wformat=2
				-Wformat-security
				-Wlogical-op
				-Wmisleading-indentation
				-Wmissing-include-dirs
				-Wnon-virtual-dtor
				-Wno-missing-declarations
				-Wno-unknown-pragmas
				-Wnon-virtual-dtor
				-Wnull-dereference
				-Wold-style-cast
				-Woverloaded-virtual
				-Wpointer-arith
				-Wpedantic
				-Wredundant-decls
				-Wshadow
				-Wold-style-cast
				-Woverloaded-virtual
				-Wpedantic
				-Wsign-conversion
				-Wstrict-aliasing
				-Wswitch-enum
				-Wtrampolines
				-Wuninitialized
				-Wunused
				-Wunused-but-set-variable
				-Wunused-result
				-Wuseless-cast
				-Wzero-as-null-pointer-constant
			>
			$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CUDA_COMPILER_ID:Clang>>:
				-Werror
				-Wall
				-Wextra
				-Wcast-align
				-Wcast-qual
				-Wconversion
				-Wdouble-promotion
				-Wformat=2
				-Wmissing-include-dirs
				-Wno-unknown-pragmas
				-Wnon-virtual-dtor
				-Wnull-dereference
				-Wmove
				-Wold-style-cast
				-Woverloaded-virtual
				-Wpointer-arith
				-Wpedantic
				-Wredundant-decls
				-Wshadow
				-Wold-style-cast
				-Woverloaded-virtual
				-Wpedantic
				-Wsign-conversion
				-Wstrict-aliasing
				-Wswitch-enum
				-Wuninitialized
				-Wunused
				-Wunused-result
				-Wzero-as-null-pointer-constant
            >
			$<$<CXX_COMPILER_ID:Intel>:  # also IntelLLVM, XL (ibm), XLClang (ibm)
				-Werror
				-Wall
				-Wextra
				-diag-disable=remark
				-diag-error:3846
				-diag-disable=1011  # disables warning missing return at the end of non-void function
				-wd161
				-Wabi
				-Warray-bounds
				-Wcast-qual
				-Wchar-subscripts
				-Wcomment
				-Wdeprecated
				-Wenum-compare
				-Wextra-tokens
				-Wformat
				-Wformat=2
				-Wformat-security
				-Wic-pointer
				-Wnarrowing
				-Wno-return-type
				-Wnon-virtual-dtor
				-Wnonnull
				-Wmaybe-uninitialized
				-Wmain
				-Wmissing-declarations
				-Wmissing-prototypes
				-Wmultichar
				-Woverloaded-virtual
				-Woverflow
				-Wparentheses
				-Wpointer-arith
				-Wpointer-sign
				-Wreorder
				-Wreturn-type
				-Wsign-compare
				-Wsequence-point
				-Wshadow
				-Wshorten-64-to-32
				-Wno-strict-aliasing
				-Wstrict-prototypes
				-Wtrigraphs
				-Wuninitialized
				-Wunused
				-Wunused-but-set-variable
				-Wunused-function
				-Wunused-parameter
				-Wunused-variable
				-Wwrite-strings
			>
			$<$<OR:$<CXX_COMPILER_ID:PGI>,$<CXX_COMPILER_ID:NVHPC>>:
				-Werror
				-Wall
				-Wcast-qual
				-Wformat=2
				-Wshadow
				-Wunused
			>
			$<$<CXX_COMPILER_ID:MSVC>:
				/WX
				/W4
			>
	)
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
