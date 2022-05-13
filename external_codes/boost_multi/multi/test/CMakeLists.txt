# -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#[=[Multi Test suite can be run like this:
  mkdir -p build
  cd build
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

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")

	include(FindCUDA/select_compute_arch)
	cuda_detect_installed_gpus(INSTALLED_GPU_CCS_1)
	string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
	string(
		REPLACE " "
				";"
				INSTALLED_GPU_CCS_3
				"${INSTALLED_GPU_CCS_2}"
	)
	string(
		REPLACE "."
				""
				CUDA_ARCH_LIST
				"${INSTALLED_GPU_CCS_3}"
	)
	set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})

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
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
	endif()
	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})
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
			PRIVATE # -fstrict-aliasing -Wstrict-aliasing
					-Werror
					-Wall
					-Wshadow
					-Wsign-conversion
					-Wunused
					-Wnon-virtual-dtor
					-Woverloaded-virtual
					-Wcast-qual
					-Wconversion
					-Wformat=2
					$<$<CXX_COMPILER_ID:GNU>:
					-Werror
					-Wall
					-Wextra
					-Wpedantic
					-Wdouble-promotion
					-Wduplicated-branches
					-Wduplicated-cond
					-Wcast-align
					-Wcast-qual
					-Wconversion
					-Wformat-truncation
					-Wformat=2
					-Wlogical-op
					-Wmisleading-indentation
					-Wno-missing-declarations
					-Wmissing-include-dirs
					-Wnon-virtual-dtor
					-Wnull-dereference
					-Wold-style-cast
					-Woverloaded-virtual
					-Wpointer-arith
					-Wshadow
					-Wredundant-decls
					-Wsign-conversion
					-Wswitch-enum
					-Wtrampolines
					-Wuninitialized
					-Wunused
					-Wunused-but-set-variable
					-Wunused-result
					-Wuseless-cast
					-Wzero-as-null-pointer-constant
					>
					$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
					-Wextra
					-Wnon-virtual-dtor
					-Wpedantic
					-Wmove
					-Wdouble-promotion
					-Wcast-align
					-Wold-style-cast
					-Wnull-dereference
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
					$<$<CXX_COMPILER_ID:PGI>: # PGI 20  doesn't recognize the -Wextra option
					--diag_suppress=implicit_return_from_non_void_function>
					$<$<CXX_COMPILER_ID:MSVC>:
					/WX
					/W4>
		)
	endif()
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
