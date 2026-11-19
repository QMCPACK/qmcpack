cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

#project(
#   boost-multi-adaptors-fftw-test
#   VERSION 0.1
#   LANGUAGES CXX
#)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(CMAKE_COMPILE_WARNING_AS_ERROR "Compile warnings as errors" ON)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
endif()

if(ENABLE_CIRCLE)
find_package(BLAS REQUIRED)
find_path(
   BLAS_INCLUDE_DIRS
   cblas.h
   /usr/include
   /usr/local/include
   $ENV{BLAS_HOME}/include
)
include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})
link_libraries(${BLAS_LIBRARIES})
endif()

enable_testing()

include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
  combinations.cpp
# copy.cpp
  core.cpp
  direct_vs_transpose.cpp
  fft.cpp
  moved.cpp
  shift.cpp
  so_shift.cpp
# thrust.cpp
  transpose.cpp
  transpose_square.cpp
)

foreach(TEST_FILE ${TEST_SRCS})

	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		# target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()

	target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})

	target_link_libraries(${TEST_EXE} PUBLIC multi)
	target_link_libraries(${TEST_EXE} PRIVATE multi_fftw3)  # PkgConfig::FFTW)

	if(NOT DEFINED ENABLE_CIRCLE)
		target_compile_options(
			${TEST_EXE}
			PRIVATE
				$<$<AND:$<CXX_COMPILER_ID:GNU>,$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>,$<NOT:$<CUDA_COMPILER_ID:Clang>>>:
					-Werror
					-Wall
					-Wextra
					-fno-common
					-Wpedantic
					-Wformat-truncation
					-Wparentheses
					-fstack-usage
				>
				$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
					-Wall
					-Wextra
					-Wpedantic
					# -Wno-exit-time-destructors  # bug in cmake 4, boost,
					# -Wno-poison-system-directories  # bug in cmake 4, boost,
					-Wno-documentation-unknown-command  # array_ref.hpp uses \mathrm{} etc. in doc comments, not a recognized Doxygen/clang command
					#-fno-common
					#-Wmove
					#-Wno-error=\#warnings
				>
				$<$<CXX_COMPILER_ID:Intel>:
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
					-Werror
					-diag-error:3846
				>
				$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
					#-Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wdouble-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
					#--Werror=cross-execution-space-call
					#--Werror ext-lambda-captures-this
					#--use_fast_math
					#--extra-device-vectorization
					#--restrict
					-Xcudafe=--display_error_number -Xcudafe=--diag_suppress=20011 -Xcudafe=--diag_suppress=20014 -Xcudafe=--diag_suppress=20011 -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier -Xcudafe=--diag_error=invalid_error_tag
				>
				$<$<CXX_COMPILER_ID:MSVC>:
					$<$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:
						/GS      # Enables buffer security check (can help detect stack-based buffer overflows)
						/Wall    # /volatile:iso /WX /Zc:wchar_t /Zc:forScope /Zc:inline /W4
						/wd4371  # layout of class may have changed from a previous version of the compiler due to better packing of member
						/wd4514  # unreferenced inline function has been removed
						/wd4623  # default constructor was implicitly defined as deleted #TODO(correaa) investigate
						/wd4625  # copy constructor was implicitly defined as deleted
						/wd4626  # assignment operator was implicitly defined as deleted #TODO(correaa) investigate
						/wd4710  # function not inlined
						/wd4711  # function selected for automatic inline expansion
						/wd4820  # bytes padding added after base class
						/wd4866  # compiler may not enforce left-to-right evaluation order for call
						/wd4848  # support for attribute [[msvc::no_unique_address]] in C++17 and earlier is a vendor extension
						/wd4868  # compiler may not enforce left-to-right evaluation order in braced initializer list
						/wd5026  # move constructor was implicitly defined as deleted
						/wd5027  # move assignment operator was implicitly defined as deleted
						/wd5045  # Compiler will insert Spectre mitigation for memory load if /Qspectre switch specified
						# /wd5246  # '_Elems': the initialization of a subobject should be wrapped in braces
						/wd4464  # relative include path contains '..' (triggered by CUDA Toolkit's own cooperative_groups headers)
						/wd5039  # pointer or reference to potentially throwing function passed to 'extern "C"' function under -EHc (from Windows SDK's winbase.h)
					>
					$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
						-Xcompiler=/GS
						-Xcompiler=/Wall
						-Xcompiler=/wd4371
						-Xcompiler=/wd4514
						-Xcompiler=/wd4623
						-Xcompiler=/wd4625
						-Xcompiler=/wd4626
						-Xcompiler=/wd4710
						-Xcompiler=/wd4711
						-Xcompiler=/wd4820
						-Xcompiler=/wd4866
						-Xcompiler=/wd4848
						-Xcompiler=/wd4868
						-Xcompiler=/wd5026
						-Xcompiler=/wd5027
						-Xcompiler=/wd5045
						-Xcompiler=/wd4464
						-Xcompiler=/wd4619  # #pragma warning: there is no warning number 'N' (from CUDA Toolkit's own cuda/std/cstddef)
						-Xcompiler=/wd5031  # #pragma warning(pop): likely mismatch, popping warning state pushed in different file (from CUDA Toolkit's cuda_fp16.hpp/cuda_bf16.hpp)
						-Xcompiler=/wd5032  # #pragma warning(push) with no corresponding #pragma warning(pop) (nvcc host-code splitting desyncs cl.exe's pragma stack across boost/multi/array_ref.hpp)
						-Xcompiler=/wd4365  # signed/unsigned mismatch (from CUDA Toolkit's own cuda_fp8.hpp)
						-Xcompiler=/wd4061  # enumerator in switch of enum is not explicitly handled by a case label (from CUDA Toolkit's own cuda/std/__cuda/ensure_current_device.h switching over cudaError)
						-Xcompiler=/wd5220  # non-static data member with a volatile qualified type no longer implies (from CUDA Toolkit's own nvtx3/nvtxDetail/nvtxImpl.h)
						-Xcompiler=/wd4296  # '<': expression is always false (from CUDA Toolkit's own cuda/std/__type_traits/is_signed.h)
						-Xcompiler=/wd5246  # the initialization of a subobject should be wrapped in braces (from CUDA Toolkit's own cuda/std/__exception/cuda_error.h)
						-Xcompiler=/wd4555  # result of expression not used (from Windows SDK's winnt.h)
						-Xcompiler=/wd4191  # unsafe conversion from 'FARPROC' to function pointer type (from CUDA Toolkit's own nvtx3/nvtxDetail/nvtxInit.h)
						-Xcompiler=/wd5039  # pointer or reference to potentially throwing function passed to 'extern "C"' function under -EHc (from Windows SDK's winbase.h)
					>
				>
		)
	endif()

	add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
endforeach()
