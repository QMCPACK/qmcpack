# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.11)
project(
	boost-multi-adaptor-cuda-thrust-test
	VERSION 0.1
	LANGUAGES CXX CUDA
)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--display_error_number\"")

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

enable_testing()

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")
include(CTest)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
	# temporary workaround to detect CUDA arch, for between the deprecation of FindCUDA in CMake 3.10 and the availability of FindCudaToolkit in CMake 3.17 # https://stackoverflow.com/a/68223399/225186
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
	set_property(GLOBAL PROPERTY CUDA_ARCHITECTURES "${CUDA_ARCH_LIST}")
endif()

if("30;35;50;52;60;61;70;70+PTX;75;75+PTX" STREQUAL "${CUDA_ARCH_LIST}" OR "" STREQUAL "${CUDA_ARCH_LIST}")
	message(STATUS "gpu not detected ${CUDA_ARCH_LIST}")
	set(TEST_SRCS)
else()
	message(STATUS "gpu detected: ${CUDA_ARCH_LIST}")
	set(TEST_SRCS array.cu vector.cu)
endif()

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
	endif()
	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})
	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
	# if(NOT ENABLE_CUDA AND (NOT DART_COMPILER_NAME STREQUAL "nvcc")) target_compile_options  (${TEST_EXE} PRIVATE -Werror $<$<CXX_COMPILER_ID:GNU>: -Wall -Wshadow -Wsign-conversion -Wold-style-cast -Wcast-align -Wunused -Wnon-virtual-dtor
	# -Woverloaded-virtual -Wcast-qual -Wconversion -Wnull-dereference -Wformat=2 -Wextra -Wpedantic -Wuninitialized -Wformat-truncation -Wdouble-promotion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op>
	# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
	# -Wextra -Wpedantic -Wmove> $<$<CXX_COMPILER_ID:Intel>: -Wextra                            -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain
	# -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -diag-error:3846>
	# $<$<CXX_COMPILER_ID:PGI>: #PGI 20  doesn't recognize the -Wextra option --diag_suppress=implicit_return_from_non_void_function> $<$<CXX_COMPILER_ID:MSVC>: /WX /W4>) endif() if(NOT ENABLE_CUDA) target_compile_options  (${TEST_EXE}
	# PRIVATE -Werror -Wall -Wextra -fno-common $<$<CXX_COMPILER_ID:GNU>: -Wpedantic -Wformat-truncation -fstack-usage>#-Wconversion
	# $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
	# -Wpedantic -Wmove> $<$<CXX_COMPILER_ID:Intel>: -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses
	# -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Werror -diag-error:3846 > $<$<CXX_COMPILER_ID:MSVC>: /W4>) endif()
	add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()
