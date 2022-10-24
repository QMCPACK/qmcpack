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

# TODO(correaa) this will be simpler in Cmake 3.24 https://stackoverflow.com/a/72332793/225186
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
endif()

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)

set(TEST_SRCS
    array.cu
    speed.cu
    speed_algo.cu
    universal.cu
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
	endif()

	target_include_directories(${TEST_EXE}        PUBLIC ../../../../../include)

	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
    if("30;35;50;52;60;61;70;70+PTX;75;75+PTX" STREQUAL "${CUDA_ARCH_LIST}" OR "" STREQUAL "${CUDA_ARCH_LIST}" OR NOT ENABLE_GPU)
	    message(STATUS "gpu not detected ${CUDA_ARCH_LIST}")
    else()
	    message(STATUS "gpu detected: ${CUDA_ARCH_LIST}")
	    add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
    endif()
    if(ENABLE_GPU)
	    message(STATUS "gpu run forced")
	    add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
    endif()
endforeach()
