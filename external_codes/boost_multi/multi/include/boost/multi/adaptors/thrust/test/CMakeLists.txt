cmake_minimum_required(VERSION 3.18)
project(
	boost-multi-adaptor-cuda-thrust-test
	VERSION 0.1
	LANGUAGES CXX CUDA
)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
	set(CMAKE_CUDA_ARCHITECTURES native)
endif()

# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--display_error_number\"")

find_package(Boost)  #  COMPONENTS boost unit_test_framework)  # tests require Boost.Core LightweightTest

include(CTest)
enable_testing()

# find_program(MEMORYCHECK_COMMAND valgrind)
# set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")

find_package(Boost COMPONENTS timer)
if(NOT Boost_FOUND)
  message(ERROR "  Cannot find Boost.Timer. Try\n    sudo apt install libboost-timer-dev  # in Debian/Ubuntu")
endif()

set(TEST_SRCS
    array.cu
    memory_resource.cu
	run_reduce.cu
	set_identity_kernel.cu
    speed.cu
    speed_algo.cu
    universal.cu
	vector.cu
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		set_property(TARGET ${TEST_EXE} PROPERTY CUDA_STANDARD 17)
	endif()

	target_include_directories(${TEST_EXE}        PUBLIC ../../../../../include)

	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
    if(ENABLE_CUDA)
	    add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
    endif()
endforeach()

set_tests_properties(speed.cu.x PROPERTIES RUN_SERIAL TRUE)
