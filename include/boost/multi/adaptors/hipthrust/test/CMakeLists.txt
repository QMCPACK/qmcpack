cmake_minimum_required(VERSION 3.24)

project(boost-multi-adaptor-thrust-hip LANGUAGES CXX HIP)

find_package(Boost     REQUIRED COMPONENTS timer unit_test_framework)
find_package(rocthrust REQUIRED)

include(CMakePrintHelpers)

cmake_print_properties(TARGETS roc::rocthrust PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

enable_testing()

include(CTest)

# set(TEST_SRCS
# #    array.cu
# #   memory_resource.cu
# #   set_identity_kernel.cu
# #    speed.cu
# #    speed_algo.cu
# #    universal.cu
#   vector.cpp
# )

# add_executable(vector.cu.nvidia.x vector.cu)
# set_source_files_properties(vector.cu PROPERTIES LANGUAGE CUDA)

# target_link_libraries(vector.cu.nvidia.x PRIVATE Boost::unit_test_framework)
# add_test(NAME vector.cu.nvidia.x COMMAND vector.cu.nvidia.x)

add_executable(vector.hip.amd.x vector.hip)
set_source_files_properties(vector.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(vector.hip.amd.x PRIVATE Boost::unit_test_framework roc::rocthrust)
target_compile_definitions(vector.hip.amd.x PRIVATE BOOST_TEST_DYN_LINK=1)

add_test(NAME vector.hip.amd.x COMMAND $<TARGET_FILE:vector.hip.amd.x>)

add_executable(array.hip.amd.x array.hip)
set_source_files_properties(array.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(array.hip.amd.x PRIVATE Boost::unit_test_framework Boost::timer multi roc::rocthrust)
target_compile_definitions(array.hip.amd.x PRIVATE BOOST_TEST_DYN_LINK=1)

add_test(NAME array.hip.amd.x COMMAND $<TARGET_FILE:array.hip.amd.x>)

add_executable(speed.hip.amd.x speed.hip)
set_source_files_properties(speed.hip PROPERTIES LANGUAGE HIP)
target_link_libraries(speed.hip.amd.x PRIVATE Boost::unit_test_framework Boost::timer multi roc::rocthrust)
target_compile_definitions(speed.hip.amd.x PRIVATE BOOST_TEST_DYN_LINK=1)
target_compile_definitions(speed.hip.amd.x PRIVATE BOOST_TEST_MODULE="C++ Unit Tests for Multi CUDA thrust universal copy and assignment")

add_test(NAME speed.hip.amd.x COMMAND $<TARGET_FILE:speed.hip.amd.x>)

# foreach(TEST_FILE ${TEST_SRCS})
# # set(TEST_EXE "${TEST_FILE}")
#   add_executable(vector ${TEST_FILE})
# # if(ENABLE_CUDA OR DEFINED CXXCUDA)
# #   set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE HIP)
# # endif()

# # target_include_directories(${TEST_EXE}        PRIVATE ../../../../../include)

# #    if(ENABLE_CUDA)
# #     add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
# #    endif()
# endforeach()

# set_tests_properties(speed.cu.x PROPERTIES RUN_SERIAL TRUE)
