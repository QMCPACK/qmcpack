#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(container_testing_src makeRngSpdMatrix.cpp)
add_library(container_testing ${container_testing_src})
target_include_directories(container_testing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(container_testing utilities_for_test containers)

# This test is always built and works due to openmp pragmas becoming no ops
# even in versions of the code that do not use OMPAllocator.

set(SRC_DIR dualallocators)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-integration_${SRC_DIR})

set(INTEGRATION_TEST_SRC test_dual_allocators_ohmms_containers.cpp)

add_executable(${UTEST_EXE} ${INTEGRATION_TEST_SRC})
target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(${UTEST_EXE} catch_main containers platform_LA container_testing utilities_for_test)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
