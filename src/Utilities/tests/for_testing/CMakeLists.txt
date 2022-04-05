#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(TEST_NAME utilities_for_testing)
set(UTEST_EXE test_${TEST_NAME})
set(UTEST_NAME deterministic-unit_test_${TEST_NAME})

add_executable(${UTEST_EXE} test_checkMatrix.cpp test_RandomForTest.cpp)
target_link_libraries(${UTEST_EXE} catch_main_no_mpi utilities_for_test containers)
target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
