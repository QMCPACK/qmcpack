#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR} )

SET(TEST_NAME utilities_for_testing)
SET(UTEST_EXE test_${TEST_NAME})
SET(UTEST_NAME deterministic-unit_test_${TEST_NAME})

ADD_EXECUTABLE(${UTEST_EXE} test_checkMatrix.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main utilities_for_test containers)
target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
