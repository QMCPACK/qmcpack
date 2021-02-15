#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})
SET(SRC_DIR MinimalContainers)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

SET(SRCS test_ConstantSizeMatrix.cpp test_RecordArray.cpp)

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})

TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main containers)

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
