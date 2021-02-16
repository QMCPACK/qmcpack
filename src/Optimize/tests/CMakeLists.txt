#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR optimize)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

ADD_EXECUTABLE(${UTEST_EXE} test_variable_set.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcutil)

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
