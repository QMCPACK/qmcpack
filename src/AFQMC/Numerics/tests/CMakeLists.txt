#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

MESSAGE("Adding this unit test")


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR afqmc_numerics)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})


ADD_EXECUTABLE(${UTEST_EXE} test_sparse_ops.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} ${QMC_UTIL_LIBS} afqmc)

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit;afqmc")

