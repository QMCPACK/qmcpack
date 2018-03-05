#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR drivers)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})


SET(DRIVER_TEST_SRC test_vmc.cpp test_dmc.cpp test_drift.cpp)
IF(HAVE_MPI)
  SET(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_walker_control.cpp)
ENDIF()
ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})

USE_FAKE_RNG(${UTEST_EXE})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcdriver_unit qmcham qmcwfs qmcbase qmcutil qmcfakerng ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")
