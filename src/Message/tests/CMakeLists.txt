#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@ganl.gov Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR} )

SET(SRC_DIR message)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

ADD_EXECUTABLE(${UTEST_EXE} test_communciate.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} PUBLIC message catch_main)

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")

IF(HAVE_MPI)
  SET(UTEST_EXE test_${SRC_DIR}_mpi)
  SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)
  #this is dependent on the directory creation and sym linking of earlier driver tests
  SET(MPI_UTILITY_TEST_SRC test_mpi_exception_wrapper.cpp)
  ADD_EXECUTABLE(${UTEST_EXE} ${MPI_UTILITY_TEST_SRC})
  #Way too many depenedencies make for very slow test linking
  TARGET_LINK_LIBRARIES(${UTEST_EXE} PUBLIC message catch_main)
  # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
  ADD_UNIT_TEST(${UTEST_NAME} 3 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ENDIF(HAVE_MPI)
