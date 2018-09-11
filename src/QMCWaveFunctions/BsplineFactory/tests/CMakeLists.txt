#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})
SET(SRC_DIR wavefunction_bsplinefactory)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})

SET(SRCS test_BsplineSet.cpp test_BsplineDeviceCUDA.cpp test_BsplineDeviceCPU.cpp test_SplineR2RAdoptorBatched.cpp)
IF(HAVE_CUDA)
  SET(SRCS ${SRCS}
      test_BsplineCUDA.cu)
  CUDA_ADD_LIBRARY(bspline_cudatests test_BsplineCUDA.cu)
ENDIF()

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} bspline_cudatests qmc qmcwfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")
