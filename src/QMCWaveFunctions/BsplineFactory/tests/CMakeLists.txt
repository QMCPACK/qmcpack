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

SET(SRCS test_BsplineSet.cpp test_BsplineDeviceCPU.cpp test_SplineR2RAdoptor.cpp)
#test_SplineC2CAdoptorBatched.cpp)
SET(BSPLINEFACTORY_LIBS qmc qmcwfs qmcbase qmcutil)
IF(HAVE_CUDA)
  SET(SRCS ${SRCS} test_BsplineDeviceCUDA.cpp test_BsplineCUDA.cu)
  CUDA_ADD_LIBRARY(bspline_cudatests test_BsplineCUDA.cu)
  SET(BSPLINEFACTORY_LIBS ${BSPLINEFACTORY_LIBS} bspline_cudatests)
ENDIF(HAVE_CUDA)

IF(QMC_COMPLEX)
  SET(SRCS ${SRCS} test_BsplineSetComplex.cpp)
ENDIF(QMC_COMPLEX)

MESSAGE("###Adding BsplineFactoryTests with srcs: ${SRCS}")
MESSAGE("###                           with LIBS: ${BSPLINEFACTORY_LIBS}")

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} ${BSPLINEFACTORY_LIBS} ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")
