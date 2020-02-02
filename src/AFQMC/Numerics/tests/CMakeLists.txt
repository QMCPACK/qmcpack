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

MESSAGE_VERBOSE("Adding AFQMC numerics unit tests")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR afqmc_numerics)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})


#IF(ENABLE_CUDA)
#  CUDA_ADD_EXECUTABLE(${UTEST_EXE} test_dense_numerics.cpp test_ma_blas.cpp
#test_sparse_numerics.cpp test_sparse_numerics_native.cpp)
#ELSE(ENABLE_CUDA)
  ADD_EXECUTABLE(${UTEST_EXE} test_dense_numerics.cpp test_ma_blas.cpp
test_sparse_numerics.cpp test_sparse_numerics_native.cpp)
#ENDIF(ENABLE_CUDA)
TARGET_LINK_LIBRARIES(${UTEST_EXE} afqmc qmcutil ${QMC_UTIL_LIBS})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})
SET_PROPERTY(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")


