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

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR} )

SET(SRC_DIR io_hdf5)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# Directory where input file is copied to, and working directory for unit test
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/io/hdf)

SET(UTEST_SRC test_hdf_archive.cpp  test_hdf_parallel.cpp test_hdf_reshape.cpp test_hdf_hyperslab.cpp)
IF(BUILD_AFQMC)
  SET(UTEST_SRC ${UTEST_SRC} test_hdf_multiarray.cpp)
ENDIF(BUILD_AFQMC)

ADD_EXECUTABLE(${UTEST_EXE} ${UTEST_SRC})
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcio)

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
