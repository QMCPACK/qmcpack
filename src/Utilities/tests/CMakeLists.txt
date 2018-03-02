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

SET(SRC_DIR utilities)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})

ADD_EXECUTABLE(${UTEST_EXE} test_rng.cpp test_parser.cpp test_timer.cpp
                            test_prime_set.cpp test_partition.cpp test_pooled_memory.cpp
                            test_infostream.cpp test_output_manager.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")
