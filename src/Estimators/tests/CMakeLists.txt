#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////
    
    

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR estimators)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})
#SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/hamiltonians)

SET(SRCS test_accumulator.cpp test_local_energy_est.cpp test_manager.cpp)

#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple.txt" ${UTEST_DIR})
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${qmcpack_SOURCE_DIR}/pseudopotentials/BFD/C.BFD.xml" ${UTEST_DIR})


ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcwfs qmcham qmcbase qmcutil qmcdriver ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
#SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit" WORKING_DIRECTORY ${UTEST_DIR})
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")

