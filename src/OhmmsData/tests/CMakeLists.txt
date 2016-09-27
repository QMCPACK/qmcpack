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

SET(SRC_DIR ohmmsdata)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})

# Directory where input file is copied to, and working directory for unit test
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/xml_test)

ADD_EXECUTABLE(${UTEST_EXE} test_xml.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcutil ${QMC_UTIL_LIBS} )


EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple.xml" ${UTEST_DIR})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}" WORKING_DIRECTORY  ${UTEST_DIR})
#SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit" WORKING_DIRECTORY  ${UTEST_DIR})
