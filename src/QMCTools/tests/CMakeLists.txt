#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers
#//
#// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
#//
#// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR tools)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/${SRC_DIR})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

SET(SRCS test_qmcfstool.cpp)
         
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple_Sk.dat" ${UTEST_DIR})
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple_input.xml" ${UTEST_DIR})

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main fstool)
IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcparticle qmcutil containers platform_omptarget)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

