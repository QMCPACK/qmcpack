#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Miguel A. Morales, LLNL
#//
#// File created by: Fionn Malone, malone14@llnl.gov, LLNL
#//////////////////////////////////////////////////////////////////////////////////////

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR afqmc_estimators)
SET(UTEST_EXE test_${SRC_DIR})

SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/afqmc/Estimators)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

ADD_EXECUTABLE(${UTEST_EXE} test_estimators.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main afqmc)

LIST(LENGTH AFQMC_UNIT_TEST_INPUTS NUM_PAIRS)
MATH(EXPR ENDP "${NUM_PAIRS}-1")
FOREACH(I RANGE 0 ${ENDP} 2)
    LIST(GET AFQMC_UNIT_TEST_INPUTS ${I} HAMIL_FILE)
    MATH(EXPR J "${I}+1")
    LIST(GET AFQMC_UNIT_TEST_INPUTS ${J} WFN_FILE)
    GET_FILENAME_COMPONENT(HAMIL ${HAMIL_FILE} NAME_WE)
    GET_FILENAME_COMPONENT(WFN ${WFN_FILE} NAME_WE)
    SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_${HAMIL}_${WFN})
    ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}"
                  "--hamil ${HAMIL_FILE}" "--wfn ${WFN_FILE}")
    SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
    SET_PROPERTY(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")
ENDFOREACH()
