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

set(SRC_DIR afqmc_estimators)
set(UTEST_EXE test_${SRC_DIR})

set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${UTEST_EXE} test_estimators.cpp)
target_link_libraries(${UTEST_EXE} catch_main afqmc)

list(LENGTH AFQMC_UNIT_TEST_INPUTS NUM_PAIRS)
math(EXPR ENDP "${NUM_PAIRS}-1")
foreach(I RANGE 0 ${ENDP} 2)
  list(GET AFQMC_UNIT_TEST_INPUTS ${I} HAMIL_FILE)
  math(EXPR J "${I}+1")
  list(GET AFQMC_UNIT_TEST_INPUTS ${J} WFN_FILE)
  get_filename_component(HAMIL ${HAMIL_FILE} NAME_WE)
  get_filename_component(WFN ${WFN_FILE} NAME_WE)
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_${HAMIL}_${WFN})
  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}> "--hamil ${HAMIL_FILE}" "--wfn ${WFN_FILE}")
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  set_property(
    TEST ${UTEST_NAME}
    APPEND
    PROPERTY LABELS "afqmc")
endforeach()
