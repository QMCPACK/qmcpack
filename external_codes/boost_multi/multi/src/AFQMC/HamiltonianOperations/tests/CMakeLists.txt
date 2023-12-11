#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Miguel A. Morales, LLNL
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR afqmc_hamiltonian_operations)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${UTEST_EXE} test_hamiltonian_operations.cpp)
target_link_libraries(${UTEST_EXE} catch_main afqmc)
if(QMC_COMPLEX)
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_uc_wfn_rhf)
  add_unit_test(
    ${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>
    "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/ham_chol_uc.h5"
    "--wfn ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/wfn_rhf.h5")
  #ELSE()
  #  SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_wfn_rhf)
  #  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}"
  #      "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/ham_chol.h5"
  #      "--wfn ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/wfn_rhf.h5")
endif()
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
set_property(
  TEST ${UTEST_NAME}
  APPEND
  PROPERTY LABELS "afqmc")
