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

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR afqmc_hamiltonian_operations)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/afqmc/HamiltonianOperations)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

ADD_EXECUTABLE(${UTEST_EXE} test_hamiltonian_operations.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main afqmc)
IF (QMC_COMPLEX)
  SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_uc_wfn_rhf)
  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}"
                "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/ham_chol_uc.h5"
                "--wfn ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/wfn_rhf.h5")
#ELSE()
#  SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_wfn_rhf)
#  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}"
#      "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/ham_chol.h5"
#      "--wfn ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/wfn_rhf.h5")
ENDIF()
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
SET_PROPERTY(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")
