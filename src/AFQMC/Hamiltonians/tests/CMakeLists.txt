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

set(SRC_DIR afqmc_hamiltonians)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${UTEST_EXE} test_hamiltonian_factory.cpp)
target_link_libraries(${UTEST_EXE} catch_main afqmc)

set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_uc)

if(QMC_COMPLEX)
  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>
                "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/ham_chol_uc.h5")
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  set_property(
    TEST ${UTEST_NAME}
    APPEND
    PROPERTY LABELS "afqmc")
  if(NOT ENABLE_CUDA)
    set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol_sc)
    add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>
                  "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5")
    set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
    set_property(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")

    set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_thc_sc)
    add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>
                  "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_szv/ham_thc_sc.h5")
    set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
    set_property(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")
  endif()
else()
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_ham_chol)
  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>
                "--hamil ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/ham_chol.h5")
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  set_property(TEST ${UTEST_NAME} APPEND PROPERTY LABELS "afqmc")
endif()
