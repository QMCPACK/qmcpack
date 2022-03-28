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

message(VERBOSE "Adding AFQMC Slater Determinant Operations unit tests")

set(SRC_DIR afqmc_slaterdeterminantoperations)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(${UTEST_EXE} test_sdet_ops.cpp)
target_link_libraries(${UTEST_EXE} catch_main afqmc)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(
  TEST ${UTEST_NAME}
  APPEND
  PROPERTY LABELS "afqmc")
