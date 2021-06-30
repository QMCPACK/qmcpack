#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

set(SRC_DIR optimize)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(${UTEST_EXE} test_variable_set.cpp)
target_link_libraries(${UTEST_EXE} catch_main qmcutil)

add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
