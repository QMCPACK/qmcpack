#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

set(SRC_DIR spline2)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})
#SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/hamiltonians)

set(SRCS test_multi_spline.cpp)

#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple.txt" ${UTEST_DIR})
#EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/C.BFD.xml" ${UTEST_DIR})

add_executable(${UTEST_EXE} ${SRCS})
target_link_libraries(${UTEST_EXE} catch_main einspline qmcutil)

add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
