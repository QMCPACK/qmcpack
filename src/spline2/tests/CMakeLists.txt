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

set(SRC_DIR spline2)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

set(SRCS test_multi_spline.cpp)

add_executable(${UTEST_EXE} ${SRCS})
target_link_libraries(${UTEST_EXE} catch_main einspline qmcutil)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
