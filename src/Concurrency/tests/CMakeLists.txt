#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR Concurrency)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

set(SRCS test_ParallelExecutorOPENMP.cpp test_UtilityFunctionsOPENMP.cpp)

if(QMC_EXP_THREADING)
  set(SRCS ${SRCS} test_ParallelExecutorSTD.cpp)
endif(QMC_EXP_THREADING)
add_executable(${UTEST_EXE} ${SRCS})

target_link_libraries(${UTEST_EXE} catch_main)

add_unit_test(${UTEST_NAME} 1 3 $<TARGET_FILE:${UTEST_EXE}>)
