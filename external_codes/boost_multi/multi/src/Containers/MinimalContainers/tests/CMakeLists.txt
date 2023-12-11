#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(UTEST_EXE test_containers_MinimalContainers)
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

set(SRCS test_ConstantSizeMatrix.cpp test_RecordArray.cpp)

add_executable(${UTEST_EXE} ${SRCS})

target_link_libraries(${UTEST_EXE} catch_main containers)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
