#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(UTEST_EXE test_omptarget)
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

add_executable(${UTEST_EXE} test_vector.cpp test_math.cpp test_deep_copy.cpp test_class_member.cpp test_runtime_mem.cpp)
target_link_libraries(${UTEST_EXE} platform_runtime catch_main)

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)

set(UTEST_EXE test_omptarget_blas)
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

add_executable(${UTEST_EXE} test_ompBLAS.cpp)
target_link_libraries(${UTEST_EXE} catch_main platform_LA containers)

if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} platform_omptarget_LA)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
