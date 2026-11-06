#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2026 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR qmc_app)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(${UTEST_EXE} test_qmc_main.cpp)
target_link_libraries(${UTEST_EXE} catch_main qmc qmcdriver minimal_test_pools)
if(USE_OBJECT_TARGET)
  target_link_libraries(
    ${UTEST_EXE}
    qmcestimators
    qmcham
    qmcwfs
    qmcwfs_omptarget
    spline2
    spline2_omptarget
    qmcparticle
    qmcparticle_omptarget
    qmcutil
    platform_omptarget_LA)
endif()

add_unit_test(BASENAME ${UTEST_NAME} TEST_BINARY $<TARGET_FILE:${UTEST_EXE}>)
add_unit_test(BASENAME ${UTEST_NAME} PROCS 2 THREADS 3 TEST_BINARY $<TARGET_FILE:${UTEST_EXE}>)

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
