#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

# Also includes tests for code OhmmsApp directory

SET(SRC_DIR qmcapp)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

ADD_EXECUTABLE(${UTEST_EXE} test_project_data.cpp test_rng_control.cpp)

TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmc)
IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcestimators qmcdriver qmcham qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget)
  IF(BUILD_AFQMC)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} afqmc)
  ENDIF(BUILD_AFQMC)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")

# temporary label until asan leaks are resolved
IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
ENDIF()
