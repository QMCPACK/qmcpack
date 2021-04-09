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

SET(SRC_DIR estimators)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

SET(SRCS test_accumulator.cpp test_local_energy_est.cpp FakeOperatorEstimator.cpp EstimatorManagerBaseTest.cpp EstimatorManagerNewTest.cpp test_manager.cpp test_EstimatorManagerNew.cpp test_trace_manager.cpp SpinDensityTesting.cpp test_SpinDensityInput.cpp test_SpinDensityNew.cpp)

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcutil qmcestimators_unit utilities_for_test)
IF(USE_OBJECT_TARGET)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcutil qmcestimators_unit qmcham_unit qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget utilities_for_test)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
ENDIF()

IF(HAVE_MPI)
  SET(UTEST_EXE test_${SRC_DIR}_mpi)
  SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)
  #this is dependent on the directory creation and sym linking of earlier driver tests
  SET(SRCS FakeOperatorEstimator.cpp EstimatorManagerNewTest.cpp test_manager_mpi.cpp)
  ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
  IF(USE_OBJECT_TARGET)
    TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcestimators_unit qmcham_unit qmcdriver_unit qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget utilities_for_test)
  ENDIF()
  TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcestimators_unit)
  # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
  ADD_UNIT_TEST(${UTEST_NAME} 3 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ENDIF()
