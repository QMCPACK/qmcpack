#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR drivers)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# test_drivers shouldn't be dependent on wavefunctions
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/drivers)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)

SET(DRIVER_TEST_SRC test_vmc.cpp test_dmc.cpp test_drift.cpp test_clone_manager.cpp test_fixed_node_branch.cpp test_DescentEngine.cpp test_HybridEngine.cpp)

IF (NOT QMC_CUDA)
  SET(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_vmc_driver.cpp test_dmc_driver.cpp)
ENDIF()

# legacy driver mpi, not tested in a meaningful way since the standard unit test is
# just one rank with mpi.
IF(HAVE_MPI)
  SET(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_walker_control.cpp)
ENDIF()
ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})

IF(QMC_CUDA)
  target_sources(${UTEST_EXE} PRIVATE test_QMCDriverFactory_CUDA.cpp)
ELSE(QMC_CUDA)
  target_sources(${UTEST_EXE} PRIVATE test_QMCDriverFactory.cpp)
ENDIF(QMC_CUDA)

USE_FAKE_RNG(${UTEST_EXE})
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcdriver_unit qmcfakerng)
IF(USE_OBJECT_TARGET)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcestimators_unit qmcham_unit qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
ENDIF()


IF(NOT QMC_CUDA)
  # New Driver Only test program
  # to avoid dealing with the build dependencies
  # and memory leaks brought in by the old driver tests

  SET(UTEST_EXE test_new_${SRC_DIR})
  SET(UTEST_NAME deterministic-unit_test_new_${SRC_DIR})
  SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/drivers)
  SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
  #this is dependent on the directory creation and sym linking of earlier driver tests

  SET(DRIVER_TEST_SRC SetupPools.cpp test_Crowd.cpp test_MCPopulation.cpp test_ContextForSteps.cpp test_QMCDriverInput.cpp test_QMCDriverNew.cpp test_VMCDriverInput.cpp test_VMCFactoryNew.cpp test_VMCBatched.cpp test_DMCBatched.cpp test_SimpleFixedNodeBranch.cpp test_SFNBranch.cpp test_QMCCostFunctionBatched.cpp)
  ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})
  USE_FAKE_RNG(${UTEST_EXE})
  TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcdriver_unit qmcfakerng)
  IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcestimators_unit qmcham_unit qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget)
  ENDIF()

  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
    SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
  ENDIF()


  IF(HAVE_MPI)
    SET(UTEST_EXE test_${SRC_DIR}_mpi)
    SET(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)
    SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/drivers)
    SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
    #this is dependent on the directory creation and sym linking of earlier driver tests
    SET(DRIVER_TEST_SRC SetupPools.cpp test_WalkerControl.cpp)
    ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})
    USE_FAKE_RNG(${UTEST_EXE})
    TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcdriver_unit qmcfakerng)
    IF(USE_OBJECT_TARGET)
    TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcestimators_unit qmcham_unit qmcwfs qmcparticle qmcutil containers platform_device platform_omptarget)
    ENDIF()
    # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
    ADD_UNIT_TEST(${UTEST_NAME} 3 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
    SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  ENDIF(HAVE_MPI)
ENDIF(NOT QMC_CUDA)
