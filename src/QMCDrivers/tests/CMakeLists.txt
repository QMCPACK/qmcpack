#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
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

IF(HAVE_MPI)
  SET(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_walker_control.cpp)
ENDIF()
ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})

USE_FAKE_RNG(${UTEST_EXE})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcdriver_unit qmcham_unit qmcwfs qmcbase qmcutil qmcfakerng ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})


IF(NOT QMC_CUDA AND ENABLE_SOA)
# New Driver Only test program
# to avoid dealing with the build dependencies
# and memory leaks brought in by the old driver tests

SET(UTEST_EXE test_new_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_new_${SRC_DIR})

# test_drivers shouldn't be dependent on tests/wavefunctions 
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/drivers)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)

SET(DRIVER_TEST_SRC SetupPools.cpp test_Crowd.cpp test_MCPopulation.cpp test_ContextForSteps.cpp test_QMCDriverInput.cpp test_QMCDriverNew.cpp test_VMCDriverInput.cpp test_VMCFactoryNew.cpp test_VMCBatched.cpp test_DMCBatched.cpp)

IF(HAVE_MPI)
  SET(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_WalkerControlMPI.cpp)
ENDIF()

ADD_EXECUTABLE(${UTEST_EXE} ${DRIVER_TEST_SRC})

USE_FAKE_RNG(${UTEST_EXE})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcdriver_unit qmcham_unit qmcwfs qmcbase qmcutil qmcfakerng ${QMC_UTIL_LIBS} ${MPI_LIBRARY})


ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
ENDIF(NOT QMC_CUDA AND ENABLE_SOA)
