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

set(SRC_DIR drivers)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# test_drivers shouldn't be dependent on wavefunctions
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)

set(DRIVER_TEST_SRC
    test_TauParams.cpp
    test_vmc.cpp
    test_dmc.cpp
    test_drift.cpp
    test_clone_manager.cpp
    test_fixed_node_branch.cpp
    test_DescentEngine.cpp
    test_HybridEngine.cpp)

if(NOT QMC_CUDA)
  set(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_vmc_driver.cpp test_dmc_driver.cpp)
endif()

# legacy driver mpi, not tested in a meaningful way since the standard unit test is
# just one rank with mpi.
if(HAVE_MPI)
  set(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_walker_control.cpp)
endif()
add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})

if(QMC_CUDA)
  target_sources(${UTEST_EXE} PRIVATE test_QMCDriverFactory_CUDA.cpp)
else(QMC_CUDA)
  target_sources(${UTEST_EXE} PRIVATE test_QMCDriverFactory.cpp)
endif(QMC_CUDA)

use_fake_rng(${UTEST_EXE})
target_link_libraries(${UTEST_EXE} catch_main qmcdriver_unit)
if(USE_OBJECT_TARGET)
  target_link_libraries(
    ${UTEST_EXE}
    qmcestimators_unit
    qmcham_unit
    qmcwfs
    qmcparticle
    qmcutil
    platform_omptarget_LA)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

if(NOT QMC_CUDA)
  # New Driver Only test program
  # to avoid dealing with the build dependencies
  # and memory leaks brought in by the old driver tests

  set(UTEST_EXE test_new_${SRC_DIR})
  set(UTEST_NAME deterministic-unit_test_new_${SRC_DIR})
  set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
  #this is dependent on the directory creation and sym linking of earlier driver tests

  set(DRIVER_TEST_SRC
      SetupPools.cpp
      test_Crowd.cpp
      test_MCPopulation.cpp
      test_ContextForSteps.cpp
      test_QMCDriverInput.cpp
      test_QMCDriverNew.cpp
      test_VMCDriverInput.cpp
      test_VMCFactoryNew.cpp
      test_VMCBatched.cpp
      test_DMCBatched.cpp
      test_SimpleFixedNodeBranch.cpp
      test_SFNBranch.cpp
      test_QMCCostFunctionBatched.cpp
      test_QMCCostFunctionBase.cpp
      test_WFOptDriverInput.cpp)
  add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})
  use_fake_rng(${UTEST_EXE})
  target_link_libraries(${UTEST_EXE} catch_main qmcdriver_unit)
  if(USE_OBJECT_TARGET)
    target_link_libraries(
      ${UTEST_EXE}
      qmcestimators_unit
      qmcham_unit
      qmcwfs
      qmcparticle
      qmcutil
      platform_omptarget_LA)
  endif()

  add_unit_test(${UTEST_NAME} 1 8 $<TARGET_FILE:${UTEST_EXE}>)
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

  if(HAVE_MPI)
    set(UTEST_EXE test_${SRC_DIR}_mpi)
    set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)
    set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
    #this is dependent on the directory creation and sym linking of earlier driver tests
    set(DRIVER_TEST_SRC SetupPools.cpp test_WalkerControl.cpp)
    add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})
    use_fake_rng(${UTEST_EXE})
    target_link_libraries(${UTEST_EXE} catch_main qmcdriver_unit)
    if(USE_OBJECT_TARGET)
      target_link_libraries(
        ${UTEST_EXE}
        qmcestimators_unit
        qmcham_unit
        qmcwfs
        qmcparticle
        qmcutil
        platform_omptarget_LA)
    endif()
    # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
    add_unit_test(${UTEST_NAME} 3 1 $<TARGET_FILE:${UTEST_EXE}>)
    set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  endif(HAVE_MPI)
endif(NOT QMC_CUDA)
