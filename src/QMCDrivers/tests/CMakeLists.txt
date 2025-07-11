#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(PLATFORM_PREFIX "")
if(ENABLE_CUDA)
    set(PLATFORM_PREFIX "nv_gpu")
else()
  set(PLATFORM_PREFIX "cpu")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/reference_prefix.h.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/reference_prefix.h
               @ONLY)
set(SRC_DIR drivers)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

# test_drivers shouldn't be dependent on wavefunctions
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)

set(DRIVER_TEST_SRC
    test_WalkerLogCollector.cpp
    test_TauParams.cpp
    test_vmc.cpp
    test_dmc.cpp
    test_drift.cpp
    test_clone_manager.cpp
    test_fixed_node_branch.cpp
    test_DescentEngine.cpp
    test_HybridEngine.cpp
    test_EngineHandle.cpp
    test_LMYEngine.cpp
    test_vmc_driver.cpp
    test_dmc_driver.cpp
    test_SimpleFixedNodeBranch.cpp
    test_QMCDriverFactory.cpp
    test_WalkerLogManager.cpp)

# legacy driver mpi, not tested in a meaningful way since the standard unit test is
# just one rank with mpi.
if(HAVE_MPI)
  set(DRIVER_TEST_SRC ${DRIVER_TEST_SRC} test_walker_control.cpp)
endif()
add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})

target_link_libraries(${UTEST_EXE} catch_main qmcdriver minimal_test_pools)
target_include_directories(${UTEST_EXE} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

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

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

# New Driver Only test program
# to avoid dealing with the build dependencies
# and memory leaks brought in by the old driver tests

set(UTEST_EXE test_new_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_new_${SRC_DIR})
#this is dependent on the directory creation and sym linking of earlier driver tests

set(DRIVER_TEST_SRC
    SetupPools.cpp
    test_Crowd.cpp
    test_MCPopulation.cpp
    test_QMCDriverInput.cpp
    test_QMCDriverNew.cpp
    test_VMCDriverInput.cpp
    test_VMCBatched.cpp
    test_DMCBatched.cpp
    test_SFNBranch.cpp
    test_QMCCostFunctionBatched.cpp
    test_QMCCostFunctionBase.cpp
    test_LinearMethod.cpp
    test_ConjugateGradient.cpp
    test_Eigensolver.cpp
    test_WFOptDriverInput.cpp)
add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})
target_link_libraries(${UTEST_EXE} catch_main qmcdriver minimal_test_pools)
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

add_unit_test(${UTEST_NAME} 1 8 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

if(HAVE_MPI)
  set(UTEST_EXE test_new_${SRC_DIR}_mpi)
  #this is dependent on the directory creation and sym linking of earlier driver tests
  set(DRIVER_TEST_SRC SetupPools.cpp test_WalkerControl.cpp test_QMCDriverNew.cpp)
  add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})
  target_link_libraries(${UTEST_EXE} catch_main qmcdriver mock_walker_elements)
  if(USE_OBJECT_TARGET)
    target_link_libraries(
      ${UTEST_EXE}
      qmcestimators
      qmcham
      qmcwfs
      qmcparticle
      qmcwfs_omptarget
      qmcparticle_omptarget
      qmcutil
      platform_omptarget_LA)
  endif()
  foreach(NUM_RANKS 1 2 3 4 16)
    set(UTEST_NAME deterministic-unit_${UTEST_EXE}-r${NUM_RANKS})
    add_unit_test(${UTEST_NAME} ${NUM_RANKS} 1 $<TARGET_FILE:${UTEST_EXE}>)
    set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  endforeach()

  set(UTEST_EXE test_new_${SRC_DIR}_integration)
  #this is dependent on the directory creation and sym linking of earlier driver tests
  set(DRIVER_TEST_SRC test_DMCBatched_integration.cpp)
  add_executable(${UTEST_EXE} ${DRIVER_TEST_SRC})
  set_property(
    TARGET ${UTEST_EXE}
    APPEND
    PROPERTY COMPILE_DEFINITIONS REFERENCE_PREFIX=${PLATFORM_PREFIX})
  target_link_libraries(${UTEST_EXE} catch_main qmcdriver mock_walker_elements)
  target_include_directories(${UTEST_EXE} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
  if(USE_OBJECT_TARGET)
    target_link_libraries(
      ${UTEST_EXE}
      qmcestimators
      qmcham
      qmcwfs
      qmcparticle
      qmcwfs_omptarget
      qmcparticle_omptarget
      qmcutil
      platform_omptarget_LA)
  endif()

  foreach(NUM_RANKS 3)
    set(UTEST_NAME deterministic-unit_${UTEST_EXE}-r${NUM_RANKS}-c4)
    add_unit_test(${UTEST_NAME} ${NUM_RANKS} 4 $<TARGET_FILE:${UTEST_EXE}>)
    set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
    set(UTEST_HDF_RANDOM_SEED_FILE rank_count_${NUM_RANKS}_seeds_for_DMCBatched.random.h5)
    set(UTEST_HDF_RANDOM_SEED ${CMAKE_SOURCE_DIR}/src/QMCDrivers/tests/reference/${UTEST_HDF_RANDOM_SEED_FILE})
    maybe_symlink(${UTEST_HDF_RANDOM_SEED} ${UTEST_DIR}/${UTEST_HDF_RANDOM_SEED_FILE})
    #yes the cmake devs never read https://www.cs.utexas.edu/~EWD/ewd08xx/EWD831.PDF
    math(EXPR RANK_LIMIT "${NUM_RANKS} - 1" OUTPUT_FORMAT DECIMAL)
    maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)
    foreach(REF_RANK RANGE ${RANK_LIMIT})
      set(UTEST_PPHL_FILE ${PLATFORM_PREFIX}_rank_${REF_RANK}of${NUM_RANKS}_dmc_vem_test.dat)
      set(UTEST_PPHL_REF ${CMAKE_SOURCE_DIR}/src/QMCDrivers/tests/reference/${UTEST_PPHL_FILE})
      maybe_symlink(${UTEST_PPHL_REF} ${UTEST_DIR}/${UTEST_PPHL_FILE})
    endforeach()
  endforeach()

endif(HAVE_MPI)
