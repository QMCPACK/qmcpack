#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(SRC_DIR estimators)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})
set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)

set(SRCS
    test_accumulator.cpp
    test_local_energy_est.cpp
    FakeOperatorEstimator.cpp
    EstimatorManagerBaseTest.cpp
    EstimatorManagerNewTest.cpp
    test_manager.cpp
    test_EstimatorManagerNew.cpp
    test_trace_manager.cpp
    EstimatorTesting.cpp
    test_SpinDensityInput.cpp
    test_SpinDensityNew.cpp
    test_InputSection.cpp
    test_EstimatorManagerInput.cpp
    )

# OneBodyDensityMatrices are not tested with legacy CUDA
# DiracDeterminantsCUDA cannot be copied
if(NOT QMC_CUDA)
  set(SRCS ${SRCS} test_MomentumDistribution.cpp test_OneBodyDensityMatricesInput.cpp test_OneBodyDensityMatrices.cpp)
endif()

add_executable(${UTEST_EXE} ${SRCS})
use_fake_rng(${UTEST_EXE})
target_link_libraries(${UTEST_EXE} catch_main qmcutil qmcestimators_unit utilities_for_test)
if(USE_OBJECT_TARGET)
  target_link_libraries(
    ${UTEST_EXE}
    qmcutil
    qmcestimators_unit
    qmcham_unit
    qmcwfs
    qmcparticle
    qmcutil
    platform_omptarget_LA
    utilities_for_test)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)

if(HAVE_MPI)
  set(UTEST_EXE test_${SRC_DIR}_mpi)
  set(UTEST_NAME deterministic-unit_test_${SRC_DIR}_mpi)

  #this is dependent on the directory creation and sym linking of earlier driver tests
  set(SRCS FakeOperatorEstimator.cpp EstimatorManagerNewTest.cpp test_manager_mpi.cpp)
  add_executable(${UTEST_EXE} ${SRCS})
  if(USE_OBJECT_TARGET)
    target_link_libraries(
      ${UTEST_EXE}
      qmcestimators_unit
      qmcham_unit
      qmcdriver_unit
      qmcwfs
      qmcparticle
      qmcutil
      platform_omptarget_LA
      utilities_for_test)
  endif()
  target_link_libraries(${UTEST_EXE} catch_main qmcestimators_unit)
  # Right now the unified driver mpi tests are hard coded for 3 MPI ranks
  add_unit_test(${UTEST_NAME} 3 1 $<TARGET_FILE:${UTEST_EXE}>)
endif()
