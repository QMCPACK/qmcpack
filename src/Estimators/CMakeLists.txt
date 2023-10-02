#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

####################################
# Estimators which accumulate observables
####################################

set(QMCEST_SRC
    InputSection.cpp
    CSEnergyEstimator.cpp
    LocalEnergyEstimator.cpp
    RMCLocalEnergyEstimator.cpp
    EstimatorInput.cpp
    EstimatorManagerInput.cpp
    ScalarEstimatorInputs.cpp
    SpinDensityInput.cpp
    MomentumDistributionInput.cpp
    EstimatorManagerBase.cpp
    EstimatorManagerNew.cpp
    EstimatorManagerCrowd.cpp
    CollectablesEstimator.cpp
    OperatorEstBase.cpp
    SpinDensityNew.cpp
    MomentumDistribution.cpp
    OneBodyDensityMatricesInput.cpp
    OneBodyDensityMatrices.cpp
    MagnetizationDensity.cpp
    MagnetizationDensityInput.cpp
    PerParticleHamiltonianLoggerInput.cpp
    PerParticleHamiltonianLogger.cpp)

####################################
# create libqmcestimators
####################################
if(USE_OBJECT_TARGET)
  add_library(qmcestimators OBJECT ${QMCEST_SRC})
else()
  add_library(qmcestimators ${QMCEST_SRC})
endif()

target_include_directories(qmcestimators PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcestimators PUBLIC containers qmcham qmcparticle qmcutil)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
