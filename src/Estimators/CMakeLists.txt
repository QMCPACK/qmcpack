#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////
    
####################################
# Estimators which accumulate observables
####################################

SET(QMCEST_SRC
  CSEnergyEstimator.cpp
  LocalEnergyEstimator.cpp
  RMCLocalEnergyEstimator.cpp
  SpinDensityInput.cpp
  EstimatorManagerBase.cpp
  EstimatorManagerNew.cpp
  EstimatorManagerCrowd.cpp
  CollectablesEstimator.cpp
  OperatorEstBase.cpp
  SpinDensityNew.cpp)

####################################
# create libqmcestimators
####################################
IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcestimators OBJECT ${QMCEST_SRC})
  ADD_LIBRARY(qmcestimators_unit OBJECT ${QMCEST_SRC})
ELSE()
  ADD_LIBRARY(qmcestimators ${QMCEST_SRC})
  ADD_LIBRARY(qmcestimators_unit ${QMCEST_SRC})
ENDIF()
USE_FAKE_RNG(qmcestimators_unit)

TARGET_INCLUDE_DIRECTORIES(qmcestimators PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_INCLUDE_DIRECTORIES(qmcestimators_unit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_LINK_LIBRARIES(qmcestimators PUBLIC containers qmcham qmcparticle qmcutil)
TARGET_LINK_LIBRARIES(qmcestimators_unit PUBLIC containers qmcham_unit qmcparticle qmcutil)

add_subdirectory(tests)

# TARGET_LINK_LIBRARIES(qmcest PUBLIC qmcham)
# TARGET_LINK_LIBRARIES(qmcdriver_unit PUBLIC qmcham_unit)
# TARGET_LINK_LIBRARIES(qmcdriver PRIVATE Math::BLAS_LAPACK)
# TARGET_LINK_LIBRARIES(qmcdriver_unit PRIVATE Math::BLAS_LAPACK)
# IF(BUILD_LMYENGINE_INTERFACE)
#   TARGET_LINK_LIBRARIES(qmcdriver PRIVATE formic_utils)
#   TARGET_LINK_LIBRARIES(qmcdriver_unit PRIVATE formic_utils)
# ENDIF()
