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

set(QMCEST_SRC
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
if(USE_OBJECT_TARGET)
  add_library(qmcestimators OBJECT ${QMCEST_SRC})
  add_library(qmcestimators_unit OBJECT ${QMCEST_SRC})
else()
  add_library(qmcestimators ${QMCEST_SRC})
  add_library(qmcestimators_unit ${QMCEST_SRC})
endif()
use_fake_rng(qmcestimators_unit)

target_include_directories(qmcestimators PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(qmcestimators_unit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcestimators PUBLIC containers qmcham qmcparticle qmcutil)
target_link_libraries(qmcestimators_unit PUBLIC containers qmcham_unit qmcparticle qmcutil)

add_subdirectory(tests)

# TARGET_LINK_LIBRARIES(qmcest PUBLIC qmcham)
# TARGET_LINK_LIBRARIES(qmcdriver_unit PUBLIC qmcham_unit)
# TARGET_LINK_LIBRARIES(qmcdriver PRIVATE Math::BLAS_LAPACK)
# TARGET_LINK_LIBRARIES(qmcdriver_unit PRIVATE Math::BLAS_LAPACK)
# IF(BUILD_LMYENGINE_INTERFACE)
#   TARGET_LINK_LIBRARIES(qmcdriver PRIVATE formic_utils)
#   TARGET_LINK_LIBRARIES(qmcdriver_unit PRIVATE formic_utils)
# ENDIF()
