//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TEST_ESTIMATOR_MANAGER_INPUT
#define QMCPLUSPLUS_TEST_ESTIMATOR_MANAGER_INPUT

#include "EstimatorInputDelegates.h"
#include "StructureFactorInput.h"

namespace qmcplusplus::testing
{

template<typename T>
struct ExpectedEstimatorInputName;

template<>
struct ExpectedEstimatorInputName<EnergyDensityInput>
{
  using Type = EnergyDensityInput;
  std::string name{"EDcell"};
};

template<>
struct ExpectedEstimatorInputName<OneBodyDensityMatricesInput>
{
  using Type = OneBodyDensityMatricesInput;
  std::string name{"OneBodyDensityMatrices"};
};

template<>
struct ExpectedEstimatorInputName<MomentumDistributionInput>
{
  using Type = MomentumDistributionInput;
  std::string name{"nofk"};
};

template<>
struct ExpectedEstimatorInputName<StructureFactorInput>
{
  using Type = StructureFactorInput;
  std::string name{"sk1"};
};

} // namespace qmcplusplus::testing

#endif
