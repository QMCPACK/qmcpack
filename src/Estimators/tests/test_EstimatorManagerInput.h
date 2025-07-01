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
struct ExpectedEstimatorInputNameType;

template<>
struct ExpectedEstimatorInputNameType<EnergyDensityInput>
{
  using Type = EnergyDensityInput;
  std::string name{"EDcell"};
  std::string type{"EnergyDensity"};
};

template<>
struct ExpectedEstimatorInputNameType<OneBodyDensityMatricesInput>
{
  using Type = OneBodyDensityMatricesInput;
  std::string name{"OneBodyDensityMatrices"};
  std::string type{"OneBodyDensityMatrices"};
};

template<>
struct ExpectedEstimatorInputNameType<MomentumDistributionInput>
{
  using Type = MomentumDistributionInput;
  std::string name{"nofk"};
  std::string type{"MomentumDistribution"};
};

template<>
struct ExpectedEstimatorInputNameType<StructureFactorInput>
{
  using Type = StructureFactorInput;
  std::string name{"sk1"};
  std::string type{"StructureFactor"};
};

} // namespace qmcplusplus::testing

#endif
