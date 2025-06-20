//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "MakeEstimatorManager.h"
#include "QMCHamiltonian.h"
#include "TrialWaveFunction.h"
#include <optional>
#include <memory>
#include <Communicate.h>
#include "Estimators/EstimatorInputDelegates.h"

namespace qmcplusplus
{

UPtr<EstimatorManagerNew> makeEstimatorManager(const std::optional<EstimatorManagerInput>& global_emi,
                                               const std::optional<EstimatorManagerInput>& driver_emi,
                                               ParticleSet& pset_primary,
                                               TrialWaveFunction& twf_primary,
                                               QMCHamiltonian& ham_primary,
                                               const EstimatorManagerNew::PSPool& pset_pool,
                                               Communicate* comm)
{ // This is done so that the application level input structures reflect the actual input to the code.
  // While the actual simulation objects still take singular input structures at construction.
  auto make_estimator_manager_input = [](auto& global_emi, auto& local_emi) -> EstimatorManagerInput {
    if (global_emi.has_value() && local_emi.has_value())
      return {global_emi.value(), local_emi.value()};
    else if (global_emi.has_value())
      return {global_emi.value()};
    else if (local_emi.has_value())
      return {local_emi.value()};
    else
      return {};
  };

  auto estimator_manager = std::make_unique<EstimatorManagerNew>(ham_primary, comm);
  estimator_manager->constructEstimators(make_estimator_manager_input(global_emi, driver_emi), pset_primary,
                                         twf_primary, ham_primary, pset_pool);
  return estimator_manager;
}

} // namespace qmcplusplus
