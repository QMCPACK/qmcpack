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
#include <optional>

namespace qmcplusplus
{

UPtr<EstimatorManagerNew> makeEstimatorManager(const std::optional<EstimatorManagerInput>& global_emi,
                                               const std::optional<EstimatorManagerInput>& driver_emi)
{ // This is done so that the application level input structures reflect the actual input to the code.
  // While the actual simulation objects still take singular input structures at construction.
  auto makeEstimatorManagerInput = [](auto& global_emi, auto& local_emi) -> EstimatorManagerInput {
    if (global_emi.has_value() && local_emi.has_value())
      return {global_emi.value(), local_emi.value()};
    else if (global_emi.has_value())
      return {global_emi.value()};
    else if (local_emi.has_value())
      return {local_emi.value()};
    else
      return {};
  };

  auto estimator_manager = std::make_unique<EstimatorManagerNew>(*primaryH, comm);
  estimator_manager->constructEstimators(makeEstimatorManagerInput(global_emi, driver_emi), qmc_system, *primaryPsi,
                                         *primaryH, particle_pool.getPool());
  return estimator_manager;
}

} // namespace qmcplusplus
