//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "DMCFactoryNew.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "EstimatorInputDelegates.h"
#include "Concurrency/OpenMP.h"
#include "Estimators/EstimatorManagerNew.h"

namespace qmcplusplus
{
std::unique_ptr<QMCDriverInterface> DMCFactoryNew::create(const ProjectData& project_data,
                                                          const std::optional<EstimatorManagerInput> global_emi,
                                                          WalkerConfigurations& wc,
                                                          MCPopulation&& pop,
                                                          const ParticleSetPool::PoolType& pset_pool,
                                                          Communicate* comm)
{
  app_summary() << "\n========================================"
                   "\n  Reading DMC driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  DMCDriverInput dmcdriver_input;
  try
  {
    qmcdriver_input.readXML(input_node_);
    dmcdriver_input.readXML(input_node_);
  }
  catch (const std::exception& e)
  {
    throw UniformCommunicateError(e.what());
  }

  // This is done so that the application level input structures reflect the actual input to the code.
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

  // In my opinion unwrapping qmcdriver_input_ here and not in QMCDriver new doesn't make sense.
  auto estimator_manager_ =
      std::make_unique<EstimatorManagerNew>(comm,
                                            makeEstimatorManagerInput(global_emi,
                                                                      qmcdriver_input.get_estimator_manager_input()),
                                            pop.get_golden_hamiltonian(), pop.get_golden_electrons(),
                                            pset_pool, pop.get_golden_twf());

  auto qmc = std::make_unique<DMCBatched>(project_data, std::move(qmcdriver_input), estimator_manager,
                                          std::move(dmcdriver_input), wc, std::move(pop), comm);
  // This can probably be eliminated completely since we only support PbyP
  qmc->setUpdateMode(dmc_mode_ & 1);
  return qmc;
}
} // namespace qmcplusplus
