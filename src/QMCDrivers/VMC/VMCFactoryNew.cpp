//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "VMCFactoryNew.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "EstimatorInputDelegates.h"
#include "Concurrency/Info.hpp"
#include "Estimators/EstimatorManagerNew.h"

namespace qmcplusplus
{
std::unique_ptr<QMCDriverInterface> VMCFactoryNew::create(const ProjectData& project_data,
                                                          const std::optional<EstimatorManagerInput>& global_emi,
                                                          WalkerConfigurations& wc,
                                                          MCPopulation&& pop,
							  const ParticleSetPool::PoolType& pset_pool,
                                                          SampleStack& samples,
                                                          Communicate* comm)
{
  app_summary() << "\n========================================"
                   "\n  Reading VMC driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  VMCDriverInput vmcdriver_input;
  try
  {
    qmcdriver_input.readXML(input_node_);
    vmcdriver_input.readXML(input_node_);
  }
  catch (const std::exception& e)
  {
    throw UniformCommunicateError(e.what());
  }

  std::unique_ptr<QMCDriverInterface> qmc;

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
					    pset_pool,
                                            pop.get_golden_twf());


  if (vmc_mode_ == 0 || vmc_mode_ == 1) //(0,0,0) (0,0,1)
  {
    qmc = std::make_unique<VMCBatched>(project_data, std::move(qmcdriver_input), estimator_manager, std::move(vmcdriver_input),
                                       wc, std::move(pop), samples, comm);
  }
  else
  {
    throw std::runtime_error("VMCFactoryNew does not support VMC_MODE");
  }

  // TODO: I believe this is redundant and could become a bug.
  qmc->setUpdateMode(vmc_mode_ & 1);
  return qmc;
}
} // namespace qmcplusplus
