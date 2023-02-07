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

namespace qmcplusplus
{
std::unique_ptr<QMCDriverInterface> VMCFactoryNew::create(const ProjectData& project_data,
                                                          const std::optional<EstimatorManagerInput>& global_emi,
                                                          WalkerConfigurations& wc,
                                                          MCPopulation&& pop,
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

  if (vmc_mode_ == 0 || vmc_mode_ == 1) //(0,0,0) (0,0,1)
  {
    qmc = std::make_unique<VMCBatched>(project_data, std::move(qmcdriver_input), global_emi, std::move(vmcdriver_input),
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
