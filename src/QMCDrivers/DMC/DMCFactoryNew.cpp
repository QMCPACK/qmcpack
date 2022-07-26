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

namespace qmcplusplus
{
QMCDriverInterface* DMCFactoryNew::create(const ProjectData& project_data,
                                          const std::optional<EstimatorManagerInput> global_emi,
                                          WalkerConfigurations& wc,
                                          MCPopulation&& pop,
                                          Communicate* comm)
{
#if defined(QMC_CUDA)
  comm->barrier_and_abort("DMC batched driver is not supported by legacy CUDA builds.");
#endif

  app_summary() << "\n========================================"
                   "\n  Reading DMC driver XML input section"
                   "\n========================================"
                << std::endl;

  QMCDriverInput qmcdriver_input;
  qmcdriver_input.readXML(input_node_);
  DMCDriverInput dmcdriver_input;
  dmcdriver_input.readXML(input_node_);
  QMCDriverInterface* qmc = new DMCBatched(project_data, std::move(qmcdriver_input), global_emi,
                                           std::move(dmcdriver_input), wc, std::move(pop), comm);
  // This can probably be eliminated completely since we only support PbyP
  qmc->setUpdateMode(dmc_mode_ & 1);
  return qmc;
}
} // namespace qmcplusplus
