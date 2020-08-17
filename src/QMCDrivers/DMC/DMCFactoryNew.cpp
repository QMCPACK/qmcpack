//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/DMC/DMCFactoryNew.h"
#include "QMCDrivers/DMC/DMCBatched.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{
QMCDriverInterface* DMCFactoryNew::create(MCPopulation& pop,
                                          TrialWaveFunction& psi,
                                          QMCHamiltonian& h,
                                          WaveFunctionPool& wf_pool,
                                          Communicate* comm)
{
  QMCDriverInput qmcdriver_input(qmc_counter_);
  qmcdriver_input.readXML(input_node_);
  DMCDriverInput dmcdriver_input;
  dmcdriver_input.readXML(input_node_);
  QMCDriverInterface* qmc = new DMCBatched(std::move(qmcdriver_input), std::move(dmcdriver_input), pop, psi, h, wf_pool, comm);
  // This can probably be eliminated completely since we only support PbyP
  qmc->setUpdateMode(dmc_mode_ & 1);
  return qmc;
}
} // namespace qmcplusplus
