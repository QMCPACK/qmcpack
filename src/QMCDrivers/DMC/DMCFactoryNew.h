//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: DMCFactor.h
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMCFACTORYNEW_H
#define QMCPLUSPLUS_DMCFACTORYNEW_H
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCApp/WaveFunctionPool.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
class ParticleSetPool;
class HamiltonianPool;
class MCPopulation;
  
class DMCFactoryNew
{
private:
  const int dmc_mode_;
  xmlNodePtr input_node_;
  const int qmc_counter_;

  
public:
  DMCFactoryNew(xmlNodePtr cur, const int dmc_mode, const int qmc_counter)
      : dmc_mode_(dmc_mode), input_node_(cur), qmc_counter_(qmc_counter)
  {}

  QMCDriverInterface* create(MCPopulation& pop,
                             TrialWaveFunction& psi,
                             QMCHamiltonian& h,
                             WaveFunctionPool& wf_pool,
                             Communicate* comm);
};
} // namespace qmcplusplus

#endif
