//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VMCFACTORYNEW_H
#define QMCPLUSPLUS_VMCFACTORYNEW_H
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCApp/WaveFunctionPool.h"
#include "Message/Communicate.h"


namespace qmcplusplus
{
class ParticleSetPool;
class HamiltonianPool;
class MCPopulation;

class VMCFactoryNew
{
private:
  const int vmc_mode_;
  // const ?
  xmlNodePtr input_node_;
  const int qmc_counter_;

public:
  VMCFactoryNew(xmlNodePtr cur, const int vmode, const int qmc_counter)
      : vmc_mode_(vmode), input_node_(cur), qmc_counter_(qmc_counter)
  {}

  QMCDriverInterface* create(MCPopulation&& pop,
                             TrialWaveFunction& psi,
                             QMCHamiltonian& h,
                             ParticleSetPool& ptclpool,
                             HamiltonianPool& hpool,
                             WaveFunctionPool& ppool,
                             Communicate* comm);
};
} // namespace qmcplusplus

#endif
