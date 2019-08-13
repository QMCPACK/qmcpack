//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
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
  int VMCMode;
  xmlNodePtr myNode;

public:
  VMCFactoryNew(int vmode, xmlNodePtr cur) : VMCMode(vmode), myNode(cur) {}

  QMCDriverInterface* create(MCPopulation& pop,
                             TrialWaveFunction& psi,
                             QMCHamiltonian& h,
                             ParticleSetPool& ptclpool,
                             HamiltonianPool& hpool,
                             WaveFunctionPool& ppool,
                             Communicate* comm);
};
} // namespace qmcplusplus

#endif
