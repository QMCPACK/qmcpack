//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DMC_FACTORY_H
#define QMCPLUSPLUS_DMC_FACTORY_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCApp/HamiltonianPool.h"

namespace qmcplusplus
{
class DMCFactory
{
private:
  bool PbyPUpdate, GPU;
  xmlNodePtr myNode;

public:
  DMCFactory(bool pbyp, bool gpu, xmlNodePtr cur) : PbyPUpdate(pbyp), GPU(gpu), myNode(cur) {}

  QMCDriver* create(MCWalkerConfiguration& w,
                    TrialWaveFunction& psi,
                    QMCHamiltonian& h,
                    HamiltonianPool& hpool,
                    WaveFunctionPool& ppool,
                    Communicate* comm);
};
} // namespace qmcplusplus

#endif
