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


#ifndef QMCPLUSPLUS_VMC_FACTORY_H
#define QMCPLUSPLUS_VMC_FACTORY_H
#include "QMCDrivers/QMCDriver.h"

namespace qmcplusplus
{
class ParticleSetPool;
class HamiltonianPool;

class VMCFactory
{
private:
  unsigned long VMCMode;
  xmlNodePtr myNode;

public:
  VMCFactory(unsigned long vmode, xmlNodePtr cur) : VMCMode(vmode), myNode(cur) {}

  std::unique_ptr<QMCDriverInterface> create(const ProjectData& project_data,
                                             MCWalkerConfiguration& w,
                                             TrialWaveFunction& psi,
                                             QMCHamiltonian& h,
                                             Communicate* comm,
                                             bool enable_profiling);
};
} // namespace qmcplusplus

#endif
