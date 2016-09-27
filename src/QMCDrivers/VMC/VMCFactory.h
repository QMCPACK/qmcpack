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

struct VMCFactory
{
  int VMCMode;
  xmlNodePtr myNode;
  VMCFactory(int vmode, xmlNodePtr cur):VMCMode(vmode), myNode(cur) {}
  QMCDriver* create(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                    ParticleSetPool& ptclpool, HamiltonianPool& hpool, WaveFunctionPool& ppool);
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCFactory.h,v $   $Author: jnkim $
 * $Revision: 1.2 $   $Date: 2006/04/05 00:49:59 $
 * $Id: DMCFactory.h,v 1.2 2006/04/05 00:49:59 jnkim Exp $
 ***************************************************************************/
