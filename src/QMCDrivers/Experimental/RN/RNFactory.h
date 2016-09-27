//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_RN_FACTORY_H
#define QMCPLUSPLUS_RN_FACTORY_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCApp/HamiltonianPool.h"

namespace qmcplusplus
{
struct RNFactory
{
  bool PbyPUpdate;
  xmlNodePtr myNode;
  RNFactory(bool pbyp, xmlNodePtr cur):PbyPUpdate(pbyp), myNode(cur) {}
  QMCDriver* create(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool);
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCFactory.h,v $   $Author: jnkim $
 * $Revision: 2435 $   $Date: 2007-12-07 08:22:22 -0600 (Fri, 07 Dec 2007) $
 * $Id: DMCFactory.h 2435 2007-12-07 14:22:22Z jnkim $
 ***************************************************************************/
