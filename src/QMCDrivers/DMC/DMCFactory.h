//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_DMC_FACTORY_H
#define QMCPLUSPLUS_DMC_FACTORY_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCApp/HamiltonianPool.h"

namespace qmcplusplus
{
struct DMCFactory
{
  bool PbyPUpdate, GPU;
  xmlNodePtr myNode;
  DMCFactory(bool pbyp, bool gpu, xmlNodePtr cur) :
    PbyPUpdate(pbyp), myNode(cur), GPU(gpu)
  { }

  QMCDriver* create(MCWalkerConfiguration& w,
                    TrialWaveFunction& psi,
                    QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool);
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCFactory.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 **************************************************************************/
