//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_EE_FACTORY_H
#define QMCPLUSPLUS_EE_FACTORY_H
#include "QMCDrivers/QMCDriver.h"

namespace qmcplusplus
{
class ParticleSetPool;
class HamiltonianPool;

struct EEFactory
{
  int VMCMode;
  xmlNodePtr myNode;
  EEFactory(int vmode, xmlNodePtr cur):VMCMode(vmode), myNode(cur) {}
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
