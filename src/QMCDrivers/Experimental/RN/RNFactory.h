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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
