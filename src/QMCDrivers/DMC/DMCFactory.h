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
#ifndef QMCPLUSPLUS_DMC_FACTORY_H
#define QMCPLUSPLUS_DMC_FACTORY_H
#include "QMCDrivers/QMCDriver.h" 

namespace qmcplusplus {

  class HamiltonianPool;
  class RandomNumberControl;

  struct DMCFactory {
    bool PbyPUpdate;
    xmlNodePtr myNode;
    DMCFactory(bool pbyp, xmlNodePtr cur):PbyPUpdate(pbyp), myNode(cur){}
    QMCDriver* create(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, 
        RandomNumberControl& rc, HamiltonianPool& hpool);
  };
}

#endif
/***************************************************************************
 * $RCSfile: DMCFactory.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
