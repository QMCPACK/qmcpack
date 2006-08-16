//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DMC_NONLOCAL_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_NONLOCAL_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/DMC/DMCUpdateBase.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
namespace qmcplusplus {

  class DMCNonLocalUpdate: public DMCUpdateBase {

  public:

    /// Constructor.
    DMCNonLocalUpdate(ParticleSet& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCNonLocalUpdate();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);

    bool put(xmlNodePtr cur);

  private:

    NonLocalTOperator nonLocalOps;

    /// Copy Constructor (disabled)
    DMCNonLocalUpdate(const DMCNonLocalUpdate& a): DMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCNonLocalUpdate& operator=(const DMCNonLocalUpdate&) { return *this;}

  };

  class DMCNonLocalUpdatePbyP: public DMCUpdateBase {

  public:

    /// Constructor.
    DMCNonLocalUpdatePbyP(ParticleSet& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCNonLocalUpdatePbyP();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);

    bool put(xmlNodePtr cur);

  private:

    NonLocalTOperator nonLocalOps;
    /// Copy Constructor (disabled)
    DMCNonLocalUpdatePbyP(const DMCNonLocalUpdatePbyP& a): DMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCNonLocalUpdatePbyP& operator=(const DMCNonLocalUpdatePbyP&) { return *this;}

  };

}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
