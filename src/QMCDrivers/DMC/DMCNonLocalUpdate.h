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
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus {

  class DMCNonLocalUpdate: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCNonLocalUpdate(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCNonLocalUpdate();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:
    /// Copy Constructor (disabled)
    DMCNonLocalUpdate(const DMCNonLocalUpdate& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCNonLocalUpdate& operator=(const DMCNonLocalUpdate&) { return *this;}

  };

  class DMCNonLocalUpdatePbyP: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCNonLocalUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCNonLocalUpdatePbyP();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:
    /// Copy Constructor (disabled)
    DMCNonLocalUpdatePbyP(const DMCNonLocalUpdatePbyP& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCNonLocalUpdatePbyP& operator=(const DMCNonLocalUpdatePbyP&) { return *this;}
    vector<NewTimer*> myTimers;

  };

  class DMCNonLocalUpdatePbyPFast: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCNonLocalUpdatePbyPFast(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCNonLocalUpdatePbyPFast();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:
    /// Copy Constructor (disabled)
    DMCNonLocalUpdatePbyPFast(const DMCNonLocalUpdatePbyP& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCNonLocalUpdatePbyPFast& operator=(const DMCNonLocalUpdatePbyP&) { return *this;}
    vector<NewTimer*> myTimers;

  };

}

#endif
/***************************************************************************
 * $RCSfile: DMCNonLocalUpdate.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
