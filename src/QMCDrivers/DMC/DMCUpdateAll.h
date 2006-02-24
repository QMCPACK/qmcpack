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
#ifndef QMCPLUSPLUS_DMC_ALLPARTICLE_UPDATE_H
#define QMCPLUSPLUS_DMC_ALLPARTICLE_UPDATE_H
#include "QMCDrivers/DMC/DMCUpdateBase.h"
namespace qmcplusplus {

  class DMCUpdateAllWithRejection: public DMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdateAllWithRejection(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdateAllWithRejection();

    void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
    void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);
    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);
  private:
    /// Copy Constructor (disabled)
    DMCUpdateAllWithRejection(const DMCUpdateAllWithRejection& a): DMCUpdateBase(a){}
    /// Copy operator (disabled).
    DMCUpdateAllWithRejection& operator=(const DMCUpdateAllWithRejection&) { return *this;}
  };

  class DMCUpdateAllWithKill: public DMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdateAllWithKill(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdateAllWithKill();

    void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
    void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);
    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);
  private:
    /// Copy Constructor (disabled)
    DMCUpdateAllWithKill(const DMCUpdateAllWithKill& a): DMCUpdateBase(a){}
    /// Copy operator (disabled).
    DMCUpdateAllWithKill& operator=(const DMCUpdateAllWithKill&) { return *this;}
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
