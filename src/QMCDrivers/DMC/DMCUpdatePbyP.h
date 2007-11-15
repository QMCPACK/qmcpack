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
#ifndef QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#define QMCPLUSPLUS_DMC_UPDATE_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus {

  class DMCUpdatePbyPWithRejection: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdatePbyPWithRejection(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdatePbyPWithRejection();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:

    /// Copy Constructor (disabled)
    DMCUpdatePbyPWithRejection(const DMCUpdatePbyPWithRejection& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCUpdatePbyPWithRejection& operator=(const DMCUpdatePbyPWithRejection&) { return *this;}

  };

  class DMCUpdatePbyPWithKill: public QMCUpdateBase {

  public:

    /// Constructor.
    DMCUpdatePbyPWithKill(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCUpdatePbyPWithKill();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

  private:

    /// Copy Constructor (disabled)
    DMCUpdatePbyPWithKill(const DMCUpdatePbyPWithKill& a): QMCUpdateBase(a){ }
    /// Copy operator (disabled).
    DMCUpdatePbyPWithKill& operator=(const DMCUpdatePbyPWithKill&) { return *this;}

  };
}

#endif
/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
