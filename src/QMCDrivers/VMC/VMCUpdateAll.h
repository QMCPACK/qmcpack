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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_VMC_UPDATEALL_H
#define QMCPLUSPLUS_VMC_UPDATEALL_H
#include "QMCDrivers/QMCUpdateBase.h"

namespace qmcplusplus {

  /** @ingroup QMCDrivers  ParticleByParticle
   *@brief Implements the VMC algorithm using particle-by-particle move. 
   */
  class VMCUpdateAll: public QMCUpdateBase {
  public:
    /// Constructor.
    VMCUpdateAll(ParticleSet& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);

    ~VMCUpdateAll();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);
 
    bool put(xmlNodePtr cur);

  private:
    ///sub steps
    int nSubSteps;
    /// Copy Constructor (disabled)
    VMCUpdateAll(const VMCUpdateAll& a): QMCUpdateBase(a) { }
    /// Copy operator (disabled).
    VMCUpdateAll& operator=(const VMCUpdateAll&) { return *this;}
  };

  /** @ingroup QMCDrivers  ParticleByParticle
   *@brief Implements the VMC algorithm using particle-by-particle move with the drift equation. 
   */
  class VMCUpdateAllWithDrift: public QMCUpdateBase {
  public:
    /// Constructor.
    VMCUpdateAllWithDrift(ParticleSet& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);

    ~VMCUpdateAllWithDrift();

    void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end);
 
  private:
    /// Copy Constructor (disabled)
    VMCUpdateAllWithDrift(const VMCUpdateAllWithDrift& a): QMCUpdateBase(a) { }
    /// Copy operator (disabled).
    VMCUpdateAllWithDrift& operator=(const VMCUpdateAllWithDrift&) { return *this;}
  };
}

#endif
/***************************************************************************
 * $RCSfile: VMCUpdateAll.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCUpdateAll.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $ 
 ***************************************************************************/
