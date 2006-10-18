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
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_TESTING_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_TESTING_H
#include "QMCDrivers/QMCDriver.h" 
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  //class DMCPbyPUpdate;
  class DMCUpdateBase;

  /** @ingroup QMCDrivers ParticleByParticle 
   *@brief Implements the DMC algorithm using particle-by-particle move. 
   */
  class DMCPbyP: public QMCDriver {
  public:
    /// Constructor.
    DMCPbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    ///destructor
    ~DMCPbyP();

    bool run();
    bool put(xmlNodePtr cur);

  private:

    ///Index to determine what to do when node crossing is detected
    IndexType KillNodeCrossing;
    ///Total number of accepted steps per block
    IndexType nAcceptTot;
    ///Total number of rejected steps per block
    IndexType nRejectTot;
    ///Interval between branching
    IndexType BranchInterval;
    ///Interval between branching
    IndexType NonLocalMoveIndex;
    ///hdf5 file name for Branch conditions
    string BranchInfo;
    ///input string to determine kill walkers or not
    string KillWalker;
    ///input string to determine swap walkers among mpi processors
    string SwapWalkers;
    ///input string to determine to use reconfiguration
    string Reconfiguration;
    ///input string to determine to use nonlocal move
    string NonLocalMove;
    ///update engine
    DMCUpdateBase *Mover;
    /// Copy Constructor (disabled)
    DMCPbyP(const DMCPbyP& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    DMCPbyP& operator=(const DMCPbyP&) { return *this;}

    bool dmcWithReconfiguration();
    bool dmcWithBranching();

  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
