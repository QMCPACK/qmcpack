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
/**@file DMCMoveAll.h
 * @brief Declaration of DMCMoveAll
 */
#ifndef QMCPLUSPLUS_DMC_MOVEALL_H
#define QMCPLUSPLUS_DMC_MOVEALL_H

#include "QMCDrivers/QMCDriver.h" 
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  class DMCUpdateBase;

  /** @ingroup QMCDrivers  WalkerByWalker
   *@brief implements the DMC algorithm using walker-by-walker move. 
   */
  class DMCMoveAll: public QMCDriver {

  public:
    /// Constructor.
    DMCMoveAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
    ///destructor
    ~DMCMoveAll();

    bool run();
    bool put(xmlNodePtr q);

  private:

    DMCUpdateBase *Mover;
    IndexType KillNodeCrossing;
    IndexType BranchInterval;
    ///Interval between branching
    IndexType NonLocalMoveIndex;
    string Reconfiguration;
    string KillWalker;
    ///input string to determine to use nonlocal move
    string NonLocalMove;
    /// Copy Constructor (disabled)
    DMCMoveAll(const DMCMoveAll& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    DMCMoveAll& operator=(const DMCMoveAll&) { return *this;}

    bool dmcWithBranching();
    bool dmcWithReconfiguration();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
