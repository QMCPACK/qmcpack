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
#ifndef OHMMS_QMC_REPTATION_H
#define OHMMS_QMC_REPTATION_H

#include "QMC/QMCDriver.h" 
#include <deque>
namespace ohmmsqmc {

  /** Implements the DMC algorithm. */
  class ReptationMC: public QMCDriver {

  public:

    /// Constructor.
    ReptationMC(MCWalkerConfiguration& w, 
	       TrialWaveFunction& psi, 
	       QMCHamiltonian& h, 
	       xmlNodePtr q);

    bool run();
    bool put(xmlNodePtr q);

  protected:

    ///boolean for the direction. true, if head is moving
    bool MoveHead;

    ///boolean for using bounce algorithm. true, if bounce algorithm of D. Ceperley
    bool UseBounce;

    ///the number of the beads that will be cut
    int  NumBeads;

    ///paths
    std::deque<MCWalkerConfiguration::Walker_t*> Polymer;

    std::vector<MCWalkerConfiguration::Walker_t*> WalkerRepository;

  private:

    /// Copy Constructor (disabled)
    ReptationMC(const ReptationMC& a): QMCDriver(a) { }

    /// Copy operator (disabled).
    ReptationMC& operator=(const ReptationMC&) { return *this;}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
