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
#ifndef OHMMS_QMC_VMC_H
#define OHMMS_QMC_VMC_H
#include "QMC/QMCDriver.h" 
namespace ohmmsqmc {

  /** Implements the VMC algorithm. 
   */
  class VMC: public QMCDriver {
  public:
    /// Constructor.
    VMC(MCWalkerConfiguration& w, 
	TrialWaveFunction& psi, 
	QMCHamiltonian& h, 
	xmlNodePtr q);

    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    VMC(const VMC& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMC& operator=(const VMC&) { return *this;}

    ///temporary storage for drift
    ParticleSet::ParticlePos_t drift;

    ///temporary storage for random displacement
    ParticleSet::ParticlePos_t deltaR;

    void advanceWalkerByWalker();

  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
