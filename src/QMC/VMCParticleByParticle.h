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
#ifndef OHMMS_QMC_VMC_PARTICLEBYPARTCLE_H
#define OHMMS_QMC_VMC_PARTICLEBYPARTCLE_H
#include "QMC/QMCDriver.h" 
namespace ohmmsqmc {

  /** Implements the VMC algorithm using particle-by-particle move. */
  class VMCParticleByParticle: public QMCDriver {
  public:
    /// Constructor.
    VMCParticleByParticle(MCWalkerConfiguration& w, 
			  TrialWaveFunction& psi, 
			  QMCHamiltonian& h, 
			  xmlNodePtr q);
    bool run();
    bool put(xmlNodePtr cur);
 
  private:
    /// Copy Constructor (disabled)
    VMCParticleByParticle(const VMCParticleByParticle& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    VMCParticleByParticle& operator=(const VMCParticleByParticle&) { return *this;}
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
