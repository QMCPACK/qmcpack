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
#ifndef OHMMS_QMC_DMC_PARTICLEBYPARTCLE_H
#define OHMMS_QMC_DMC_PARTICLEBYPARTCLE_H
#include "QMCDrivers/QMCDriver.h" 
namespace ohmmsqmc {

  /** @ingroup QMCDrivers ParticleByParticle 
   *@brief Implements the DMC algorithm using particle-by-particle move. 
   */
  class DMCParticleByParticle: public QMCDriver {
  public:
    /// Constructor.
    DMCParticleByParticle(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    bool run();
    bool put(xmlNodePtr cur);
 
    void setBranchInfo(const string& aname);
  private:

    ///Column index for Populaton
    IndexType PopIndex;
    ///Column index for E_T
    IndexType EtrialIndex;
    ///hdf5 file name for Branch conditions
    string BranchInfo;

    /// Copy Constructor (disabled)
    DMCParticleByParticle(const DMCParticleByParticle& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    DMCParticleByParticle& operator=(const DMCParticleByParticle&) { return *this;}
  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
