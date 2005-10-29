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
#include "Utilities/OhmmsInfo.h"
namespace ohmmsqmc {

  /** @ingroup QMCDrivers ParticleByParticle 
   *@brief Implements the DMC algorithm using particle-by-particle move. 
   */
  class DMCParticleByParticle: public QMCDriver {
  public:
    /// Constructor.
    DMCParticleByParticle(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
			  QMCHamiltonian& h);
    ///destructor
    ~DMCParticleByParticle();

    bool run();
    bool put(xmlNodePtr cur);
 
    //void setBranchInfo(const string& aname);
  private:

    ///Index to determine what to do when node crossing is detected
    IndexType KillNodeCrossing;
    ///Column index for Populaton
    IndexType PopIndex;
    ///Column index for E_T
    IndexType EtrialIndex;
    ///Total number of accepted steps per block
    IndexType nAcceptTot;
    ///Total number of rejected steps per block
    IndexType nRejectTot;
    ///Total number of the steps when all the particle moves are rejected.
    IndexType nAllRejected;
    ///Total number of node crossings per block
    IndexType nNodeCrossing;
    ///hdf5 file name for Branch conditions
    string BranchInfo;
    ///input string to determine kill walkers or not
    string KillWalker;
    ///input string to determine swap walkers among mpi processors
    string SwapWalkers;
    ParticleSet::ParticleGradient_t G, dG;
    ParticleSet::ParticleLaplacian_t L, dL;

    /// Copy Constructor (disabled)
    DMCParticleByParticle(const DMCParticleByParticle& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    DMCParticleByParticle& operator=(const DMCParticleByParticle&) { return *this;}

    void advanceRejectNodeCrossing(int nat);
    void advanceKillNodeCrossing(int nat);

  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
