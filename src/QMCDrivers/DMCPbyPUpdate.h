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
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_UPDATE_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_UPDATE_H
#include "Particle/MCWalkerConfiguration.h" 
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
namespace qmcplusplus {

  /** @ingroup QMCDrivers ParticleByParticle 
   *@brief Implements the DMC algorithm using particle-by-particle move. 
   */
  class DMCPbyPUpdate: public QMCTraits {

  public:

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::iterator WalkerIter_t;
    typedef Walker_t::Buffer_t              Buffer_t;
    typedef SimpleFixedNodeBranch           BranchEngineType;

    /// Constructor.
    DMCPbyPUpdate(ParticleSet& w, TrialWaveFunction& psi, 
        QMCHamiltonian& h, RandomGenerator_t& rg);
    ///destructor
    ~DMCPbyPUpdate();

    void initialize(WalkerIter_t it, WalkerIter_t it_end);
    void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);

    void advanceRejectNodeCrossing(WalkerIter_t it, WalkerIter_t it_end);
    void advanceKillNodeCrossing(WalkerIter_t it, WalkerIter_t it_end);
    void benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip);

    void resetRun(BranchEngineType* brancher);
    void resetBlock();
    void resetEtrial(RealType et);

    ///counter for number of moves accepted
    IndexType nAccept;
    ///counter for number of moves /rejected
    IndexType nReject; 
    ///Total number of the steps when all the particle moves are rejected.
    IndexType nAllRejected;
    ///Total number of node crossings per block
    IndexType nNodeCrossing;

  private:
    ///number of particles
    IndexType NumPtcl;
    ///timestep
    RealType Tau;

    ///Time-step factor \f$ 1/(2\Tau)\f$
    RealType m_oneover2tau;

    ///Time-step factor \f$ \sqrt{\Tau}\f$
    RealType m_sqrttau;

    ///walker ensemble
    ParticleSet& W;

    ///trial function
    TrialWaveFunction& Psi;

    ///Hamiltonian
    QMCHamiltonian& H;

    RandomGenerator_t& RandomGen;

    ///branch engine
    BranchEngineType* branchEngine;

    ///temporary storage for random displacement
    ParticleSet::ParticlePos_t deltaR;
    ParticleSet::ParticleGradient_t G, dG;
    ParticleSet::ParticleLaplacian_t L, dL;

    /// Copy Constructor (disabled)
    DMCPbyPUpdate(const DMCPbyPUpdate& a): W(a.W),Psi(a.Psi),H(a.H), RandomGen(a.RandomGen){ }
    /// Copy operator (disabled).
    DMCPbyPUpdate& operator=(const DMCPbyPUpdate&) { return *this;}

  };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
