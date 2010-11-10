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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file QMCUpdateBase
 * @brief Declare QMCUpdateBase class
 */
#ifndef QMCPLUSPLUS_QMCUPDATE_BASE_H
#define QMCPLUSPLUS_QMCUPDATE_BASE_H
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
//#define ENABLE_COMPOSITE_ESTIMATOR
//#include "Estimators/CompositeEstimators.h"
#include "Estimators/EstimatorManager.h"

namespace qmcplusplus
  {

  /** @ingroup QMC
   * @brief Base class for update methods for each step
   *
   * QMCUpdateBase provides the common functions to update all the walkers for each time step.
   * Derived classes should implement advanceWalkers to complete a step.
   */
  class QMCUpdateBase: public QMCTraits
    {

    public:

      typedef MCWalkerConfiguration::Walker_t Walker_t;
      typedef MCWalkerConfiguration::iterator WalkerIter_t;
      typedef SimpleFixedNodeBranch           BranchEngineType;

      ///MaxAge>0 indicates branch is done
      IndexType MaxAge;
      ///counter for number of moves accepted
      IndexType nAccept;
      ///counter for number of moves rejected
      IndexType nReject;
      ///Total number of the steps when all the particle moves are rejected.
      IndexType nAllRejected;
      ///Total number of node crossings per block
      IndexType nNodeCrossing;
      ///Total numer of non-local moves accepted
      IndexType NonLocalMoveAccepted;
      /// Constructor.
      QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                    RandomGenerator_t& rg);
      ///destructor
      virtual ~QMCUpdateBase();

      inline RealType acceptRatio() const
        {
          return static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject);
        }

      /** reset the QMCUpdateBase parameters
       * @param brancher engine which handles branching
       *
       * Update time-step variables to move walkers
       */
      void resetRun(BranchEngineType* brancher, EstimatorManager* est);

      inline void setTau(RealType i)
      {
        SpeciesSet tspecies(W.getSpeciesSet());
        int massind=tspecies.addAttribute("mass");
        RealType mass = tspecies(massind,0);

        RealType oneovermass = 1.0/mass;
        RealType oneoversqrtmass = std::sqrt(oneovermass);
//       Tau=brancher->getTau();
//       assert (Tau==i);
        m_tauovermass = i/mass;
        m_oneover2tau = 0.5/(m_tauovermass);
        m_sqrttau = std::sqrt(m_tauovermass);
      }

      ///** start a run */
      void startRun(int blocks, bool record);
      /** stop a run */
      void stopRun();
      /** reset the trial energy */
      void resetEtrial(RealType et);

      /** prepare to start a block
       * @param steps number of steps within the block
       */
      void startBlock(int steps);

      /** stop a block
       */
      void stopBlock(bool collectall=true);

      /** set the multiplicity of the walkers to branch */
      void setMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

      /** set the multiplicity of the walkers to branch */
      void setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

      /** initialize Walker buffers for PbyP update
       */
      virtual void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

      /** initalize Walker for walker update
       */
      void initWalkers(WalkerIter_t it, WalkerIter_t it_end);

      /** update Walker buffers for PbyP update
       */
      void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);

      /** simple routine to test the performance
       */
      void benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip);

      /**  process options
       */
      bool put(xmlNodePtr cur);

      inline void accumulate(WalkerIter_t it, WalkerIter_t it_end)
      {
        Estimators->accumulate(W,it,it_end);
      }

      /** advance walkers executed at each step
       *
       * Derived classes implement how to move walkers and accept/reject
       * moves.
       */
      virtual void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)=0;
      
      virtual void setLogEpsilon(RealType eps) {};
      virtual void advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone, vector<RandomGenerator_t*>& rng){};

    protected:
      ///update particle-by-particle
      bool UpdatePbyP;
      ///use T-moves
      bool UseTMove;
      ///number of particles
      IndexType NumPtcl;
      ///number of steps per measurement
      int nSubSteps;
      ///timestep
      RealType Tau;
      ///Time-step factor \f$ 1/(2\Tau)\f$
      RealType m_oneover2tau;
      ///Time-step factor \f$ \sqrt{\Tau}\f$
      RealType m_sqrttau;
      ///tau/mass
      RealType m_tauovermass;
      ///maximum displacement^2
      RealType m_r2max;
      ///walker ensemble
      MCWalkerConfiguration& W;
      ///trial function
      TrialWaveFunction& Psi;
      ///Hamiltonian
      QMCHamiltonian& H;
      ///random number generator
      RandomGenerator_t& RandomGen;
      ///branch engine
      BranchEngineType* branchEngine;
      ///estimator
      EstimatorManager* Estimators;
      ///parameters
      ParameterSet myParams;
      ///non local operator
      NonLocalTOperator nonLocalOps;
      ///temporary storage for drift
      ParticleSet::ParticlePos_t drift;
      ///temporary storage for random displacement
      ParticleSet::ParticlePos_t deltaR;
      ///storage for differential gradients for PbyP update
      ParticleSet::ParticleGradient_t G, dG;
      ///storage for differential laplacians for PbyP update
      ParticleSet::ParticleLaplacian_t L, dL;

      /** evaluate the ratio of scaled velocity and velocity
       * @param g gradient
       * @param gscaled scaled gradient
       * @return the ratio
       */
      RealType getNodeCorrection(const ParticleSet::ParticleGradient_t& g, ParticleSet::ParticlePos_t& gscaled);

      /// Copy Constructor (disabled)
      QMCUpdateBase(const QMCUpdateBase& a): W(a.W),Psi(a.Psi),H(a.H), RandomGen(a.RandomGen) { }
      /// Copy operator (disabled).
      QMCUpdateBase& operator=(const QMCUpdateBase&)
      {
        return *this;
      }
    };
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: QMCUpdateBase.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
