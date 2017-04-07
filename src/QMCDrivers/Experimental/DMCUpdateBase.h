//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file DMCUpdateBase
 * @brief Declare DMCUpdateBase class
 */
#ifndef QMCPLUSPLUS_DMC_UPDATE_BASE_H
#define QMCPLUSPLUS_DMC_UPDATE_BASE_H
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"

namespace qmcplusplus
{

/** @ingroup DMC
 * @brief Base class for update methods for each step
 *
 * DMCUpdateBase provides the common functions to update all the walkers for each time step.
 * Derived classes should implement advanceWalkers to complete a step.
 */
class DMCUpdateBase: public QMCTraits
{

public:

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::iterator WalkerIter_t;
  typedef SimpleFixedNodeBranch           BranchEngineType;

  ///MaxAge>0 indicates branch is done
  IndexType MaxAge;
  ///counter for number of moves accepted
  IndexType nAccept;
  ///counter for number of moves /rejected
  IndexType nReject;
  ///Total number of the steps when all the particle moves are rejected.
  IndexType nAllRejected;
  ///Total number of node crossings per block
  IndexType nNodeCrossing;
  ///Total numer of non-local moves accepted
  IndexType NonLocalMoveAccepted;

  /// Constructor.
  DMCUpdateBase(ParticleSet& w, TrialWaveFunction& psi,
                QMCHamiltonian& h, RandomGenerator_t& rg);
  ///destructor
  virtual ~DMCUpdateBase();

  inline RealType acceptRatio() const
  {
    return static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject);
  }
  /** reset the DMCUpdateBase parameters
   * @param brancher engine which handles branching
   *
   * Update time-step variables to move walkers
   */
  void resetRun(BranchEngineType* brancher);

  /** reset the trial energy */
  void resetEtrial(RealType et);

  /** prepare to start a block
   */
  void startBlock();

  /** set the multiplicity of the walkers to branch */
  void setMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

  /** initialize Walker buffers for PbyP update
   */
  void initWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** update Walker buffers for PbyP update
   */
  void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** simple routine to test the performance
   */
  void benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip);

  /** advance walkers executed at each step
   *
   * Derived classes implement how to move walkers and accept/reject
   * moves.
   */
  virtual void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end)=0;


protected:
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

  ///random number generator
  RandomGenerator_t& RandomGen;

  ///branch engine
  BranchEngineType* branchEngine;

  ///temporary storage for drift
  ParticleSet::ParticlePos_t drift;

  ///temporary storage for random displacement
  ParticleSet::ParticlePos_t deltaR;

  ///storage for differential gradients and laplacians for PbyP update
  ParticleSet::ParticleGradient_t G, dG;
  ParticleSet::ParticleLaplacian_t L, dL;

  /// Copy Constructor (disabled)
  DMCUpdateBase(const DMCUpdateBase& a): W(a.W),Psi(a.Psi),H(a.H), RandomGen(a.RandomGen) { }
  /// Copy operator (disabled).
  DMCUpdateBase& operator=(const DMCUpdateBase&)
  {
    return *this;
  }
};
}

#endif
