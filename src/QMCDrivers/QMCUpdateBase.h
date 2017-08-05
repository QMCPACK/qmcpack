//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
#include "Estimators/EstimatorManagerBase.h"
namespace qmcplusplus
{

class TraceManager;
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
#ifdef MIXED_PRECISION
  typedef TinyVector<OHMMS_PRECISION_FULL, DIM>         mPosType;
#else
  typedef PosType                         mPosType;
#endif

  ///If true, terminate the simulation
  bool BadState;
  ///number of steps per measurement
  int nSubSteps;
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
  ///timestep
  RealType Tau;
  ///use Drift
  bool UseDrift;

  /// Constructor.
  QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                RandomGenerator_t& rg);
  ///Alt Constructor.
  QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h,
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
  void resetRun(BranchEngineType* brancher, EstimatorManagerBase* est);

  void resetRun(BranchEngineType* brancher, EstimatorManagerBase* est, TraceManager* traces);

  inline RealType getTau()
  {
    //SpeciesSet tspecies(W.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //return m_tauovermass*mass;
    return Tau;
  }

  inline void setTau(RealType t)
  {
    //SpeciesSet tspecies(W.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //RealType oneovermass = 1.0/mass;
    //RealType oneoversqrtmass = std::sqrt(oneovermass);
//  //     Tau=brancher->getTau();
//  //     assert (Tau==i);
    //m_tauovermass = i/mass;
    Tau=t;
    m_tauovermass = t*MassInvS[0];
    m_oneover2tau = 0.5/(m_tauovermass);
    m_sqrttau = std::sqrt(m_tauovermass);
  }

  inline void getLogs(std::vector<RealType>& logs)
  {
    Psi.getLogs(logs);
  }

  inline void set_step(int step)
  {
    W.current_step = step;
  }



  ///** start a run */
  void startRun(int blocks, bool record);
  /** stop a run */
  void stopRun();
  void stopRun2();
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

  inline void setMultiplicity(Walker_t& awalker) const
  {
    CONSTEXPR RealType onehalf(0.5);
    CONSTEXPR RealType cone(1);
    RealType M=awalker.Weight;
    if (awalker.Age>MaxAge)
      M = std::min(onehalf,M);
    else
      if (awalker.Age > 0) M = std::min(cone,M);
    awalker.Multiplicity = M + RandomGen();
  }

  /** set the multiplicity of the walkers to branch */
  void setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

  /** initialize Walker buffers for PbyP update
   */
  virtual void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

  /** initalize Walker for walker update
   */
  virtual void initWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** update Walker buffers for PbyP update
   */
  void updateWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /** simple routine to test the performance
   */
  void benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip);

  /**  process options
   */
  virtual bool put(xmlNodePtr cur);

  inline void accumulate(WalkerIter_t it, WalkerIter_t it_end)
  {
    Estimators->accumulate(W,it,it_end);
  }

  /** advance walkers executed at each step
   *
   * Derived classes implement how to move walkers and accept/reject
   * moves.
   */
  virtual void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool recompute);

  ///move a walker
  virtual void advanceWalker(Walker_t& thisWalker, bool recompute)=0;

  virtual RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios)
  {
    return 0.0;
  };
//       virtual RealType advanceWalkerForCSEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios, std::vector<RealType>& weights, std::vector<RealType>& logs ) {return 0.0;};
  virtual void setLogEpsilon(RealType eps) {};
//       virtual void advanceCSWalkers(std::vector<TrialWaveFunction*>& pclone, std::vector<MCWalkerConfiguration*>& wclone, std::vector<QMCHamiltonian*>& hclone, std::vector<RandomGenerator_t*>& rng, std::vector<RealType>& c_i){};

  ///normalization offset for cs type runs.
  RealType csoffset;

//       virtual void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
//     , std::vector<MCWalkerConfiguration*>& wclone
//     , std::vector<QMCHamiltonian*>& hclone
//     , std::vector<RandomGenerator_t*>& rng
//     , std::vector<RealType>& ratio_i_0){};
  int RMC_checkIndex(int N, int NMax)
  {
    if(N<0)
      return N+NMax;
    else
      if (N>=NMax)
        return N-NMax;
      else
        return N;
  }

  void RMC_checkWalkerBounds(WalkerIter_t& it, WalkerIter_t first, WalkerIter_t last)
  {
    if (it>=last)
      it-=(last-first);
    else
      if (it<first)
        it+=(last-first);
  }

  inline RealType logBackwardGF(const ParticleSet::ParticlePos_t& displ)
  {
    RealType logGb=0.0;
    for(int iat=0; iat<W.getTotalNum(); ++iat)
    {
      RealType mass_over_tau = 1.0/(SqrtTauOverMass[iat]*SqrtTauOverMass[iat]);
      logGb += 0.5*dot(displ[iat],displ[iat])*mass_over_tau;
    }
    return -logGb;
  }

public:
  ///traces
  TraceManager* Traces;
protected:
  ///update particle-by-particle
  bool UpdatePbyP;
  ///use T-moves
  int UseTMove;
  ///number of particles
  IndexType NumPtcl;
  ///Time-step factor \f$ 1/(2\tau)\f$
  RealType m_oneover2tau;
  ///Time-step factor \f$ \sqrt{\tau}\f$
  RealType m_sqrttau;
  ///tau/mass
  RealType m_tauovermass;
  ///maximum displacement^2
  RealType m_r2max;
  ///walker ensemble
  MCWalkerConfiguration& W;
  ///trial function
  TrialWaveFunction& Psi;
  ///guide function
  TrialWaveFunction& Guide;
  ///Hamiltonian
  QMCHamiltonian& H;
  ///random number generator
  RandomGenerator_t& RandomGen;
  ///branch engine
  BranchEngineType* branchEngine;
  ///estimator
  EstimatorManagerBase* Estimators;
  ///parameters
  ParameterSet myParams;
  ///1/Mass per species
  std::vector<RealType> MassInvS;
  ///1/Mass per particle
  std::vector<RealType> MassInvP;
  ///sqrt(tau/Mass) per particle
  std::vector<RealType> SqrtTauOverMass;
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

  ///copy constructor
  QMCUpdateBase(const QMCUpdateBase& a);

  /** a VMC step to randomize awalker
   */
  void randomize(Walker_t& awalker);

private:

  ///set default parameters
  void setDefaults();
  /// Copy operator (disabled).
  QMCUpdateBase& operator=(const QMCUpdateBase&)
  {
    return *this;
  }
  ///
  NewTimer* InitWalkersTimer;
};
}

#endif
