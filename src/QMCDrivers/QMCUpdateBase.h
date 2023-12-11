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


/** @file
 * @brief Declare QMCUpdateBase class
 */
#ifndef QMCPLUSPLUS_QMCUPDATE_BASE_H
#define QMCPLUSPLUS_QMCUPDATE_BASE_H
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "GreenFunctionModifiers/DriftModifierBase.h"
#include "SimpleFixedNodeBranch.h"
#include "DriverDebugChecks.h"
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
class QMCUpdateBase : public QMCTraits
{
public:
  using Walker_t         = MCWalkerConfiguration::Walker_t;
  using WalkerIter_t     = MCWalkerConfiguration::iterator;
  using BranchEngineType = SimpleFixedNodeBranch;
#ifdef MIXED_PRECISION
  using mPosType    = TinyVector<OHMMS_PRECISION_FULL, DIM>;
  using mTensorType = Tensor<OHMMS_PRECISION_FULL, DIM>;
#else
  using mPosType    = PosType;
  using mTensorType = TensorType;
#endif

  ///number of steps per measurement
  int nSubSteps;
  /// determine additional checks for debugging purpose
  DriverDebugChecks debug_checks_ = DriverDebugChecks::ALL_OFF;
  std::string debug_checks_str_;
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
  ///spin mass
  RealType spinMass;
  ///use Drift
  bool UseDrift;

  /// Constructor.
  QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomBase<FullPrecRealType>& rg);
  ///Alt Constructor.
  QMCUpdateBase(MCWalkerConfiguration& w,
                TrialWaveFunction& psi,
                TrialWaveFunction& guide,
                QMCHamiltonian& h,
                RandomBase<FullPrecRealType>& rg);
  ///destructor
  virtual ~QMCUpdateBase();

  inline RealType acceptRatio() const
  {
    return static_cast<RealType>(nAccept) / static_cast<RealType>(nAccept + nReject);
  }

  /** reset the QMCUpdateBase parameters
   * @param brancher engine which handles branching
   *
   * Update time-step variables to move walkers
   */
  void resetRun(BranchEngineType* brancher,
                EstimatorManagerBase* est,
                TraceManager* traces,
                const DriftModifierBase* driftmodifer);

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
    Tau           = t;
    m_tauovermass = t * MassInvS[0];
    m_oneover2tau = 0.5 / (m_tauovermass);
    m_sqrttau     = std::sqrt(m_tauovermass);
  }

  inline RealType getSpinMass() { return spinMass; }

  inline void setSpinMass(RealType m) { spinMass = m; }

  inline void getLogs(std::vector<RealType>& logs) { Psi.getLogs(logs); }

  inline void set_step(int step) { W.current_step = step; }


  ///** start a run */
  void startRun(int blocks, bool record);
  /** stop a run */
  void stopRun();
  void stopRun2();

  /** prepare to start a block
   * @param steps number of steps within the block
   */
  void startBlock(int steps);

  /** stop a block
   */
  void stopBlock(bool collectall = true);

  /** set the multiplicity of the walkers to branch */
  void setMultiplicity(WalkerIter_t it, WalkerIter_t it_end);

  inline void setMultiplicity(Walker_t& awalker) const
  {
    constexpr RealType onehalf(0.5);
    constexpr RealType cone(1);
    RealType M = awalker.Weight;
    if (awalker.Age > MaxAge)
      M = std::min(onehalf, M);
    else if (awalker.Age > 0)
      M = std::min(cone, M);
    awalker.Multiplicity = M + RandomGen();
  }

  /** initialize Walker buffers for PbyP update
   */
  virtual void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);

  /** initialize Walker for walker update
   */
  virtual void initWalkers(WalkerIter_t it, WalkerIter_t it_end);

  /**  process options
   */
  virtual bool put(xmlNodePtr cur);

  inline void accumulate(WalkerIter_t it, WalkerIter_t it_end) { Estimators->accumulate(W, it, it_end); }

  /** advance walkers executed at each step
   *
   * Derived classes implement how to move walkers and accept/reject
   * moves.
   */
  virtual void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool recompute);

  ///move a walker
  virtual void advanceWalker(Walker_t& thisWalker, bool recompute) = 0;

  virtual RealType advanceWalkerForEE(Walker_t& w1,
                                      std::vector<PosType>& dR,
                                      std::vector<int>& iats,
                                      std::vector<int>& rs,
                                      std::vector<RealType>& ratios)
  {
    return 0.0;
  };

  ///normalization offset for cs type runs.
  RealType csoffset;

  //       virtual void estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
  //     , std::vector<MCWalkerConfiguration*>& wclone
  //     , std::vector<QMCHamiltonian*>& hclone
  //     , std::vector<RandomGenerator*>& rng
  //     , std::vector<RealType>& ratio_i_0){};
  int RMC_checkIndex(int N, int NMax)
  {
    if (N < 0)
      return N + NMax;
    else if (N >= NMax)
      return N - NMax;
    else
      return N;
  }

  void RMC_checkWalkerBounds(WalkerIter_t& it, WalkerIter_t first, WalkerIter_t last)
  {
    if (it >= last)
      it -= (last - first);
    else if (it < first)
      it += (last - first);
  }

  inline RealType logBackwardGF(const ParticleSet::ParticlePos& displ)
  {
    RealType logGb = 0.0;
    for (int iat = 0; iat < W.getTotalNum(); ++iat)
    {
      RealType mass_over_tau = 1.0 / (SqrtTauOverMass[iat] * SqrtTauOverMass[iat]);
      logGb += 0.5 * dot(displ[iat], displ[iat]) * mass_over_tau;
    }
    return -logGb;
  }

public:
  ///traces
  TraceManager* Traces;

protected:
  ///update particle-by-particle
  bool UpdatePbyP;
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
  RandomBase<FullPrecRealType>& RandomGen;
  ///branch engine, stateless reference to the one in QMCDriver
  const BranchEngineType* branchEngine;
  ///drift modifer, stateless reference to the one in QMCDriver
  const DriftModifierBase* DriftModifier;
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
  ///temporary storage for drift
  ParticleSet::ParticlePos drift;
  ///temporary storage for random displacement
  ParticleSet::ParticlePos deltaR;
  ///temporart storage for spin displacement
  ParticleSet::ParticleScalar deltaS;
  ///storage for differential gradients for PbyP update
  ParticleSet::ParticleGradient G, dG;
  ///storage for differential laplacians for PbyP update
  ParticleSet::ParticleLaplacian L, dL;

  /** evaluate the ratio of scaled velocity and velocity
   * @param g gradient
   * @param gscaled scaled gradient
   * @return the ratio
   */
  RealType getNodeCorrection(const ParticleSet::ParticleGradient& g, ParticleSet::ParticlePos& gscaled);

  ///copy constructor (disabled)
  QMCUpdateBase(const QMCUpdateBase&) = delete;

  /// check logpsi and grad and lap against values computed from scratch
  static void checkLogAndGL(ParticleSet& pset, TrialWaveFunction& twf, const std::string_view location);

private:
  ///set default parameters
  void setDefaults();
  /// Copy operator (disabled).
  QMCUpdateBase& operator=(const QMCUpdateBase&) { return *this; }
  ///
  NewTimer& initWalkers_timer_;
};
} // namespace qmcplusplus

#endif
