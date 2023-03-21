//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Refactored from: SimpleFixedNodeBranch.cpp
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCH_H
#define QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCH_H

#include <array>
#include <Configuration.h>
#include "OhmmsData/ParameterSet.h"
#include "Estimators/BlockHistogram.h"
#include "Estimators/accumulators.h"
#include "type_traits/template_types.hpp"
#include "Particle/Walker.h"
#include "QMCDrivers/Crowd.h"
#include "DMC/DMCRefEnergy.h"
#include <bitset>

namespace qmcplusplus
{
/** Manages the state of QMC sections and handles population control for DMCs
 *
 * \todo: Remove duplicate reading of Driver XML section with own copies of input
 *       parameters.
 * \todo: Rename, it is the only branching class so its name is too much
 * \todo: Use normal types for data members, don't be clever,
 *       the parameter enums violate KISS and make debugging annoying
 * \todo: Remove as much state as possible.
 *
 * QMCDriver object owns a SFNBranch to keep track of the
 * progress of a qmc section. It implements several methods to control the
 * population and trial energy during a DMC and evaluate the properties of
 * a population, e.g., energy, variance, population etc.
  * \see {http://qmcpack.cmscc.org/qmc-basics}
 *
 * Steps in SFNB states machine
 * 1. Construction (gets global walker number (rank or section wide?)
 * 3. put(reads driver XML node yet again)
 * 4. InitParam
 *   a. If TargetWalkers isn't known
 *      aa. allreduce and updates MCMW globalWalkers.
 *      bb. sets target walkers to whatever current total active walkers is.
 *   b. If not a restart
 *      aa. saves fixW and killWalker to internal params, otherwise just discards.
 *      bb. updates SFNB copy of MAX/MINWALKRS from walker controller, 
 *          these were set in constructer but I guess thats ony if this is a restart
 * 7. updateParam after each walker control branch (outside SFNBranch)
 *   a. Not first iter during warmup then call WalkerController branch.
 *      else call WC doNotBranch, returns pop_now
 *   b. copy a bunch a state from WC to SFNB (should be with respect to pop_now - number of released nodes)
 *   c. If using taueff update that based on acceptance ration and current tau.
 *   d. If not warmup calculate ETRIAL based on EREF and feedback * log(TargetWalkers) - log(pop_now) 
 *   e. set WC's TrialEnergy
 */
class SFNBranch : public QMCTraits
{
public:
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;

  /*! enum for booleans
   * \since 2008-05-05
   */
  enum
  {
    B_DMC = 0 /**< 1 for dmc, 0 for anything else */
    ,
    B_DMCSTAGE = 1 /**< 1 for main, 0 for wamrup  */
    ,
    B_POPCONTROL = 2 /**< 1 for the standard dmc, 0 for the comb method */
    ,
    B_USETAUEFF = 3 /**< 1 to use taueff accordning to JCP 93, 0 to use tau */
    ,
    B_CLEARHISTORY = 4 /**< 1 to clear the history */
    ,
    B_KILLNODES = 5 /**< 1 to kill walkers when a node crossing is detected */
    ,
    B_RESTART = 6 /**< 1 if restarting */
    ,
    B_RMC = 7 /**< 1 for rmc, 0 for anything else */
    ,
    B_RMCSTAGE = 8 /**< 1 for main, 0 for warmup */
    ,
    B_MODE_MAX = 10 /**< size of BranchMode */
  };

  using BranchModeType = std::bitset<B_MODE_MAX>;

  /*! enum for iParam std::bitset<B_IPARAM_MAX>
   * When introducing a new iParam, check if B_IPARAM_MAX is sufficiently large. Use multiples of 8
   * Why?  Much easier to use bool flags.  Are these ever serialized?
   */
  enum
  {
    B_WARMUPSTEPS = 0,      /**< warmup steps, valid when BranchMode[D_DMCSTAGE] == 0 */
    B_ENERGYUPDATEINTERVAL, /**< frequency of the trial energy updates, default 1 */
    B_COUNTER,              /**< counter for tracking object state */
    B_TARGETWALKERS,        /**< target total number of walkers per mpi group */
    B_BRANCHINTERVAL,       /**< interval between branch, see population control */
    B_IPARAM_MAX            /**< size of iParam */
  };

  /** input parameters of integer types
   * \since 2008-05-05
   */
  using IParamType = TinyVector<int, B_IPARAM_MAX>;

  /** enum for vParam 
   *
   *  Easy serialization is a relatively minor concern compared to the 
   *  annoyance this causes elsewhere.
   */
  enum class SimpleBranchVectorParameter
  {
    TAU = 0,
    TAUEFF, // effective time step
    ETRIAL, // Trial energy
    EREF,   // Center of the branching cutoff energy window
    ENOW,   // weighted average energy of the population in the current step
    BRANCHMAX,
    BRANCHCUTOFF,
    BRANCHFILTER,
    SIGMA2,
    SIGMA_BOUND,
    FEEDBACK,
    FILTERSCALE,
    VPARAM_MAX = 17 // four extra, why? Sloppy or undocumented hack?
  };
  using SBVP = SimpleBranchVectorParameter;

  /** controlling parameters of full precision real type
   *
   * Mostly internal
   */
  template<typename PAR_ENUM>
  struct VParams : public std::array<FullPrecRealType, static_cast<size_t>(PAR_ENUM::VPARAM_MAX)>
  {
    using Base = std::array<FullPrecRealType, static_cast<size_t>(PAR_ENUM::VPARAM_MAX)>;
    FullPrecRealType& operator[](PAR_ENUM sbvp) { return Base::operator[](static_cast<size_t>(sbvp)); }
    const FullPrecRealType& operator[](PAR_ENUM sbvp) const { return Base::operator[](static_cast<size_t>(sbvp)); }
  };
  using VParamType = VParams<SBVP>;

  ///Constructor
  SFNBranch(RealType tau, RealType feedback, DMCRefEnergyScheme);

  ///copy constructor
  SFNBranch(const SFNBranch& abranch) = delete;

  ~SFNBranch();

  inline bool phaseChanged(RealType psi0) const
  {
// TODO: remove ifdef
#if defined(QMC_COMPLEX)
    return false;
#else
    return std::cos(psi0) < std::numeric_limits<RealType>::epsilon();
#endif
  }

  /** increment QMCCounter
   *
   * QMCCounter is the number of times any QMC section is processed.
   */
  inline void advanceQMCCounter() { iParam[B_COUNTER]++; }
  inline void regressQMCCounter() { iParam[B_COUNTER]--; }

  /** initialize branching parameters
   * @param pop Population of Walkers
   * @param fixW true, if reconfiguration with the fixed number of walkers is used
   * @param killwalker 
   * @return number of copies to make in case targetwalkers changed
   */
  int initParam(const MCPopulation& population, FullPrecRealType ene, FullPrecRealType var, bool fixW, bool killwalker);

  /** return the bare branch weight
   *
   * This is equivalent to calling branchWeight(enew,eold,1.0,1.0)
   */
  inline RealType branchWeightBare(RealType enew, RealType eold) const
  {
    return std::exp(vParam[SBVP::TAUEFF] * (vParam[SBVP::ETRIAL] - 0.5 * (enew + eold)));
  }

  /** return the bare branch weight with a filtering using an energy window
   *
   * Cutoff values are set by the variance
   */
  inline RealType branchWeight(FullPrecRealType enew, FullPrecRealType eold) const
  {
    FullPrecRealType taueff_ = vParam[SBVP::TAUEFF] * 0.5;
    FullPrecRealType x       = std::max(vParam[SBVP::EREF] - enew, vParam[SBVP::EREF] - eold);
    if (x > vParam[SBVP::BRANCHMAX])
      taueff_ = 0.0;
    else if (x > vParam[SBVP::BRANCHCUTOFF])
      taueff_ *= (1.0 - (x - vParam[SBVP::BRANCHCUTOFF]) * vParam[SBVP::BRANCHFILTER]);
    return std::exp(taueff_ * (vParam[SBVP::ETRIAL] * 2.0 - enew - eold));
  }

  inline RealType symLinkAction(RealType logGf, RealType logGb, RealType enew, RealType eold) const
  {
    RealType driftaction = -0.5 * (logGf + logGb);
    //RealType energyaction =
    RealType taueff_ = vParam[SBVP::TAUEFF] * 0.5;
    RealType x       = std::max(vParam[SBVP::EREF] - enew, vParam[SBVP::EREF] - eold);
    if (x > vParam[SBVP::BRANCHMAX])
      taueff_ = 0.0;
    else if (x > vParam[SBVP::BRANCHCUTOFF])
      taueff_ *= (1.0 - (x - vParam[SBVP::BRANCHCUTOFF]) * vParam[SBVP::BRANCHFILTER]);
    RealType energyaction = taueff_ * (enew + eold);
    return driftaction + energyaction;
  }

  inline RealType symLinkActionBare(RealType logGf, RealType logGb, RealType enew, RealType eold) const
  {
    RealType driftaction  = -0.5 * (logGf + logGb);
    RealType taueff_      = vParam[SBVP::TAUEFF] * 0.5;
    RealType energyaction = taueff_ * (enew + eold);
    // RealType wavefunctionaction= -psinew + psiold;
    return driftaction + energyaction;
  }

  inline RealType DMCLinkAction(RealType enew, RealType eold) const
  {
    RealType taueff_ = vParam[SBVP::TAUEFF] * 0.5;
    RealType x       = std::max(vParam[SBVP::EREF] - enew, vParam[SBVP::EREF] - eold);
    if (x > vParam[SBVP::BRANCHMAX])
      taueff_ = 0.0;
    else if (x > vParam[SBVP::BRANCHCUTOFF])
      taueff_ *= (1.0 - (x - vParam[SBVP::BRANCHCUTOFF]) * vParam[SBVP::BRANCHFILTER]);
    return taueff_ * (enew + eold);
  }
  /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
   * @param enew new energy
   * @param eold old energy
   * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
   * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
   */
  inline RealType branchWeight(RealType enew, RealType eold, RealType scnew, RealType scold) const
  {
    FullPrecRealType s1 = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - enew) * scnew;
    FullPrecRealType s0 = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - eold) * scold;
    return std::exp(vParam[SBVP::TAUEFF] * 0.5 * (s1 + s0));
  }

  /** return the branch weight according to JCP1993 Umrigar et al. Appendix A
   * @param enew new energy
   * @param eold old energy
   * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
   * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
   * @param p acceptance ratio
   */
  inline RealType branchWeight(RealType enew, RealType eold, RealType scnew, RealType scold, RealType p) const
  {
    FullPrecRealType s1 = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - enew) * scnew;
    FullPrecRealType s0 = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - eold) * scold;
    return std::exp(vParam[SBVP::TAUEFF] * (p * 0.5 * (s1 - s0) + s0));
    //return std::exp(TauEff*(p*0.5*(sp-sq)+sq));
  }

  inline RealType getEref() const { return vParam[SBVP::EREF]; }
  inline RealType getEtrial() const { return vParam[SBVP::ETRIAL]; }
  inline RealType getTau() const { return vParam[SBVP::TAU]; }
  inline RealType getTauEff() const { return vParam[SBVP::TAUEFF]; }

  int getWarmupToDoSteps() const { return WarmUpToDoSteps; }
  /** perform branching
   * @param iter current step
   * @param w Walker configuration
   */
  void updateParamAfterPopControl(const MCDataType<FullPrecRealType>& wc_ensemble_prop, int Nelec);

  bool put(xmlNodePtr cur);

  /** create map between the parameter name and variables */
  void registerParameters();

  ///finalize the simulation
  void printStatus() const;

  friend std::ostream& operator<<(std::ostream& os, SFNBranch::VParamType& rhs);

private:
  ///set branch cutoff, max, filter
  void setBranchCutoff(FullPrecRealType variance,
                       FullPrecRealType targetSigma,
                       FullPrecRealType maxSigma,
                       int Nelec = 0);

  BranchModeType BranchMode;

  IParamType iParam;

  VParamType vParam;
  /// number of remaning steps in warmup, [0, iParam[B_WARMUPSTEPS]]
  int WarmUpToDoSteps;
  /// number of remaning steps in before adjusting ETRIAL, [0, iParam[B_ENERGYUPDATEINTERVAL]]
  int EtrialUpdateToDoSteps;
  ///save xml element
  xmlNodePtr myNode;
  /// collect energy and variance history
  DMCRefEnergy ref_energy_collector;
  ///a simple accumulator for energy
  accumulator_set<RealType> R2Accepted;
  ///a simple accumulator for energy
  accumulator_set<RealType> R2Proposed;
  ///a simple accumulator for reptation's center slice
  accumulator_set<RealType> R2Center;
  /////histogram of populations
  //BlockHistogram<RealType> DMCEnergyHist;
  ///scheme of branching cutoff
  std::string branching_cutoff_scheme;
  ///set of parameters
  ParameterSet m_param;
  ///string parameters
  std::vector<std::string> sParam;
};

std::ostream& operator<<(std::ostream& os, SFNBranch::VParamType& rhs);

} // namespace qmcplusplus
#endif
