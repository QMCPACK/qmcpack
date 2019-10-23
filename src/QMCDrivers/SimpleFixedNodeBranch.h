//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file SimpleFixedNodeBranch.h
 * @brief declare a handler of DMC branching
 *
 */
#ifndef QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H
#define QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H

#include <array>
#include <Configuration.h>
#include <OhmmsData/ParameterSet.h>
#include <Particle/MCWalkerConfiguration.h>
#include <Estimators/BlockHistogram.h>
#include <Estimators/accumulators.h>
#include "type_traits/template_types.hpp"
#include "Particle/Walker.h"
#include "QMCDrivers/WalkerControlBase.h"
#include "QMCDrivers/Crowd.h"
#include <Utilities/NewTimer.h>
#include <bitset>

namespace qmcplusplus
{
class EstimatorManagerBase;

/** Manages the state of QMC sections and handles population control for DMCs
 *
 * TODO: Remove duplicate reading of Driver XML section with own copies of input
 *       parameters.
 * TODO: Rename, it is the only branching class so its name is too much
 * TODO: Use normal types for data members, don't be clever to be clever
 *       the parameter enums violate KISS and make debugging annoying
 *
 * QMCDriver object owns a SimpleFixedNodeBranch to keep track of the
 * progress of a qmc section. It implements several methods to control the
 * population and trial energy during a DMC and evaluate the properties of
 * a population, e.g., energy, variance, population etc.
 * It owns WalkerController (pointer to a WalkerControlBase object) which
 * manages the population (killing and duplicating walkers) and
 * load balancing among multiple MPI tasks.
 * \see {http://qmcpack.cmscc.org/qmc-basics}
 */
struct SimpleFixedNodeBranch : public QMCTraits
{
  typedef SimpleFixedNodeBranch ThisType;
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

  /** booleans to set the branch modes
   * \since 2008-05-05
   */
  typedef std::bitset<B_MODE_MAX> BranchModeType;
  BranchModeType BranchMode;

  /*! enum for iParam std::bitset<B_IPARAM_MAX>
   * \since 2008-05-05
   *
   * When introducing a new iParam, check if B_IPARAM_MAX is sufficiently large. Use multiples of 8
   * Why?  Much easier to use bool flags.  Are these ever serialized?
   */
  enum
  {
    B_WARMUPSTEPS = 0 /**< warmup steps, valid when BranchMode[D_DMCSTAGE] == 0 */
    ,
    B_ENERGYUPDATEINTERVAL = 1 /**< frequency of the trial energy updates, default 1 */
    ,
    B_COUNTER = 2 /**< counter for tracking object state */
    ,
    B_TARGETWALKERS = 3 /**< target total number of walkers per mpi group */
    ,
    B_MAXWALKERS = 4 /**< maximum number of walkers per node */
    ,
    B_MINWALKERS = 5 /**< minimum number of walkers per node */
    ,
    B_BRANCHINTERVAL = 6 /**< interval between branch, see population control */
    ,
    B_IPARAM_MAX = 8 /**< size of iParam */
  };

  /** input parameters of integer types
   * \since 2008-05-05
   */
  typedef TinyVector<int, B_IPARAM_MAX> IParamType;
  IParamType iParam;

  /** enum for vParam 
   *
   *  Easy serialization is a relatively minor concern compared to the 
   *  annoyance this causes elsewhere.
   */
  enum class SimpleBranchVectorParameter
  {
    TAU = 0,
    TAUEFF,
    ETRIAL,
    EREF,
    ENOW,
    BRANCHMAX,
    BRANCHCUTOFF,
    BRANCHFILTER,
    SIGMA2,
    ACC_ENERGY,
    ACC_SAMPLES,
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
  VParamType vParam;

  /** number of remaning steps for a specific tasks
   *
   * set differently for BranchMode[B_DMCSTAGE]
   */
  int ToDoSteps;
  ///Feed*log(N)
  FullPrecRealType logN;
  ///save xml element
  xmlNodePtr myNode;
  ///WalkerController
  std::unique_ptr<WalkerControlBase> WalkerController;
  ///Backup WalkerController for mixed DMC
  std::unique_ptr<WalkerControlBase> BackupWalkerController;
  
  ///TODO: Should not be raw pointer 
  EstimatorManagerBase* MyEstimator;
  ///a simple accumulator for energy
  accumulator_set<FullPrecRealType> EnergyHist;
  ///a simple accumulator for variance
  accumulator_set<FullPrecRealType> VarianceHist;
  ///a simple accumulator for energy
  accumulator_set<RealType> R2Accepted;
  ///a simple accumulator for energy
  accumulator_set<RealType> R2Proposed;
  ///a simple accumulator for reptation's center slice
  accumulator_set<RealType> R2Center;
  /////histogram of populations
  //BlockHistogram<RealType> PopHist;
  /////histogram of populations
  //BlockHistogram<RealType> DMCEnergyHist;
  ///root name
  std::string RootName;
  ///scheme of branching cutoff
  std::string branching_cutoff_scheme;
  ///set of parameters
  ParameterSet m_param;
  ///string parameters
  std::vector<std::string> sParam;

  /// Used for the average scaling
  FullPrecRealType ScaleSum;
  unsigned long ScaleNum;
  //@TODO move these to private
  ///LogJacob
  RealType LogJacobRef;
  ///LogNorm
  std::vector<RealType> LogNorm;

  ///Releasednode
  bool RN;

  ///Constructor
  SimpleFixedNodeBranch(RealType tau, int nideal);

  ///copy constructor
  SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch);

  ~SimpleFixedNodeBranch() {}

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

  /** get the EstimatorManager */
  EstimatorManagerBase* getEstimatorManager() { return MyEstimator; }

  /** set the EstimatorManager
   * @param est estimator created by the first QMCDriver
   * */
  void setEstimatorManager(EstimatorManagerBase* est) { MyEstimator = est; }

  /** initialize  the WalkerController
   * @param w Walkers
   * @param tau timestep
   * @param fixW true, if reconfiguration with the fixed number of walkers is used
   * @return number of copies to make in case targetwalkers changed
   */
  int initWalkerController(MCWalkerConfiguration& mcwc, bool fixW, bool killwalker);

  /** initialize  the WalkerController
   * @param w Walkers
   * @param tau timestep
   * @param fixW true, if reconfiguration with the fixed number of walkers is used
   * @return number of copies to make in case targetwalkers changed
   */
  int initWalkerController(MCPopulation& pop, bool fixW, bool killwalker);

  /** initialize reptile stats
   *
   *
   */
  void initReptile(MCWalkerConfiguration& w);

  /** determine trial and reference energies
   */
  void checkParameters(MCWalkerConfiguration& w);

  /** determine trial and reference energies
   */
  void checkParameters(const int global_walkers, RefVector<MCPWalker>& walkers);

  /** return the bare branch weight
   *
   * This is equivalent to calling branchWeight(enew,eold,1.0,1.0)
   */
  inline RealType branchWeightBare(RealType enew, RealType eold) const
  {
    return std::exp(vParam[SBVP::TAUEFF] * (vParam[SBVP::ETRIAL] - 0.5 * (enew + eold)));
  }

  inline RealType branchWeightReleasedNode(RealType enew, RealType eold, RealType eref) const
  {
    if (BranchMode[B_DMCSTAGE])
      return std::exp(vParam[SBVP::TAU] * (eref - 0.5 * (enew + eold)));
    else
      return 1.0;
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

  /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
   * @param enew new energy
   * @param eold old energy
   * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
   * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
   */
  inline RealType branchWeightTau(RealType enew, RealType eold, RealType scnew, RealType scold, RealType taueff)
  {
    ScaleSum += scnew + scold;
    ScaleNum += 2;
    FullPrecRealType scavg = (ScaleNum > 10000) ? ScaleSum / (RealType)ScaleNum : 1.0;
    FullPrecRealType s1    = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - enew) * scnew / scavg;
    FullPrecRealType s0    = (vParam[SBVP::ETRIAL] - vParam[SBVP::EREF]) + (vParam[SBVP::EREF] - eold) * scold / scavg;
    return std::exp(taueff * 0.5 * (s1 + s0));
  }

  inline RealType getEref() const { return vParam[SBVP::EREF]; }
  inline RealType getEtrial() const { return vParam[SBVP::ETRIAL]; }
  inline RealType getTau() const { return vParam[SBVP::TAU]; }
  inline RealType getTauEff() const { return vParam[SBVP::TAUEFF]; }

  /** perform branching
   * @param iter current step
   * @param w Walker configuration
   */
  void branch(int iter, MCWalkerConfiguration& w);

  /** perform branching
   * @param iter current step
   * @param w Walker configuration
   */
  void branch(int iter, UPtrVector<Crowd>& crowds, MCPopulation& population);

  /** update RMC counters and running averages.
   * @param iter the iteration
   * @param w the walker ensemble
   * @param clones of the branch engine for OpenMP threads
   */
  void collect(int iter, MCWalkerConfiguration& w);

  /** restart averaging
   * @param counter Counter to determine the cummulative average will be reset.
   */
  void flush(int counter);

  /** reset the internal parameters */
  void reset();

  /** reset the internal parameters
   * @return new target population over old target population
   */
  int resetRun(xmlNodePtr cur);

  bool put(xmlNodePtr cur);

  /** write the state
   * @param fname name of the configuration file
   * @param overwrite NOT USED
   */
  void write(const std::string& fname, bool overwrite = true);

  void read(const std::string& fname);

  /** create map between the parameter name and variables */
  void registerParameters();

  ///start a run
  void start(const std::string& froot, bool append = false);
  ///finalize the simulation
  void finalize(MCWalkerConfiguration& w);
  ///finalize the simulation
  void finalize(const int global_walkers, RefVector<MCPWalker>& walkers);

  void setRN(bool rn);

private:
  ///default constructor (disabled)
  SimpleFixedNodeBranch() {}

  ///set branch cutoff, max, filter
  void setBranchCutoff(FullPrecRealType variance,
                       FullPrecRealType targetSigma,
                       FullPrecRealType maxSigma,
                       int Nelec = 0);
};

std::ostream& operator<<(std::ostream& os, SimpleFixedNodeBranch::VParamType& rhs);

} // namespace qmcplusplus
#endif
