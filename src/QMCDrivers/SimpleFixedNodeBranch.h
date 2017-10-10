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

#include <Configuration.h>
#include <OhmmsData/ParameterSet.h>
#include <Particle/MCWalkerConfiguration.h>
#include <Estimators/BlockHistogram.h>
#include <Estimators/accumulators.h>
#include <Utilities/NewTimer.h>
#include <bitset>

#ifdef HAVE_ADIOS
#include <adios.h>
#endif

namespace qmcplusplus
{

class WalkerControlBase;
class EstimatorManagerBase;

/** Manages the state of QMC sections and handles population control for DMCs
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
struct SimpleFixedNodeBranch: public QMCTraits
{

  typedef SimpleFixedNodeBranch ThisType;

  /*! enum for booleans
   * \since 2008-05-05
   */
  enum
  {
    B_DMC=0              /**< 1 for dmc, 0 for anything else */
      , B_DMCSTAGE=1     /**< 1 for main, 0 for wamrup  */
      , B_POPCONTROL=2   /**< 1 for the standard dmc, 0 for the comb method */
      , B_USETAUEFF=3    /**< 1 to use taueff accordning to JCP 93, 0 to use tau */
      , B_CLEARHISTORY=4 /**< 1 to clear the history */
      , B_KILLNODES=5    /**< 1 to kill walkers when a node crossing is detected */
      , B_RESTART=6      /**< 1 if restarting */
      , B_RMC=7          /**< 1 for rmc, 0 for anything else */
      , B_RMCSTAGE=8    /**< 1 for main, 0 for warmup */
      , B_MODE_MAX=10     /**< size of BranchMode */
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
   */
  enum
  {
    B_WARMUPSTEPS=0              /**< warmup steps, valid when BranchMode[D_DMCSTAGE] == 0 */
      , B_ENERGYUPDATEINTERVAL=1 /**< frequency of the trial energy updates, default 1 */
      , B_COUNTER=2              /**< counter for tracking object state */
      , B_TARGETWALKERS=3        /**< target total number of walkers per mpi group */
      , B_MAXWALKERS=4           /**< maximum number of walkers per node */
      , B_MINWALKERS=5           /**< minimum number of walkers per node */
      , B_BRANCHINTERVAL=6       /**< interval between branch, see population control */
      , B_IPARAM_MAX=8           /**< size of iParam */
  };

  /** input parameters of integer types
   * \since 2008-05-05
   */
  typedef TinyVector<int,B_IPARAM_MAX> IParamType;
  IParamType iParam;

  /*! enum for vParam */
  enum
  {
    B_TAU=0, B_TAUEFF , B_ETRIAL , B_EREF
    , B_ENOW, B_BRANCHMAX, B_BRANCHCUTOFF, B_BRANCHFILTER
    , B_SIGMA2, B_ACC_ENERGY, B_ACC_SAMPLES, B_FEEDBACK
    , B_FILTERSCALE, B_VPARAM_MAX=17
  };

  /** controlling parameters of real type
   *
   * Mostly internal
   */
  typedef TinyVector<EstimatorRealType,B_VPARAM_MAX> VParamType;
  VParamType vParam;

  /** number of remaning steps for a specific tasks
   *
   * set differently for BranchMode[B_DMCSTAGE]
   */
  int ToDoSteps;
  ///Feed*log(N)
  EstimatorRealType logN;
  ///save xml element
  xmlNodePtr myNode;
  ///WalkerController
  WalkerControlBase* WalkerController;
  ///Backup WalkerController for mixed DMC
  WalkerControlBase* BackupWalkerController;
  ///EstimatorManager
  EstimatorManagerBase*  MyEstimator;
  ///a simple accumulator for energy
  accumulator_set<EstimatorRealType> EnergyHist;
  ///a simple accumulator for variance
  accumulator_set<EstimatorRealType> VarianceHist;
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
  EstimatorRealType ScaleSum;
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
  inline void advanceQMCCounter()
  {
    iParam[B_COUNTER]++;
  }
  inline void regressQMCCounter()
  {
    iParam[B_COUNTER]--;
  }

  /** get the EstimatorManager */
  EstimatorManagerBase* getEstimatorManager()
  {
    return MyEstimator;
  }

  /** set the EstimatorManager
   * @param est estimator created by the first QMCDriver
   * */
  void setEstimatorManager(EstimatorManagerBase* est)
  {
    MyEstimator=est;
  }

  /** initialize  the WalkerController
   * @param w Walkers
   * @param tau timestep
   * @param fixW true, if reconfiguration with the fixed number of walkers is used
   * @return number of copies to make in case targetwalkers changed
   */
  int initWalkerController(MCWalkerConfiguration& w, bool fixW, bool killwalker);
  //void initWalkerController(MCWalkerConfiguration& w, RealType tau, bool fixW=false, bool killwalker=false);

  /** initialize reptile stats
   *
   *
   */
  void initReptile(MCWalkerConfiguration& w);

  /** determine trial and reference energies
   */
  void checkParameters(MCWalkerConfiguration& w);

  /** return the bare branch weight
   *
   * This is equivalent to calling branchWeight(enew,eold,1.0,1.0)
   */
  inline RealType branchWeightBare(RealType enew, RealType eold) const
  {
    return std::exp(vParam[B_TAUEFF]*(vParam[B_ETRIAL]-0.5*(enew+eold)));
  }

  inline RealType branchWeightReleasedNode(RealType enew, RealType eold, RealType eref) const
  {
    if(BranchMode[B_DMCSTAGE])
      return std::exp(vParam[B_TAU]*(eref-0.5*(enew+eold)));
    else
      return 1.0;
  }

  /** return the bare branch weight with a filtering using an energy window
   *
   * Cutoff values are set by the variance
   */
  inline RealType branchWeight(EstimatorRealType enew, EstimatorRealType eold) const
  {
    EstimatorRealType taueff_=vParam[B_TAUEFF]*0.5;
    EstimatorRealType x=std::max(vParam[B_EREF]-enew,vParam[B_EREF]-eold);
    if(x>vParam[B_BRANCHMAX])
      taueff_=0.0;
    else if(x>vParam[B_BRANCHCUTOFF])
      taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
    return std::exp(taueff_*(vParam[B_ETRIAL]*2.0-enew-eold));
  }

  inline RealType symLinkAction(RealType logGf, RealType logGb, RealType enew, RealType eold) const
  {
    RealType driftaction = -0.5*(logGf+logGb);
    //RealType energyaction =
    RealType taueff_=vParam[B_TAUEFF]*0.5;
    RealType x=std::max(vParam[B_EREF]-enew,vParam[B_EREF]-eold);
    if(x>vParam[B_BRANCHMAX])
      taueff_=0.0;
    else if(x>vParam[B_BRANCHCUTOFF])
      taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
    RealType energyaction = taueff_*(enew+eold);
    return driftaction+energyaction;
  }

  inline RealType symLinkActionBare(RealType logGf, RealType logGb, RealType enew, RealType eold) const
  {
    RealType driftaction = -0.5*(logGf+logGb);
    RealType taueff_=vParam[B_TAUEFF]*0.5;
    RealType energyaction = taueff_*(enew+eold);
    // RealType wavefunctionaction= -psinew + psiold;
    return driftaction+energyaction;
  }

  inline RealType DMCLinkAction(RealType enew, RealType eold) const
  {
    RealType taueff_=vParam[B_TAUEFF]*0.5;
    RealType x=std::max(vParam[B_EREF]-enew,vParam[B_EREF]-eold);
    if(x>vParam[B_BRANCHMAX])
      taueff_=0.0;
    else if(x>vParam[B_BRANCHCUTOFF])
      taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
    return taueff_*(enew+eold);
  }
  /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
   * @param enew new energy
   * @param eold old energy
   * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
   * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
   */
  inline RealType branchWeight(RealType enew, RealType eold, RealType scnew, RealType scold) const
  {
    EstimatorRealType s1=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-enew)*scnew;
    EstimatorRealType s0=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-eold)*scold;
    return std::exp(vParam[B_TAUEFF]*0.5*(s1+s0));
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
    EstimatorRealType s1=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-enew)*scnew;
    EstimatorRealType s0=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-eold)*scold;
    return std::exp(vParam[B_TAUEFF]*(p*0.5*(s1-s0)+s0));
    //return std::exp(TauEff*(p*0.5*(sp-sq)+sq));
  }

  /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
   * @param enew new energy
   * @param eold old energy
   * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
   * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
   */
  inline RealType branchWeightTau(RealType enew, RealType eold, RealType scnew, RealType scold,
                                  RealType taueff)
  {
    ScaleSum += scnew + scold;
    ScaleNum +=2;
    EstimatorRealType scavg = (ScaleNum > 10000) ? ScaleSum/(RealType)ScaleNum : 1.0;
    EstimatorRealType s1=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-enew)*scnew/scavg;
    EstimatorRealType s0=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-eold)*scold/scavg;
    return std::exp(taueff*0.5*(s1+s0));
  }

  inline RealType getEref() const
  {
    return vParam[B_EREF];
  }
  inline RealType getEtrial() const
  {
    return vParam[B_ETRIAL];
  }
  inline RealType getTau() const
  {
    return vParam[B_TAU];
  }
  inline RealType getTauEff() const
  {
    return vParam[B_TAUEFF];
  }

  inline void setTrialEnergy(RealType etot, RealType wtot=1.0)
  {
    vParam[B_EREF]=vParam[B_ETRIAL]=etot/wtot;
    //Eref=Etrial=etot/wtot;
  }

  /** perform branching
   * @param iter current step
   * @param w Walker configuration
   */
  void branch(int iter, MCWalkerConfiguration& w);

  /** perform branching
   * @param iter the iteration
   * @param w the walker ensemble
   * @param clones of the branch engine for OpenMP threads
   */
  void branch(int iter, MCWalkerConfiguration& w, std::vector<ThisType*>& clones);

  /** update RMC counters and running averages.
   * @param iter the iteration
   * @param w the walker ensemble
   * @param clones of the branch engine for OpenMP threads
   */

  void collect(int iter, MCWalkerConfiguration& w);
  void collect(int iter, MCWalkerConfiguration& w, std::vector<ThisType*>& clones);

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
  void write(const std::string& fname, bool overwrite=true);

#ifdef HAVE_ADIOS
  void save_energy();
#endif

  void read(const std::string& fname);

  /** create map between the parameter name and variables */
  void registerParameters();

  ///start a run
  void start(const std::string& froot, bool append=false);
  ///finalize the simulation
  void finalize(MCWalkerConfiguration& w);

  void setRN(bool rn);


//     void storeConfigsForForwardWalking(MCWalkerConfiguration& w);
//     void clearConfigsForForwardWalking( );
//     void debugFWconfig();
//     WalkerControlBase* getWalkerController()
//     {
//       return WalkerController;
//     }

private:
  ///default constructor (disabled)
  SimpleFixedNodeBranch() {}

  ///disable use by external users
  //void write(hid_t grp, bool append=false);
  //void read(hid_t grp);

  ///set branch cutoff, max, filter
  void setBranchCutoff(EstimatorRealType variance, EstimatorRealType targetSigma, EstimatorRealType maxSigma, int Nelec=0);
};

}
#endif

