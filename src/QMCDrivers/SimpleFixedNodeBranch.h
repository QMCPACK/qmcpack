//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
/**@file SimpleFixedNodeBranch.h
 * @brief declare a handler of DMC branching
 *
 */
#ifndef QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H
#define QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Estimators/BlockHistogram.h"
#include "Estimators/accumulators.h"
#include "Utilities/NewTimer.h"
#include <bitset>

namespace qmcplusplus {

  class WalkerControlBase;
  class EstimatorManager;

  /** Implements branching algorithm for the fixed-node Diffusion Monte Carlo
   *
   * Calculates the Branching Green's function
   * \f[
   * G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])
   * \f]
   * which is used for the weight and mulitplicity of each walker
   * \f[ Weight =  G_{branch} \f]
   * \f[ Mulitplicity =  G_{branch} + \nu \f]
   * and to update the energy offset
   * \f[ E_T = <E_G> - feed \log \left( \frac{P(t)}{P_0} \right) \f]
   * where \f$P(t)\f$ is the current population, \f$P_0\f$ is the 
   * ideal population, and \f$<E_G>\f$ is an estimate of the
   * local energy. 
   */
  class SimpleFixedNodeBranch: public QMCTraits 
  {

    friend class BranchIO;

    /*! enum for booleans 
     * \since 2008-05-05
     */
    enum  
    {
      B_DMC=0, B_DMCSTAGE=1, B_POPCONTROL=2, B_USETAUEFF=3, 
      B_CLEARHISTORY=4,
      B_MODE_MAX=8
    };

    /** booleans to set the branch modes
     * \since 2008-05-05
     *
     * - BranchMode[D_DMC]  1 for dmc, 0 for anything else
     * - BranchMode[D_DMCSTAGE]  1 for main, 0 for wamrup 
     * - BranchMode[D_POPCONTROL] 1 for the standard dmc, 0 for the comb method
     * - BranchMode[B_USETAUEFF]  1 to use taueff accordning to JCP 93, 0 to use tau
     */
    typedef bitset<B_MODE_MAX> BranchModeType;
    BranchModeType BranchMode;

    /*! enum for iParam 
     * \since 2008-05-05
     */
    enum  
    {
      B_WARMUPSTEPS, B_ENERGYUPDATEINTERVAL, B_COUNTER, 
      B_TARGETWALKERS,  B_MAXWALKERS, B_MINWALKERS, B_BRANCHINTERVAL,
      B_IPARAM_MAX
    };

    /** input parameters of integer types 
     * \since 2008-05-05
     *
     * - iParam[B_WARMUPSTEPS] warmup steps, valid when BranchMode[D_DMCSTAGE] == 0
     * - iParam[B_ENERGYUPDATEINTERVAL] frequency of the trial energy updates, default 1 
     *   -- iParam[B_ENERGYUPDATEINTERVAL]  = 1 for the warmup
     * - iParam[B_COUNTER] counter for tracking object state
     * - iParam[B_TARGETWALKERS] target total number of walkers per mpi group
     * - iParam[B_MAXWALKERS]  maximum number of walkers per node
     * - iParam[B_MINWALKERS]  minimum number of walkers per node
     */
    typedef TinyVector<int,B_IPARAM_MAX> IParamType;
    IParamType iParam; 

    /*! enum for vParam */
    enum  
    {
      B_TAU, B_TAUEFF, B_ETRIAL, B_EREF, B_ENOW,
      B_BRANCHMAX, B_BRANCHCUTOFF, B_BRANCHFILTER, B_SIGMA, 
      B_ACC_ENERGY, B_ACC_SAMPLES, 
      B_VPARAM_MAX
    };

    /** controlling parameters of real type
     *
     * Mostly internal
     */
    typedef TinyVector<RealType,B_VPARAM_MAX> VParamType;
    VParamType vParam;

    /** number of remaning steps for a specific tasks
     *
     * set differently for BranchMode[B_DMCSTAGE] 
     */
    int ToDoSteps;
    ////the timestep
    //RealType Tau;
    /////the effective timestep
    //RealType TauEff;
    ///feedback parameter to control the population
    RealType Feedback;
    ///energy offset to control branching
    //RealType Eref;
    /////actual trial energy
    //RealType Etrial;
    ///Feed*log(N)
    RealType logN;
    ///save xml element
    xmlNodePtr myNode;
    ///WalkerController
    WalkerControlBase* WalkerController;
    ///EstimatorManager
    EstimatorManager*  MyEstimator;
    ///a simple accumulator for energy
    accumulator_set<RealType> EnergyHist;
    ///a simple accumulator for energy
    accumulator_set<RealType> R2Accepted;
    ///a simple accumulator for energy
    accumulator_set<RealType> R2Proposed;
    ///histogram of populations
    BlockHistogram<RealType> PopHist;
    ///histogram of populations
    BlockHistogram<RealType> DMCEnergyHist;
    ///root name
    string RootName;
    ///set of parameters
    ParameterSet m_param;
    ///string parameters
    vector<string> sParam;

    public:

    typedef SimpleFixedNodeBranch ThisType;

    //@TODO move these to private
    ///LogJacob
    RealType LogJacobRef;
    ///LogNorm
    vector<RealType> LogNorm;

    ///Constructor
    SimpleFixedNodeBranch(RealType tau, int nideal);

    ///copy constructor
    SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch);

    inline bool phaseChanged(RealType psi0, RealType psi1) const 
    {
      return abs(psi0-psi1) > numeric_limits<RealType>::epsilon();
    }

    /** increment QMCCounter
     *
     * QMCCounter is the number of times any QMC section is processed.
     */
    inline void advanceQMCCounter() { iParam[B_COUNTER]++;}

    /** get the EstimatorManager */
    EstimatorManager* getEstimatorManager()
    {
      return MyEstimator;
    }

    /** set the EstimatorManager
     * @param est estimator created by the first QMCDriver
     * */
    void setEstimatorManager(EstimatorManager* est)
    {
      MyEstimator=est;
    }

    /** initialize  the WalkerController 
     * @param fixW true, if reconfiguration with the fixed number of walkers is used
     */
    void initWalkerController(RealType tau, bool fixW=false);

    /** return the bare branch weight
     *
     * This is equivalent to calling branchWeight(enew,eold,1.0,1.0)
     */
    inline RealType branchWeightBare(RealType enew, RealType eold) const 
    { 
      return std::exp(vParam[B_TAUEFF]*(vParam[B_ETRIAL]-0.5*(enew+eold)));
    }

    /** return the bare branch weight with a filtering using an energy window
     *
     * Cutoff values are set by the variance
     */
    inline RealType branchWeight(RealType enew, RealType eold) const 
    { 
      RealType taueff_=vParam[B_TAUEFF]*0.5;
      RealType x=std::max(vParam[B_EREF]-enew,vParam[B_EREF]-eold);
      if(x>vParam[B_BRANCHMAX])
        taueff_=0.0;
      else if(x>vParam[B_BRANCHCUTOFF])
        taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
      return std::exp(taueff_*(vParam[B_ETRIAL]*2.0-enew-eold));
    }

    /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
     * @param enew new energy
     * @param eold old energy
     * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
     * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
     */
    inline RealType branchWeight(RealType enew, RealType eold, RealType scnew, RealType scold) const
    {
      RealType s1=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-enew)*scnew;
      RealType s0=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-eold)*scold;
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
      RealType s1=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-enew)*scnew;
      RealType s0=(vParam[B_ETRIAL]-vParam[B_EREF])+(vParam[B_EREF]-eold)*scold;
      return std::exp(vParam[B_TAUEFF]*(p*0.5*(s1-s0)+s0));
      //return std::exp(TauEff*(p*0.5*(sp-sq)+sq));
    }

    inline RealType getEref() const { return vParam[B_EREF];}
    inline RealType getTau() const { return vParam[B_TAU];}
    inline RealType getTauEff() const { return vParam[B_TAUEFF];}

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
    void branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones);

    /** restart averaging
     * @param counter Counter to determine the cummulative average will be reset.
     */
    void flush(int counter);

    /** reset the internal parameters */
    void reset();

    bool put(xmlNodePtr cur);

    void write(const string& fname, bool overwrite);
    void read(const string& fname);

    /** create map between the parameter name and variables */
    void registerParameters();

    ///start a run
    void start(const string& froot, bool append);
    ///finalize the simulation
    void finalize();

    private:
    ///default constructor (disabled)
    SimpleFixedNodeBranch(){}

    ///disable use by external users
    //void write(hid_t grp, bool append=false);
    //void read(hid_t grp);

  };

}
#endif
/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

