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
    /*! enum for booleans 
     * \since 2008-05-05
     */
    enum  
    {
      B_DMCSTAGE, B_POPCONTROL, B_USETAUEFF, B_CLEARHISTORY
    };

    /** booleans to set the branch modes
     * \since 2008-05-05
     *
     * - BranchMode[D_DMCSTAGE]  1 for main, 0 for wamrup 
     * - BranchMode[D_POPCONTROL] 1 for the standard dmc, 0 for the comb method
     * - BranchMode[B_USETAUEFF]  1 to use taueff accordning to JCP 93, 0 to use tau
     */
    bitset<4> BranchMode;

    /*! enum for iParam 
     * \since 2008-05-05
     */
    enum  
    {
      B_WARMUPSTEPS, B_ENERGYUPDATEINTERVAL, B_COUNTER, 
      B_TARGETWALKERS,  B_MAXWALKERS, B_MINWALKERS, B_BRANCHINTERVAL 
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
    TinyVector<int,8> iParam; 

    /*! enum for vParam */
    enum  
    {
      B_BRANCHMAX, B_BRANCHCUTOFF, B_BRANCHFILTER, B_ENERGYWINDOW
    };

    /** controlling parameters of real type
     *
     * Mostly internal
     */
    TinyVector<RealType,4> vParam;

    /** number of remaning steps for a specific tasks
     *
     * set differently for BranchMode[B_DMCSTAGE] 
     */
    int ToDoSteps;
    //the timestep
    RealType Tau;
    ///the effective timestep
    RealType TauEff;
    ///feedback parameter to control the population
    RealType Feedback;
    ///energy offset to control branching
    RealType Eref;
    ///actual trial energy
    RealType Etrial;
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

    /**  Calculates the Branching Green's function
     * @param tau effective time step
     * @param emixed mixed energy \f$(E_L(R)+E_L(R'))/2\f$
     * @param reject rejection probability
     * @return \f$G_{branch}\f$
     *
     * \f[G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])\f]
     * @note Use the rejection probability \f$q\f$ to limit \f$G_{branch}\f$
     * \f[ G_{branch} = \min \left(\frac{1}{2q},G_{branch}\right). \f]
     */
    inline RealType branchGF(RealType tau, RealType emixed, RealType reject) const { 
      //return min(0.5/(reject+1e-12),exp(-tau*(emix-Etrial)));
      return std::exp(-tau*(emixed-Etrial));
    }

    //inline RealType branchWeight(RealType tau_, RealType enew, RealType eold) const 
    //{ 
    //  RealType taueff_=tau_*0.5;
    //  RealType x=std::max(Eref-enew,Eref-eold);
    //  if(x>vParam[B_BRANCHMAX])
    //    taueff_=0.0;
    //  else if(x>vParam[B_BRANCHCUTOFF])
    //    taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
    //  return std::exp(taueff_*(Etrial*2.0-enew-eold));
    //}

    inline RealType branchWeight(RealType enew, RealType eold) const 
    { 
      RealType taueff_=TauEff*0.5;
      RealType x=std::max(Eref-enew,Eref-eold);
      if(x>vParam[B_BRANCHMAX])
        taueff_=0.0;
      else if(x>vParam[B_BRANCHCUTOFF])
        taueff_ *=(1.0-(x-vParam[B_BRANCHCUTOFF])*vParam[B_BRANCHFILTER]);
      return std::exp(taueff_*(Etrial*2.0-enew-eold));
    }

    /** return the branch weight according to JCP1993 Umrigar et al. Appendix A p=1, q=0
     * @param enew new energy
     * @param eold old energy
     * @param scnew  \f$ V_{sc}(R_{new})/V(R_{new}) \f$
     * @param scold  \f$ V_{sc}(R_{old})/V(R_{old}) \f$
     */
    inline RealType branchWeight(RealType enew, RealType eold, RealType scnew, RealType scold) const
    {
      return std::exp(TauEff*(Etrial-Eref+0.5*((Eref-enew)*scnew+(Eref-eold)*scold)));
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
      RealType sp=(Etrial-Eref)+(Eref-enew)*scnew;
      RealType sq=(Etrial-Eref)+(Eref-eold)*scold;
      return std::exp(TauEff*(p*0.5*(sp-sq)+sq));
    }

    inline RealType getEref() const { return Eref;}
    inline RealType getTau() const { return Tau;}
    inline RealType getTauEff() const { return TauEff;}

    inline void setTrialEnergy(RealType etot, RealType wtot=1.0) 
    {
      Eref=Etrial=etot/wtot;
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

