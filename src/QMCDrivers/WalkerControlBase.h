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
#ifndef QMCPLUSPLUS_WALKER_CONTROL_BASE_H
#define QMCPLUSPLUS_WALKER_CONTROL_BASE_H

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Message/MPIObjectBase.h"
#include "Message/CommOperators.h"
//#include <boost/archive/binary_oarchive.hpp>

namespace qmcplusplus {

  /** Base class to control the walkers for DMC simulations.
   *
   * The virtual functions are implemented for a serial execution with a usual birth/death
   * process. Inherited classes implement other WalkerControl algorithms by implementing
   * branch function.
   * - 2006/10/18
   * -- curData and accumData are added.
   * -- branch function handles averaging of the local energies to apply the weight over the 
   *   current walker sample correctly. 
   * -- removes an extra global reduction.
   *   --- curData should be used for global reductions
   * - 2007/11/07
   * -- added functions to dump dmc histograms
   */
  class WalkerControlBase: public QMCTraits, public MPIObjectBase
  {

    public:
    ///typedef of Walker_t
    typedef MCWalkerConfiguration::Walker_t Walker_t;

    enum {ENERGY_INDEX=0, ENERGY_SQ_INDEX, WALKERSIZE_INDEX, WEIGHT_INDEX, EREF_INDEX, 
      R2ACCEPTED_INDEX, R2PROPOSED_INDEX, LE_MAX};

    ///context id
    IndexType MyContext;
    ///number of contexts
    IndexType NumContexts;
    ///0 is default
    IndexType SwapMode;
    ///minimum number of walkers
    IndexType Nmin;
    ///maximum number of walkers
    IndexType Nmax;
    ///maximum copy per walker
    IndexType MaxCopy;
    ///current number of walkers per processor
    IndexType NumWalkers;
    ///a dummy integer
    IndexType DummyIndex;
    ///trial energy energy
    RealType trialEnergy;
    ///target average energy
    RealType targetAvg;
    ///target average variance
    RealType targetVar;
    ///target sigma to limit fluctuations of the trial energy
    RealType targetSigma;
    ///bound of the energy window
    RealType targetEnergyBound;
    ///current variance
    RealType curVar;
    ///number of particle per node
    vector<int> NumPerNode;
    ///offset of the particle index
    vector<int> OffSet;
    ///offset of the particle index for a fair distribution
    vector<int> FairOffSet;
    
    ///ensenble properties
    MCDataType<RealType> EnsembleProperty;

    ///file to save energy histogram
    ofstream* dmcStream;
    ///archive
    //boost::archive::binary_oarchive *oa;

    ///any accumulated data over a block
    vector<RealType> accumData;
    ///any temporary data
    vector<RealType> curData;
    ///temporary storage for good walkers
    vector<Walker_t*> good_w;
    ///temporary storage for copy counters
    vector<int> ncopy_w;

    /** default constructor
     *
     * Set the SwapMode to zero so that instantiation can be done
     */
    WalkerControlBase(Communicate* c);

    /** empty destructor to clean up the derived classes */
    virtual ~WalkerControlBase();

    /** start controller */
    void start();

    /** take averages and writes to a file */
    void measureProperties(int iter);

    /** set the trial energy
     */
    inline void setTrialEnergy(RealType et)
    {
      trialEnergy=et;
    }

    /** return a value accumulated during a block
     * @param i index of the data
     *
     * use enum for i, see DMCEnergyEstimator
     */
    inline RealType getValue(int i) {
      return accumData[i];
    }

    /** return a current value 
     * @param i index of the data
     *
     * use enum for i, see DMCEnergyEstimator
     */
    inline RealType getCurrentValue(int i) {
      return curData[i];
    }

    /** set the target average and variance
     */
    inline void setEnergyAndVariance(RealType e, RealType v) {
      trialEnergy=e;
      targetAvg=e;
      targetVar=v;
    }

    /** sort Walkers between good and bad and prepare branching
     */
    void sortWalkers(MCWalkerConfiguration& W);
    /** copy good walkers to W
     */
    int copyWalkers(MCWalkerConfiguration& W);

    /** reset to accumulate data */
    virtual void reset();
    
    /** perform branch and swap walkers as required */
    virtual int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

    virtual RealType getFeedBackParameter(int ngen, RealType tau) {
      return 1.0/(static_cast<RealType>(ngen)*tau);
    }

    bool put(xmlNodePtr cur);

  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

