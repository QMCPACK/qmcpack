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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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
#include "Message/Communicate.h"
#include "Message/CommOperators.h"

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
   */
  struct WalkerControlBase: public QMCTraits {

    ///typedef of Walker_t
    typedef MCWalkerConfiguration::Walker_t Walker_t;

    enum {ENERGY_INDEX=0, ENERGY_SQ_INDEX, WALKERSIZE_INDEX, WEIGHT_INDEX, EREF_INDEX, LE_MAX};

    ///0 is default
    int SwapMode;
    ///minimum number of walkers
    IndexType Nmin;
    ///maximum number of walkers
    IndexType Nmax;
    ///maximum copy per walker
    IndexType MaxCopy;
    ///current number of walkers per processor
    IndexType NumWalkers;
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
    WalkerControlBase();

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

    /** sort Walkers between good and bad and prepare branching
     */
    void sortWalkers(MCWalkerConfiguration& W);
    /** copy good walkers to W
     */
    int copyWalkers(MCWalkerConfiguration& W);

    /** empty destructor to clean up the derived classes */
    virtual ~WalkerControlBase() {}

    /** reset to accumulate data */
    virtual void reset();
    
    /** perform branch and swap walkers as required */
    virtual int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

    virtual RealType getFeedBackParameter(int ngen, RealType tau) {
      return 1.0/(static_cast<RealType>(ngen)*tau);
    }
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

