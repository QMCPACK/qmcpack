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

  /** Class to handle walker controls
   *
   * Base class to handle serial mode with branching only
   */
  struct WalkerControlBase: public QMCTraits {

    ///typedef of Walker_t
    typedef MCWalkerConfiguration::Walker_t Walker_t;

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
    ///Global energy/weight
    TinyVector<RealType,2> gEavgWgt;

    ///temporary storage for good walkers
    vector<Walker_t*> good_w;

    ///temporary storage for copy counters
    vector<int> ncopy_w;

    /** default constructor
     *
     * Set the SwapMode to zero so that instantiation can be done
     */
    WalkerControlBase(): SwapMode(0), Nmin(1), Nmax(10), MaxCopy(10) {}

    /** empty destructor to clean up the derived classes */
    virtual ~WalkerControlBase() {}
    
    /** sort Walkers between good and bad and prepare branching
     */
    void sortWalkers(MCWalkerConfiguration& W);
    /** copy good walkers to W
     */
    int copyWalkers(MCWalkerConfiguration& W);

    /** perform branch and swap walkers as required */
    virtual int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

    virtual RealType getFeedBackParameter(int ngen, RealType tau) {
      return 1.0/(static_cast<RealType>(ngen)*tau);
    }

#if defined(HAVE_MPI)

    /** collect the energies */
    inline RealType average(RealType eavg, RealType wgt) {
      gEavgWgt[0]=eavg;
      gEavgWgt[1]=wgt;
      gsum(gEavgWgt,0);
      return gEavgWgt[0]/gEavgWgt[1];
    }
#else

    /** collect the energies */
    inline RealType average(RealType eavg, RealType wgt) {
      gEavgWgt[0]=eavg; gEavgWgt[1]=wgt;
      return eavg/wgt;
    }
#endif
  };

  /** function to create WalkerControlBase or its derived class
   * @param swapmode in/out indicator to determine which controller will be created
   * @param nideal ideal number of walkers
   * @param nmax maximum number of walkers
   * @param nmin minimum number of walkers
   * @param wc pointer to current WalkerControlBase object
   * @return WalkerControlBase*
   *
   * When wc is the same as the requested controller object, only reset the
   * internal values and return wc itself.
   */
  WalkerControlBase* CreateWalkerController(bool reconfig, 
      int& swapmode, int nideal, int nmax, int nmin, WalkerControlBase* wc);
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

