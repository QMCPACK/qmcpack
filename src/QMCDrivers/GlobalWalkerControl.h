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
#ifndef OHMMS_QMC_GLOBAL_WALKER_CONTROL_H
#define OHMMS_QMC_GLOBAL_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace ohmmsqmc {

  /** Class to handle walker controls with simple global sum
   *
   * Base class to handle serial mode with branching only
   */
  struct GlobalWalkerControl: public WalkerControlBase {

    int NumSwaps;
    int MyContext;
    int NumContexts;
    int Cur_max;
    int Cur_min;
    int Cur_pop;
    vector<int> NumPerNode;
    vector<int> OffSet;
    vector<int> FairOffSet;
    /** default constructor
     *
     * Set the SwapMode to zero so that instantiation can be done
     */
    GlobalWalkerControl();

    /** perform branch and swap walkers as required */
    int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

    /** collect the energies */
    inline void collect(TinyVector<RealType,2>& eavgwgt) {
      gsum(eavgwgt,0);
    }

    void swapWalkersAsync(MCWalkerConfiguration& W);
    void swapWalkersBlocked(MCWalkerConfiguration& W);
    void swapWalkersMap(MCWalkerConfiguration& W);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

