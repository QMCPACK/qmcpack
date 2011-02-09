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
#ifndef QMCPLUSPLUS_WALKER_CONTROL_MPI_H
#define QMCPLUSPLUS_WALKER_CONTROL_MPI_H

#include "QMCDrivers/WalkerControlBase.h"


namespace qmcplusplus 
{

  class NewTimer;

  /** Class to handle walker controls with simple global sum
   *
   * Base class to handle serial mode with branching only
   */
  struct WalkerControlMPI: public WalkerControlBase 
  {
    int Cur_pop;
    int Cur_max;
    int Cur_min;
    vector<NewTimer*> myTimers;
    /** default constructor
     *
     * Set the SwapMode to zero so that instantiation can be done
     */
    WalkerControlMPI(Communicate* c=0);

    /** perform branch and swap walkers as required */
    int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

    void swapWalkersSimple(MCWalkerConfiguration& W);

    //old implementations
    void swapWalkersAsync(MCWalkerConfiguration& W);
    void swapWalkersBlocked(MCWalkerConfiguration& W);
    void swapWalkersMap(MCWalkerConfiguration& W);
  };
}
#endif
/***************************************************************************
 * $RCSfile: WalkerControlMPI.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

