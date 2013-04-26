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
#ifndef QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROL_H
#define QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct WalkerReconfiguration: public WalkerControlBase
{

  //random number [0,1)
  RealType UnitZeta;

  vector<int>      IndexCopy;
  //weight per walker
  vector<RealType> wConf;
  //comb
  vector<RealType> Zeta;
  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerReconfiguration(Communicate* c);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }

  /** return the surviving Walkers
   */
  int getIndexPermutation(MCWalkerConfiguration& W);
  int shuffleIndex(int nw);
};
}
#endif
/***************************************************************************
 * $RCSfile: WalkerReconfiguration.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

