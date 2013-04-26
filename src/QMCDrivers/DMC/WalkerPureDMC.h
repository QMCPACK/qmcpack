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
#ifndef QMCPLUSPLUS_PUREDMC_WALKER_CONTROL_H
#define QMCPLUSPLUS_PUREDMC_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct WalkerPureDMC: public WalkerControlBase
{

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerPureDMC(Communicate* c):WalkerControlBase(c) {};

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger)
  {
    return doNotBranch(iter, W);
  };

  /** return 0.0 to disable feedback method */
  RealType getFeedBackParameter(int ngen, RealType tau)
  {
    return 0.0;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile: WalkerReconfiguration.h,v $   $Author: jnkim $
 * $Revision: 2284 $   $Date: 2007-11-07 10:29:32 -0600 (Wed, 07 Nov 2007) $
 * $Id: WalkerReconfiguration.h 2284 2007-11-07 16:29:32Z jnkim $
 ***************************************************************************/

