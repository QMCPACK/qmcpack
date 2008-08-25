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
#ifndef QMCPLUSPLUS_WALKER_CONTROL_FACTORY_H
#define QMCPLUSPLUS_WALKER_CONTROL_FACTORY_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus {

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
      int& swapmode, int nideal, int nmax, int nmin, 
      WalkerControlBase* wc,
      Communicate* comm);

  /** function to create WalkerControlBase
   * @param current number of walkers
   */
  WalkerControlBase* createWalkerController(int nwtot,
      Communicate* comm, xmlNodePtr cur, bool reconfig=false);
}
#endif
/***************************************************************************
 * $RCSfile: WalkerControlFactory.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

