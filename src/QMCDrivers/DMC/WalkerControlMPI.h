//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  TimerList_t myTimers;
  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerControlMPI(Communicate* c=0);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  //current implementations
  void swapWalkersSimple(MCWalkerConfiguration& W);

};
}
#endif

