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
    
    


#ifndef QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROLMPI_H
#define QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROLMPI_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct WalkerReconfigurationMPI: public WalkerControlBase
{

  ///total number of walkers
  int TotalWalkers;
  ///starting index of the local walkers
  int FirstWalker;
  ///ending index of the local walkers
  int LastWalker;
  ///random number [0,1)
  RealType UnitZeta;
  ///random number [0,1)/number of walkers
  RealType DeltaStep;
  ///1/(total number of walkers)
  RealType nwInv;
  ///the number of extra/missing walkers
  std::vector<IndexType> dN;
  //weight per walker
  std::vector<RealType> wConf;
  //accumulated weight [0,ip) for each ip
  std::vector<RealType> wOffset;
  //local sum of the weights for each ip
  std::vector<RealType> wSum;
  //comb
  //vector<RealType> Zeta;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerReconfigurationMPI(Communicate* c=0);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, RealType trigger);

  /** return the surviving Walkers
   */
  int swapWalkers(MCWalkerConfiguration& W);


  /** send the extra walkers to other node
   * @param plus local indices of the walkers to be sent
   */
  void sendWalkers(MCWalkerConfiguration& W, const std::vector<IndexType>& plus);

  /** send the missing walkers from other node
   * @param minus local indices of the walkers to be copied
   */
  void recvWalkers(MCWalkerConfiguration& W, const std::vector<IndexType>& minus);

};
}
#endif

