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


#ifndef QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROL_H
#define QMCPLUSPLUS_RECONFIGURATION_WALKER_CONTROL_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{
/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct WalkerReconfiguration : public WalkerControlBase
{
  //random number [0,1)
  FullPrecRealType UnitZeta;

  std::vector<int> IndexCopy;
  //weight per walker
  std::vector<FullPrecRealType> wConf;
  //comb
  std::vector<FullPrecRealType> Zeta;
  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   */
  WalkerReconfiguration(Communicate* c);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger) override;

  /** return the surviving Walkers
   */
  int getIndexPermutation(MCWalkerConfiguration& W);
  int shuffleIndex(int nw);
};
} // namespace qmcplusplus
#endif
