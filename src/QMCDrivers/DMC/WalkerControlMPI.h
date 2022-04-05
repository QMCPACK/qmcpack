//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKER_CONTROL_MPI_H
#define QMCPLUSPLUS_WALKER_CONTROL_MPI_H

#include "QMCDrivers/WalkerControlBase.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
struct WalkerControlMPITest;

/** Class to handle walker controls with simple global sum
 *
 * Base class to handle serial mode with branching only
 */
struct WalkerControlMPI : public WalkerControlBase
{
  int Cur_pop;
  int Cur_max;
  int Cur_min;
  TimerList_t myTimers;
  // Number of walkers sent during the exchange
  // This semi-persistent state is only here because we keep zeroing curData
  // defensively?
  IndexType NumWalkersSent;

  /** default constructor
   *
   * \param[in] comm can not be null it is not checked.
   */
  WalkerControlMPI(Communicate* comm);

  /** creates the distribution plan
   *
   *  populates the minus and plus vectors they contain 1 copy of a partition index 
   *  for each adjustment in population to the context.
   *  \todo fix this argument salad
   *
   *  \param[in] cur_pop_pop population taking multiplicity into account
   *  \param[in] num_contexts number of MPI processes
   *  \param[in] my_context i.e this processes MPI rank
   *  \param[in/out] num_per_rank as if all walkers were copied out to multiplicity
   *  \param[out] fair_offset running population count at each partition boundary
   *  \param[out] minus list of partition indexes one occurrence for each walker removed
   *  \param[out] plus list of partition indexes one occurrence for each walker added
   */
  static void determineNewWalkerPopulation(int cur_pop,
                                           int num_contexts,
                                           int my_context,
                                           std::vector<int>& num_per_rank,
                                           std::vector<int>& fair_offset,
                                           std::vector<int>& minus,
                                           std::vector<int>& plus);

  /** legacy: perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger) override;

  /** legacy: swap implementation
   */
  void swapWalkersSimple(MCWalkerConfiguration& W);

  // Testing wrappers
  friend WalkerControlMPITest;
};

} // namespace qmcplusplus
#endif
