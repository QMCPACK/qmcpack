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
#include "QMCDrivers/WalkerElementsRef.h"
#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
struct WalkerControlMPITest;

namespace testing
{
class UnifiedDriverWalkerControlMPITest;
}

/** Datastruct of a WalkerMessage
 */
struct WalkerMessage
{
  WalkerElementsRef walker_elements;
  // i.e. MPI rank
  const int source_rank;
  const int target_rank;
  WalkerMessage(WalkerElementsRef w_elem, const int source, const int target)
      : walker_elements(w_elem), source_rank(source), target_rank(target)
  {}
};

/** needed for repressing duplicate messages, this optimization is currently unimplemented.
 */
/* inline bool isRepeatedMessage(const WalkerMessage& A, const WalkerMessage& B) */
/* { */
/*   // since all the walker references are to one queue of unique_ptrs in MCPopulation */
/*   return (&A.walker == &B.walker) && (A.target_rank == B.target_rank); */
/* } */

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
   *  \param[out] minus list of partition indexes one occurance for each walker removed
   *  \param[out] plus list of partition indexes one occurance for each walker added
   */
  static void determineNewWalkerPopulation(int cur_pop,
                                           int num_contexts,
                                           int my_context,
                                           std::vector<int>& num_per_rank,
                                           std::vector<int>& fair_offset,
                                           std::vector<int>& minus,
                                           std::vector<int>& plus);

  /** legacy: perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger);

  /** unified driver: perform branch and swap walkers as required */
  FullPrecRealType branch(int iter, MCPopulation& pop);

  /** legacy: swap implementation
   */
  void swapWalkersSimple(MCWalkerConfiguration& W);

  /** unified: swap implementation
   *
   * Walkers are transfered between ranks in order to reach a difference of no more than 1 walker.
   * \param[inout] pops rankscope population
   * \param[inout] adjust population adjustment, it's updated for so on rank adjustments can occur
   * \param[inout]    num_per_rank number of walkers per rank expanded by multiplicity.
   */
  int swapWalkersSimple(MCPopulation& pop, PopulationAdjustment& adjust, std::vector<IndexType>& num_per_rank);

  // Testing wrappers
  friend WalkerControlMPITest;
  friend testing::UnifiedDriverWalkerControlMPITest;
};

} // namespace qmcplusplus
#endif
