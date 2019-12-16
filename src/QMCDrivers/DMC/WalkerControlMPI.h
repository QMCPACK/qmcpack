//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
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


namespace qmcplusplus
{
class NewTimer;
struct WalkerControlMPITest;

namespace testing
{
struct UnifiedDriverWalkerControlMPITest;
}

struct WalkerMessage
{
  // To deal with duplicate send optimization
  RefVector<WalkerControlBase::MCPWalker> walker;
  // i.e. MPI rank
  const int source_rank;
  const int target_rank;
  int multiplicity = 1;
  size_t byteSize = 0;
  WalkerMessage(WalkerControlBase::MCPWalker& walk, const int source, const int target) : source_rank(source), target_rank(target){ walker.push_back(walk); };
};

inline bool operator==(const WalkerMessage& A, const WalkerMessage& B)
{
  // since all the walker references are to one queue of unique_ptrs in MCPopulation
  // I think this is ok... If not the DataSet.data()?
  return (&A.walker == &B.walker) && (A.target_rank == B.target_rank);
}


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
  ///Number of walkers sent during the exchange
  // Is this persistent state for any reason other than we keep zeroing curData
  // defensively?
  IndexType NumWalkersSent;

  /** default constructor
   *
   * Set the SwapMode to zero so that instantiation can be done
   * comm can not be null it is not checked.
   */
  WalkerControlMPI(Communicate* comm);

  /** creates the distribution plan
   *
   *  populates the minus and plus vectors they contain 1 copy of a partition index 
   *  for each adjustment in population to the context.
   *  \todo fix this word salad
   *
   *  \param[in] Cur_pop current population
   *  \param[in] NumContexts number of MPI processes
   *  \param[in] MyContext my MPI rank
   *  \param[in] NumPerNode current walkers per node
   *  \param[out] FairOffSet running population count at each partition boundary
   *  \param[out] minus list of partition indexes one occurance for each walker removed
   *  \param[out] plus list of partition indexes one occurance for each walker added
   */
  static void determineNewWalkerPopulation(int Cur_pop,
                                           int NumContexts,
                                           int MyContext,
                                           const std::vector<int>& NumPerNode,
                                           std::vector<int>& FairOffSet,
                                           std::vector<int>& minus,
                                           std::vector<int>& plus);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCWalkerConfiguration& W, FullPrecRealType trigger);

  /** perform branch and swap walkers as required */
  int branch(int iter, MCPopulation& pop, FullPrecRealType trigger);

  //current implementations
  void swapWalkersSimple(MCWalkerConfiguration& W);

  void swapWalkersSimple(MCPopulation& pop, PopulationAdjustment& adjust);

  // Testing wrappers
  friend WalkerControlMPITest;
  friend testing::UnifiedDriverWalkerControlMPITest;
};

} // namespace qmcplusplus
#endif
