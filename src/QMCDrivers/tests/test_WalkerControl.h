//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TEST_WALKERCONTROLMPI_H
#define QMCPLUSPLUS_TEST_WALKERCONTROLMPI_H

#include "QMCDrivers/DMC/WalkerControl.h"
#include "QMCDrivers/tests/SetupPools.h"

namespace qmcplusplus
{
namespace testing
{
/** This class is friends with WalkerControl
 *  It allows the test_WalkerControl.cpp program to be called with variable numbers of ranks and simplifies
 *  writing additional test cases, while maintain ecapsulation of WalkerControls internals.
 *
 *  Since the public of WalkerControl are massive multi object state transformations most/all unit tests
 *  are of private methods.
 *
 *  There is one per rank.
 *  \Todo Once there is only one driver type rename
 */
class UnifiedDriverWalkerControlMPITest
{
public:
  /** Initially this test class creates a population with 1 walker.
   */
  UnifiedDriverWalkerControlMPITest();

  /** Wrapper to unit test determineNewWalkerPopulation
   */
  void testNewDistribution(std::vector<int>& initial_num_per_rank, std::vector<int>& minus, std::vector<int>& plus);

  /** This test function manipulates walker multiplicity and then checks that after "population control"
   *  that the number of walkers per rank is correct. Multiplicity in normal branching steps is a set
   *  walker->Weight + rng() before swapping and amplification of walker's with multiplicity >= 2
   *
   *  \param[in]  walker_multiplicity_before   walker multiplicity by rank
   *  \param[in]  rank_counts_after            walker count by rank
   */
  void testPopulationDiff(std::vector<int>& walker_multiplicity_total, std::vector<int>& rank_counts_after);

  /** This function is intended to be called after a call to testPopulationDiff to test that correct walker_ids and
   *  parent_ids are set.
   *  \param[in] walker_ids   the expected set of walker ids on each rank
   *  \param[in] parent_ids   the expected set of parent ids on each rank.
   */
  void testWalkerIDs(std::vector<std::vector<int>> walker_ids_after, std::vector<std::vector<int>> parent_ids_after);

  int getRank() const { return dpools_.comm->rank(); }
  int getNumRanks() const { return dpools_.comm->size(); }

private:
  void reportWalkersPerRank(Communicate* c, MCPopulation& pop);

  SetupPools dpools_;
  WalkerConfigurations walker_confs;
  UPtr<MCPopulation> pop_;
  WalkerControl wc_;
};
} // namespace testing
} // namespace qmcplusplus

#endif /* QMCPLUSPLUSTEST_WALKERCONTROLMPI_H */
