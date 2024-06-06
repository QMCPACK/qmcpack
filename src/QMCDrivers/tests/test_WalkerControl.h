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
/** Once there is only one driver type rename
 */
class UnifiedDriverWalkerControlMPITest
{
public:
  /** Initially this test class creates a population with 1 walker so in a multirank test you have 1 walker per rank.
   */
  UnifiedDriverWalkerControlMPITest();
  void testPopulationDiff(std::vector<int>& rank_counts_before, std::vector<int>& rank_counts_after);
  /** This test function alters the population respecting what was there before the call
   */
  void testPopulationDiffAdditionalStage(std::vector<int>& rank_counts_before, std::vector<int>& rank_counts_after);
  void testWalkerIDs(std::vector<std::vector<int>> walker_ids_after, std::vector<std::vector<int>> parent_ids_after);
//  void makeValidWalkers();
  void testNewDistribution(std::vector<int>& minus, std::vector<int>& plus);
  int getRank() { return dpools_.comm->rank(); }
  int getNumRanks() { return dpools_.comm->size(); }
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
