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
  UnifiedDriverWalkerControlMPITest();
  void testMultiplicity(std::vector<int>& rank_counts_expanded, std::vector<int>& rank_counts_after);
  void testPopulationDiff(std::vector<int>& rank_counts_before, std::vector<int>& rank_counts_after);
  void makeValidWalkers();
  static void testNewDistribution(std::vector<int>& minus, std::vector<int>& plus);

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
