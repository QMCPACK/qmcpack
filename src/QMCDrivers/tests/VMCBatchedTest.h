//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_VMCBATCHEDTEST_H
#define QMCPLUSPLUS_VMCBATCHEDTEST_H

#include <string>
#include <vector>

namespace qmcplusplus::testing
{
struct VMCBatchedSummary
{
  int walkers                 = 0;
  int proposals               = 0;
  int resolutions             = 0;
  int accepted                = 0;
  int rejected                = 0;
  int mixed_resolutions       = 0;
  int particle_moves          = 0;
  int expected_particle_moves = 0;
  std::vector<std::string> failures;
};

class VMCBatchedTest
{
public:
  static constexpr int steps            = 2;
  static constexpr int substeps         = 1;
  static constexpr int moves_per_walker = steps * substeps;

  static VMCBatchedSummary runAllParticleNoDrift(int walkers_per_crowd = 1);
  static void testAllParticleModeValidation();
};
} // namespace qmcplusplus::testing

#endif
