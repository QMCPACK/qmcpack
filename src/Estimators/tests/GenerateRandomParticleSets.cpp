//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  Several ParticleSet functions use the global Random so we have to avoid the normal
 *  sequence of particleset state transforms and set particle positions explicitly.
 *  This function does that as well as facilitates writing out a set of positions for test reproducibility.
 */

#include "GenerateRandomParticleSets.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

namespace qmcplusplus
{
namespace testing
{
std::vector<ParticleSet> generateRandomParticleSets(ParticleSet& pset_target,
                                                    ParticleSet& pset_source,
                                                    std::vector<ParticleSet::ParticlePos>& deterministic_rs,
                                                    int num_psets,
                                                    bool generate_test_data)
{
  int nwalkers = num_psets;
  std::vector<ParticleSet> psets(num_psets, pset_target);
  if (generate_test_data)
  {
    std::cout << "Initialize psets with:\n{";
    std::vector<ParticleSet> psets;
    for (int iw = 0; iw < nwalkers; ++iw)
    {
      psets.emplace_back(pset_target);
      psets.back().randomizeFromSource(pset_source);
      std::cout << "{";
      for (auto r : psets.back().R)
        std::cout << NativePrint(r) << ",";
      std::cout << "},\n";
      psets.back().update();
    }
    std::cout << "}\n";
  }
  else
  {
    for (int iw = 0; iw < nwalkers; ++iw) {
      psets[iw].R = deterministic_rs[iw];
      psets[iw].update();
    }
  }
  return psets;
}

} // namespace testing
} // namespace qmcplusplus
