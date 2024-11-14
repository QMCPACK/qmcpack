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
 *  Some ParticleSet functions use the global Random so we need some helper functions to
 *  avoid interminant test state when multiple tests are run from a single test program.
 */

#ifndef QMCPLUSPLUS_TESTING_GENERATE_RANDOM_PARTICLE_SETS_H
#define QMCPLUSPLUS_TESTING_GENERATE_RANDOM_PARTICLE_SETS_H

#include <vector>
#include <ParticleSet.h>

namespace qmcplusplus
{

namespace testing
{
/**
 *  This function sets particle set positions from a set of Rs or writes out a set of positions for test reproducibility.
 */
std::vector<ParticleSet> generateRandomParticleSets(ParticleSet& pset_target,
                                                    ParticleSet& pset_source,
                                                    std::vector<ParticleSet::ParticlePos>& deterministic_rs,
                                                    int num_psets,
						   bool generate_test_data = false);

} // namespace testing
} // namespace qmcplusplus

#endif
