//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020  QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Particle/MCSample.h"
#include "Particle/SampleStack.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("SampleStack", "[particle]")
{
  SampleStack samples;

  const int total_num = 2; // number of particles

  // reserve storage
  int nranks = 2;
  samples.setMaxSamples(8, nranks);
  REQUIRE(samples.getMaxSamples() == 8);
  REQUIRE(samples.getNumSamples() == 0);
  REQUIRE(samples.getGlobalNumSamples() == 16);

  // increase storage
  samples.setMaxSamples(10);
  REQUIRE(samples.getMaxSamples() == 10);
  REQUIRE(samples.getNumSamples() == 0);

  using Walker_t     = ParticleSet::Walker_t;
  using WalkerList_t = std::vector<std::unique_ptr<Walker_t>>;

  WalkerList_t walker_list;

  // Add size one list
  walker_list.push_back(std::make_unique<Walker_t>(total_num));
  walker_list[0]->R[0][0] = 1.1;
  for (auto& wi : walker_list)
  {
    samples.appendSample(MCSample(*wi));
  }
  REQUIRE(samples.getNumSamples() == 1);

  Walker_t w1;
  samples.getSample(0).convertToWalker(w1);
  CHECK(w1.R[0][0] == Approx(1.1));

  // Should test that more members of the Walker are saved correctly

  samples.resetSampleCount();
  REQUIRE(samples.getNumSamples() == 0);

  // clear storage
  samples.clearEnsemble();
  REQUIRE(samples.getNumSamples() == 0);
}


} // namespace qmcplusplus
