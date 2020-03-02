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


#include "Particle/SampleStack.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("SampleStack", "")
{
  OHMMS::Controller->initialize(0, NULL);

  SampleStack samples;

  const int total_num = 2; // number of particles
  samples.setTotalNum(total_num);

  samples.setMaxSamples(8);
  REQUIRE(samples.getMaxSamples() == 8);
  REQUIRE(samples.getNumSamples() == 0);

  // reserve storage
  samples.setNumSamples(7);
  REQUIRE(samples.getNumSamples() == 0);
  REQUIRE(samples.getMaxSamples() == 7);

  SampleStack::WalkerList_t walker_list;

  // Check adding an empty list doesn't break
  samples.saveEnsemble(walker_list.begin(), walker_list.end());
  REQUIRE(samples.getNumSamples() == 0);


  // Add size one list
  walker_list.push_back(new SampleStack::Walker_t(total_num));
  walker_list[0]->R[0][0] = 1.1;
  samples.saveEnsemble(walker_list.begin(), walker_list.end());
  REQUIRE(samples.getNumSamples() == 1);

  SampleStack::Walker_t w1;
  samples.getSample(0, w1);
  REQUIRE(w1.R[0][0] == Approx(1.1));

  // Should test that more members of the Walker are saved correctly
}


} // namespace qmcplusplus
