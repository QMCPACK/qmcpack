//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <cstring>
#include "QMCDrivers/WalkerProperties.h"
#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/HDFWalkerInput_0_4.h"
#include "QMCDrivers/WalkerProperties.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
using WP = WalkerProperties::Indexes;

TEST_CASE("walker", "[particle]")
{

  MCPWalker w(1);
  REQUIRE(w.R.size() == 1);
  w.R[0] = 1.0;

  REQUIRE(w.R[0][0] == Approx(1.0));
}

/** Currently significant amounts of code assumes that the Walker by default 
 *  has "Properties" ending with the LOCALPOTENTIAL element. This opens the door to off by 1
 *  when considering the default size.
 */
TEST_CASE("walker assumptions", "[particle]")
{
  using WP = WalkerProperties::Indexes;
  MCPWalker w1(1);
  REQUIRE(w1.Properties.cols() == WP::NUMPROPERTIES);
}

TEST_CASE("walker HDF read and write", "[particle]")
{
  Communicate* c = OHMMS::Controller;

  MCPWalker w1(1);
  w1.R[0] = 1.0;

  MCPWalker w2(1);
  w2.R[0] = 0.5;

  MCWalkerConfiguration W;


  W.setName("electrons");

  W.create(1);
  W.R[0][0] = 0.0;
  W.R[0][1] = 1.0;
  W.R[0][2] = 2.0;

  // This method sets ownership to false so class does not attempt to
  // free the walker elements.
  W.copyWalkerRefs(&w1, &w2);

  REQUIRE(W.getActiveWalkers() == 2);

  std::vector<int> walker_offset(c->size() + 1);

  walker_offset[0] = 0;
  int offset       = 2;
  for (int i = 0; i < c->size(); i++)
  {
    walker_offset[i + 1] = offset;
    offset += 2;
  }

  W.setWalkerOffsets(walker_offset);

  c->setName("walker_test");
  HDFWalkerOutput hout(W, "this string apparently does nothing", c);
  hout.dump(W, 0);

  c->barrier();

  MCWalkerConfiguration W2;
  W2.setName("electrons");
  W2.create(1);

  HDFVersion version(0, 4);
  HDFWalkerInput_0_4 hinp(W2, c, version);
  bool okay = hinp.read_hdf5("walker_test");
  REQUIRE(okay);

  REQUIRE(W2.getActiveWalkers() == 2);
  for (int i = 0; i < 3; i++)
  {
    REQUIRE(W2[0]->R[0][i] == w1.R[0][i]);
    REQUIRE(W2[1]->R[0][i] == w2.R[0][i]);
  }
}

TEST_CASE("walker buffer add, update, restore", "[particle]")
{
  int num_particles = 4;
  
  UPtrVector<MCPWalker> walkers(2);
  auto createWalker = [num_particles](UPtr<MCPWalker>& walker_ptr) {
    walker_ptr    = std::make_unique<MCPWalker>(num_particles);
    walker_ptr->registerData();
    walker_ptr->DataSet.allocate();
  };
  std::for_each(walkers.begin(), walkers.end(), createWalker);
  
  walkers[0]->Properties(WP::LOGPSI) = 1.2;
  walkers[0]->Properties(WP::SIGN) = 1.3;
  walkers[0]->Properties(WP::UMBRELLAWEIGHT) = 1.4;
  walkers[0]->Properties(WP::LOCALPOTENTIAL) = 1.6;
  walkers[0]->updateBuffer();

  std::memcpy(walkers[1]->DataSet.data(), walkers[0]->DataSet.data(), walkers[0]->DataSet.size());
  walkers[1]->copyFromBuffer();
  CHECK(walkers[1]->Properties(WP::LOGPSI) == Approx(1.2));
  CHECK(walkers[1]->Properties(WP::LOCALPOTENTIAL) == Approx(1.6));
}

} // namespace qmcplusplus
