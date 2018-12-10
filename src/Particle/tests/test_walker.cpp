//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/HDFWalkerInput_0_4.h"



#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{


typedef Walker<QMCTraits, PtclOnLatticeTraits> Walker_t;

TEST_CASE("walker", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);

  Walker_t w(1);
  REQUIRE(w.R.size() == 1);
  w.R[0] = 1.0;

  REQUIRE(w.R[0][0] == Approx(1.0));
}


TEST_CASE("walker HDF read and write", "[particle]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  Walker_t w1(1);
  w1.R[0] = 1.0;

  Walker_t w2(1);
  w2.R[0] = 0.5;

  MCWalkerConfiguration W;


  W.setName("electrons");

  W.create(1);
  W.R[0][0] = 0.0;
  W.R[0][1] = 1.0;
  W.R[0][2] = 2.0;

  // This method sets ownership to false so class does not attempt to
  // free the walker elements.
  W.copyWalkerRefs(&w1,&w2);

  REQUIRE(W.getActiveWalkers() == 2);

  std::vector<int> walker_offset(c->size()+1);

  walker_offset[0] = 0;
  int offset = 2;
  for (int i = 0; i < c->size(); i++)
  {
    walker_offset[i+1] = offset;
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

  HDFVersion version(0,4);
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

}
