//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "QMCDrivers/MoveContext.h"

namespace qmcplusplus
{
using WalkerHere = Walker<QMCTraits, PtclOnLatticeTraits>;

TEST_CASE("MoveContext::loadWalker", "[particle]")
{
  WalkerHere walker(1);
  walker.R[0] = TinyVector<double, 3>(1.0, 0.0, 0.0);
  MoveContext move_context;

  move_context.loadWalker(walker);
  REQUIRE(move_context.get_positions() == TinyVector<double, 3>(1.0, 0.0, 0.0));
  
}

}
