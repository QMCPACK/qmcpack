//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/NewTimer.h"
#include <stdio.h>
#include <string>
#include <vector>

namespace qmcplusplus {

TEST_CASE("test_timer", "[utilities]")
{
  NewTimer t1("timer1");
  TimerManager.addTimer(&t1);
#if ENABLE_TIMER
  t1.start();
  REQUIRE(TimerManager.current_timer() == &t1);
  t1.stop();
  REQUIRE(TimerManager.current_timer() == NULL);
#endif
}

}
