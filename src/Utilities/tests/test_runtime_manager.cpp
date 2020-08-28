//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/RunTimeManager.h"
#include <stdio.h>
#include <string>
#include <vector>

namespace qmcplusplus
{
TEST_CASE("test_runtime_manager", "[utilities]")
{
  // Use a local version rather than the global timer_manager, otherwise
  //  changes will persist from test to test.
  RunTimeManagerClass<fake_cpu_clock> rm;
  double e = rm.elapsed();
  REQUIRE(e == Approx(1.0));
}

TEST_CASE("test_loop_timer", "[utilities]")
{
  LoopTimer<fake_cpu_clock> loop;
  double it_time = loop.get_time_per_iteration();
  REQUIRE(it_time == Approx(0.0));

  loop.start();
  loop.stop();
  it_time = loop.get_time_per_iteration();
  REQUIRE(it_time == Approx(1.0));

  fake_cpu_clock::fake_cpu_clock_increment = 2.0;
  loop.start();
  loop.stop();
  it_time = loop.get_time_per_iteration();
  REQUIRE(it_time == Approx(1.5)); // 2 iterations
  // restore value
  fake_cpu_clock::fake_cpu_clock_increment = 1.0;
}

TEST_CASE("test_loop_control", "[utilities]")
{
  // fake clock advances every time cpu_clock is called
  LoopTimer<fake_cpu_clock> loop;
  int max_cpu_secs = 9;
  RunTimeManagerClass<fake_cpu_clock> rm; // fake clock = 1
  RunTimeControl<fake_cpu_clock> rc(rm, max_cpu_secs);
  rc.runtime_padding(1.0);
  bool enough_time = rc.enough_time_for_next_iteration(loop); // fake clock = 2
  REQUIRE(enough_time);
  loop.start();                                          // fake clock = 3
  loop.stop();                                           // fake clock = 4
  enough_time = rc.enough_time_for_next_iteration(loop); // fake clock = 5
      // estimated time with margin and padding =  1.0 sec/it * 1.1 + 1.0 (pad) = 2.1
      // remaining = 9 - 4.0 = 5.0  enough time for another loop.
  REQUIRE(enough_time);

  loop.start();                                          // fake clock = 6
  loop.stop();                                           // fake clock = 7
  enough_time = rc.enough_time_for_next_iteration(loop); // fake clock = 8
      // estimated time with margin and padding = 1.0 sec/it * 1.1 + 1.0 = 2.1
      // remaining = 9 - 8.0 = 1.0  not enough time for another loop
  REQUIRE(!enough_time);

  std::string msg = rc.time_limit_message("QMC", 2);
  REQUIRE(msg.size() > 0);
}


} // namespace qmcplusplus
