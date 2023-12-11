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

using namespace std::chrono_literals;

namespace qmcplusplus
{

// Convert duration input type to nanosecond duration
template<typename T>
FakeChronoClock::duration convert_to_ns(T in)
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(in);
}

TEST_CASE("test_runtime_manager", "[utilities]")
{
  // Use a local version rather than the global timer_manager, otherwise
  //  changes will persist from test to test.
  RunTimeManager<FakeChronoClock> rm;
  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.0s);
  double e                                     = rm.elapsed();
  CHECK(e == Approx(1.0));
}

TEST_CASE("test_loop_timer", "[utilities]")
{
  LoopTimer<FakeChronoClock> loop;
  double it_time = loop.get_time_per_iteration();
  CHECK(it_time == Approx(0.0));

  loop.start();
  loop.stop();
  it_time = loop.get_time_per_iteration();
  CHECK(it_time == Approx(1.0));

  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(2.0s);
  loop.start();
  loop.stop();
  it_time = loop.get_time_per_iteration();
  CHECK(it_time == Approx(1.5)); // 2 iterations
  // restore value
  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.0s);
}

TEST_CASE("test_loop_control", "[utilities]")
{
  // fake clock advances every time "now" is called
  LoopTimer<FakeChronoClock> loop;
  int max_cpu_secs = 9;
  RunTimeManager<FakeChronoClock> rm; // fake clock = 1
  RunTimeControl<FakeChronoClock> rc(rm, max_cpu_secs, "dummy", false);
  rc.runtime_padding(1.0);
  REQUIRE(!rc.checkStop(loop)); // fake clock = 2

  loop.start(); // fake clock = 3
  loop.stop();  // fake clock = 4
  // fake clock = 5
  // estimated time with margin and padding =  1.0 sec/it * 1.1 + 1.0 (pad) = 2.1
  // remaining = 9 - 4.0 = 5.0  enough time for another loop.
  REQUIRE(!rc.checkStop(loop));

  loop.start(); // fake clock = 6
  loop.stop();  // fake clock = 7
  // fake clock = 8
  // estimated time with margin and padding = 1.0 sec/it * 1.1 + 1.0 = 2.1
  // remaining = 9 - 8.0 = 1.0  not enough time for another loop
  REQUIRE(rc.checkStop(loop));

  std::string msg = rc.generateStopMessage("QMC", 2);
  REQUIRE(msg.size() > 0);
}


} // namespace qmcplusplus
