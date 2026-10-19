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
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Utilities/RunTimeManager.h"
#include <cstdio>
#include <fstream>
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
  CHECK(msg.find("RunTimeControl takes action in QMC driver.") != std::string::npos);
  CHECK(msg.find("Insufficient time for the current block. Stopping before block 2") != std::string::npos);
}

TEST_CASE("test_loop_control_average", "[utilities]")
{
  // Advance time explicitly so only the measured iterations contribute to the estimate.
  const auto saved_increment                   = FakeChronoClock::fake_chrono_clock_increment;
  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(0s);
  RunTimeManager<FakeChronoClock> rm;
  RunTimeControl<FakeChronoClock> rc(rm, 100, "dummy", false);
  LoopTimer<FakeChronoClock> loop;

  loop.start();
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(2s);
  loop.stop();
  loop.start();
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(8s);
  loop.stop();
  CHECK(loop.get_time_per_iteration() == Approx(5.0));

  // With 36 seconds left, the average fits (5 * 1.1 + 30), but the last iteration would not.
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(54s);
  CHECK_FALSE(rc.checkStop(loop));
  // Sufficient time remaining (35.6s left > 35.5s needed).
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(0.4s);
  CHECK_FALSE(rc.checkStop(loop));
  // Insufficient time for next iteration (35.4s left < 35.5s needed).
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(0.2s);
  CHECK(rc.checkStop(loop));
  const std::string msg = rc.generateStopMessage("executeLoop", 2, "loop iteration");
  CHECK(msg.find("RunTimeControl takes action in executeLoop.") != std::string::npos);
  CHECK(msg.find("Insufficient time for the current loop iteration. Stopping before loop iteration 2") != std::string::npos);
  CHECK(msg.find("loop iteration (seconds) = 5") != std::string::npos);

  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(40s);
  CHECK(rc.checkStop(loop));
  CHECK(rc.generateStopMessage("executeLoop", 2, "loop iteration")
            .find("Time limit reached. Stopping before loop iteration 2") != std::string::npos);
  FakeChronoClock::fake_chrono_clock_increment = saved_increment;
}

TEST_CASE("test_loop_control_before_first_iteration", "[utilities]")
{
  RunTimeManager<FakeChronoClock> rm;
  LoopTimer<FakeChronoClock> loop;
  SECTION("Time limit already reached")
  {
    RunTimeControl<FakeChronoClock> rc(rm, 0, "dummy", false);
    CHECK(rc.checkStop(loop));
    CHECK(rc.generateStopMessage("executeLoop", 0, "loop iteration")
              .find("Time limit reached. Stopping before loop iteration 0") != std::string::npos);
  }
  SECTION("Insufficient shutdown allowance")
  {
    RunTimeControl<FakeChronoClock> rc(rm, 30, "dummy", false);
    CHECK(rc.checkStop(loop));
    CHECK(rc.generateStopMessage("executeLoop", 0, "loop iteration")
              .find("Insufficient time for the current loop iteration. Stopping before loop iteration 0") != std::string::npos);
  }
  SECTION("Enough time to start")
  {
    RunTimeControl<FakeChronoClock> rc(rm, 100, "dummy", false);
    CHECK_FALSE(rc.checkStop(loop));
  }
}


TEST_CASE("test_loop_control_progress_message", "[utilities]")
{
  RunTimeManager<FakeChronoClock> rm;
  RunTimeControl<FakeChronoClock> rc(rm, 100, "dummy", false);
  LoopTimer<FakeChronoClock> loop;
  loop.start();
  loop.stop();

  CHECK(rc.generateProgressMessage(loop, "QMC", 0, 1).find("average 1 secs/block after ") != std::string::npos);
}

TEST_CASE("test_loop_control_stop_file", "[utilities]")
{
  const std::string prefix    = "runtime_stop_file_test";
  const std::string stop_file = prefix + ".STOP";
  RunTimeManager<FakeChronoClock> rm;
  RunTimeControl<FakeChronoClock> rc(rm, 1000, prefix);
  LoopTimer<FakeChronoClock> loop;
  {
    std::ofstream stop(stop_file);
    CHECK(stop.good());
  }
  CHECK(rc.checkStop(loop));
  CHECK(rc.generateStopMessage("QMC", 0).find("Stop requested from the control file") != std::string::npos);
  CHECK(std::remove(stop_file.c_str()) == 0);
}

} // namespace qmcplusplus
