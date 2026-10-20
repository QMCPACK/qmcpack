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
#include "Utilities/ProjectData.h"
#include "OhmmsData/Libxml2Doc.h"
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

namespace
{
ProjectData make_project_data(const std::string& title, int max_seconds)
{
  Libxml2Document doc;
  doc.parseFromString("<project id=\"" + title + "\"><parameter name=\"max_seconds\">" + std::to_string(max_seconds) +
                      "</parameter></project>");
  ProjectData proj;
  proj.put(doc.getRoot());
  return proj;
}
} // namespace

TEST_CASE("test_loop_control", "[utilities]")
{
  // fake clock advances every time "now" is called
  LoopTimer<FakeChronoClock> loop;
  int max_cpu_secs = 9;
  RunTimeManager<FakeChronoClock> rm; // fake clock = 1
  ProjectData proj = make_project_data("dummy", max_cpu_secs);
  LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "QMC driver", "block", false);
  rc.runtime_padding(1.0);
  REQUIRE_FALSE(rc.checkStop(loop, 0)); // fake clock = 2

  loop.start(); // fake clock = 3
  loop.stop();  // fake clock = 4
  // fake clock = 5
  // estimated time with margin and padding =  1.0 sec/it * 1.1 + 1.0 (pad) = 2.1
  // remaining = 9 - 4.0 = 5.0  enough time for another loop.
  REQUIRE_FALSE(rc.checkStop(loop, 1));

  loop.start(); // fake clock = 6
  loop.stop();  // fake clock = 7
  // fake clock = 8
  // estimated time with margin and padding = 1.0 sec/it * 1.1 + 1.0 = 2.1
  // remaining = 9 - 8.0 = 1.0  not enough time for another loop
  REQUIRE(rc.checkStop(loop, 2));

  const std::string& msg = rc.getStopMessage();
  CHECK(msg.find("LoopRunTimeControl takes action in QMC driver.") != std::string::npos);
  CHECK(msg.find("Insufficient time for the current block. Stopping before block 2") != std::string::npos);
}

TEST_CASE("test_loop_control_average", "[utilities]")
{
  // Advance time explicitly so only the measured iterations contribute to the estimate.
  const auto saved_increment                   = FakeChronoClock::fake_chrono_clock_increment;
  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(0s);
  RunTimeManager<FakeChronoClock> rm;
  ProjectData proj = make_project_data("dummy", 100);
  LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "executeLoop", "loop iteration", false);
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
  CHECK_FALSE(rc.checkStop(loop, 2));
  // Sufficient time remaining (35.6s left > 35.5s needed).
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(0.4s);
  CHECK_FALSE(rc.checkStop(loop, 2));
  // Insufficient time for next iteration (35.4s left < 35.5s needed).
  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(0.2s);
  CHECK(rc.checkStop(loop, 2));
  const std::string msg = rc.getStopMessage();
  CHECK(msg.find("LoopRunTimeControl takes action in executeLoop.") != std::string::npos);
  CHECK(msg.find("Insufficient time for the current loop iteration. Stopping before loop iteration 2") != std::string::npos);
  CHECK(msg.find("loop iteration (seconds) = 5") != std::string::npos);

  FakeChronoClock::fake_chrono_clock_value += convert_to_ns(40s);
  CHECK(rc.checkStop(loop, 2));
  CHECK(rc.getStopMessage().find("Time limit reached. Stopping before loop iteration 2") != std::string::npos);
  FakeChronoClock::fake_chrono_clock_increment = saved_increment;
}

TEST_CASE("test_loop_control_before_first_iteration", "[utilities]")
{
  RunTimeManager<FakeChronoClock> rm;
  LoopTimer<FakeChronoClock> loop;
  SECTION("Time limit already reached")
  {
    ProjectData proj = make_project_data("dummy", 0);
    LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "executeLoop", "loop iteration", false);
    CHECK(rc.checkStop(loop, 0));
    CHECK(rc.getStopMessage().find("Time limit reached. Stopping before loop iteration 0") != std::string::npos);
  }
  SECTION("Insufficient shutdown allowance")
  {
    ProjectData proj = make_project_data("dummy", 30);
    LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "executeLoop", "loop iteration", false);
    CHECK(rc.checkStop(loop, 0));
    CHECK(rc.getStopMessage().find("Insufficient time for the current loop iteration. Stopping before loop iteration 0") != std::string::npos);
  }
  SECTION("Enough time to start")
  {
    ProjectData proj = make_project_data("dummy", 100);
    LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "executeLoop", "loop iteration", false);
    CHECK_FALSE(rc.checkStop(loop, 0));
  }
}


TEST_CASE("test_loop_control_progress_message", "[utilities]")
{
  RunTimeManager<FakeChronoClock> rm;
  ProjectData proj = make_project_data("dummy", 100);
  LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "QMC driver", "block", false);
  LoopTimer<FakeChronoClock> loop;
  loop.start();
  loop.stop();

  CHECK(rc.generateProgressMessage(loop, 0, 1).find("average 1 secs/block after ") != std::string::npos);
}

TEST_CASE("test_loop_control_stop_file", "[utilities]")
{
  const std::string prefix    = "runtime_stop_file_test";
  const std::string stop_file = prefix + ".STOP";
  RunTimeManager<FakeChronoClock> rm;
  ProjectData proj = make_project_data(prefix, 1000);
  LoopRunTimeControl<FakeChronoClock> rc(rm, proj, "QMC driver", "block");
  LoopTimer<FakeChronoClock> loop;
  {
    std::ofstream stop(stop_file);
    CHECK(stop.good());
  }
  CHECK(rc.checkStop(loop, 0));
  CHECK(rc.getStopMessage().find("Stop requested from the control file") != std::string::npos);
  CHECK(std::remove(stop_file.c_str()) == 0);
}

TEST_CASE("test_runtime_control_checkStop", "[utilities]")
{
  RunTimeManager<FakeChronoClock> rm;
  ProjectData proj = make_project_data("test_proj", 10);
  RunTimeControl<FakeChronoClock> rc(rm, proj, "test_scope");

  // Initially: elapsed is 0, max_seconds is 10, no STOP file
  CHECK_FALSE(rc.checkStop());
  CHECK(rc.getStopMessage().empty());
  CHECK_FALSE(rm.isStopNeeded());

  SECTION("Time limit reached")
  {
    FakeChronoClock::fake_chrono_clock_value += convert_to_ns(11s);
    CHECK(rc.checkStop());
    CHECK_FALSE(rc.getStopMessage().empty());
    CHECK(rc.getStopMessage().find("RunTimeControl takes action in test_scope.") != std::string::npos);
    CHECK(rc.getStopMessage().find("Time limit reached.") != std::string::npos);
    CHECK(rc.getStopMessage().find("Hard limit (seconds) 10") != std::string::npos);
    // checkStop does not touch need_to_stop_
    CHECK_FALSE(rm.isStopNeeded());

    rm.markStop();
    CHECK(rm.isStopNeeded());
  }

  SECTION("Stop file requested")
  {
    const std::string stop_file = proj.getTitle() + ".STOP";
    {
      std::ofstream stop(stop_file);
      CHECK(stop.good());
    }
    CHECK(rc.checkStop());
    CHECK_FALSE(rc.getStopMessage().empty());
    CHECK(rc.getStopMessage().find("RunTimeControl takes action in test_scope.") != std::string::npos);
    CHECK(rc.getStopMessage().find("Stop requested from the control file \"" + stop_file + "\"") != std::string::npos);
    // checkStop does not touch need_to_stop_
    CHECK_FALSE(rm.isStopNeeded());

    rm.markStop();
    CHECK(rm.isStopNeeded());
    CHECK(std::remove(stop_file.c_str()) == 0);
  }

  SECTION("markStop flags stop needed")
  {
    RunTimeManager<FakeChronoClock> rm2;
    CHECK_FALSE(rm2.isStopNeeded());
    rm2.markStop();
    CHECK(rm2.isStopNeeded());
  }

  SECTION("subsequent checks become no-op when isStopNeeded is true")
  {
    rm.markStop();
    CHECK(rm.isStopNeeded());

    // Even if time limit is exceeded, RunTimeControl::checkStop() is a no-op returning false
    FakeChronoClock::fake_chrono_clock_value += convert_to_ns(100s);
    CHECK_FALSE(rc.checkStop());
    CHECK(rc.getStopMessage().empty());

    // LoopRunTimeControl::checkStop() is also a no-op returning false
    LoopTimer<FakeChronoClock> loop;
    LoopRunTimeControl<FakeChronoClock> loop_rc(rm, proj, "test_loop", "block", false);
    CHECK_FALSE(loop_rc.checkStop(loop, 0));
    CHECK(loop_rc.getStopMessage().empty());
  }
}

} // namespace qmcplusplus
