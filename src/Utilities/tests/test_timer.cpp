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

#define USE_FAKE_CLOCK
#include "Utilities/NewTimer.h"
#include <stdio.h>
#include <string>
#include <vector>

namespace qmcplusplus {


// Used by fake_cpu_clock in Clock.h if USE_FAKE_CLOCK is defined
double fake_cpu_clock_increment = 1.0;
double fake_cpu_clock_value = 0.0;

class FakeTimer : public NewTimer
{
public:
  FakeTimer(const std::string& myname) : NewTimer(myname) {}

  void set_total_time(double my_total_time)
  {
    total_time = my_total_time;
  }

  void set_num_calls(long my_num_calls)
  {
    num_calls = my_num_calls;
  }

};

TEST_CASE("test_timer_stack", "[utilities]")
{
  // Use a local version rather than the global TimerManager, otherwise
  //  changes will persist from test to test.
  TimerManagerClass tm;
  NewTimer t1("timer1");
  tm.addTimer(&t1);
#if ENABLE_TIMER
#ifdef USE_STACK_TIMERS
  t1.start();
  REQUIRE(tm.current_timer() == &t1);
  t1.stop();
  REQUIRE(tm.current_timer() == NULL);
#endif
#endif
}

TEST_CASE("test_timer_flat_profile", "[utilities]")
{
  TimerManagerClass tm;
  FakeTimer t1("timer1");
  tm.addTimer(&t1);
  t1.set_total_time(1.1);
  t1.set_num_calls(2);

  TimerManagerClass::FlatProfileData p;
  tm.collate_flat_profile(NULL, p);

  REQUIRE(p.nameList.size() == 1);
  REQUIRE(p.nameList.at("timer1") == 0);
  REQUIRE(p.timeList.size() == 1);
  REQUIRE(p.timeList[0] == Approx(1.1));
  REQUIRE(p.callList.size() == 1);
  REQUIRE(p.callList[0] == 2);
}

TEST_CASE("test_timer_flat_profile_same_name", "[utilities]")
{
  TimerManagerClass tm;
  FakeTimer t1("timer1");
  tm.addTimer(&t1);
  FakeTimer t2("timer2");
  tm.addTimer(&t2);
  FakeTimer t3("timer1");
  tm.addTimer(&t3);

  fake_cpu_clock_increment = 1.1;
  t1.start();
  t1.stop();
  fake_cpu_clock_increment = 1.2;
  for (int i = 0; i < 3; i++)
  {
    t2.start();
    t2.stop();

    t3.start();
    t3.stop();
  }
  t3.start();
  t3.stop();

  TimerManagerClass::FlatProfileData p;

  tm.collate_flat_profile(NULL, p);

  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  REQUIRE(p.timeList[idx1] == Approx(5.9));
  REQUIRE(p.timeList[idx2] == Approx(3.6));

  REQUIRE(p.callList.size() == 2);
  REQUIRE(p.callList[idx1] == 5);
  REQUIRE(p.callList[idx2] == 3);
}

TEST_CASE("test_timer_nested_profile", "[utilities]")
{
  TimerManagerClass tm;
  FakeTimer t1("timer1");
  tm.addTimer(&t1);
  FakeTimer t2("timer2");
  tm.addTimer(&t2);

  fake_cpu_clock_increment = 1.1;
  t1.start();
  t2.start();
  t2.stop();
  t1.stop();

  TimerManagerClass::FlatProfileData p;
  tm.collate_flat_profile(NULL, p);

  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  REQUIRE(p.timeList[idx1] == Approx(3*fake_cpu_clock_increment));
  REQUIRE(p.timeList[idx2] == Approx(fake_cpu_clock_increment));

  TimerManagerClass::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);

  REQUIRE(p2.nameList.size() == 2);
  idx1 = p2.nameList.at("timer1");
  idx2 = p2.nameList.at("timer2/timer1");
  REQUIRE(p2.timeList.size() == 2);
  REQUIRE(p2.timeExclList.size() == 2);
  REQUIRE(p2.timeList[idx1] == Approx(3*fake_cpu_clock_increment));
  REQUIRE(p2.timeList[idx2] == Approx(fake_cpu_clock_increment));

  // Time in t1 minus time inside t2
  REQUIRE(p2.timeExclList[idx1] == Approx(2*fake_cpu_clock_increment));
  REQUIRE(p2.timeExclList[idx2] == Approx(fake_cpu_clock_increment));

  std::vector<NewTimer *> roots;
  tm.get_stack_roots(roots);
  REQUIRE(roots.size() == 1);
  REQUIRE(roots[0] == &t1);

  TimerDFS dfs(&t1);

  NewTimer *t = dfs.next();
  REQUIRE(t == &t2);
  t = dfs.next();
  REQUIRE(t == NULL);

}

TEST_CASE("test_timer_nested_profile_two_children", "[utilities]")
{
  TimerManagerClass tm;
  NewTimer t1("timer1");
  tm.addTimer(&t1);
  NewTimer t2("timer2");
  tm.addTimer(&t2);
  NewTimer t3("timer3");
  tm.addTimer(&t3);

  fake_cpu_clock_increment = 1.1;
  t1.start();
  t2.start();
  t2.stop();
  t3.start();
  t3.stop();
  t1.stop();

  TimerManagerClass::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);

  std::vector<NewTimer *> roots;
  tm.get_stack_roots(roots);
  REQUIRE(roots.size() == 1);
  REQUIRE(roots[0] == &t1);

  TimerDFS dfs(&t1);

  NewTimer *t = dfs.next();
  REQUIRE(t == &t2);
  t = dfs.next();
  REQUIRE(t == &t3);
  t = dfs.next();
  REQUIRE(t == NULL);

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp.xml");
  // To really test this, should read the file in and inspect the contents.
  // For now, it makes for quick iterations on writing the file.

}

}
