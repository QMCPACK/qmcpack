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
  NewTimer *t1 = tm.createTimer("timer1", timer_level_coarse);
#if ENABLE_TIMERS
#ifdef USE_STACK_TIMERS
  t1->start();
  REQUIRE(tm.current_timer() == t1);
  t1->stop();
  REQUIRE(tm.current_timer() == NULL);
#endif
#endif
}

TEST_CASE("test_timer_scoped", "[utilities]")
{
  TimerManagerClass tm;
  NewTimer *t1 = tm.createTimer("timer1", timer_level_coarse);
  {
    ScopedTimer st(t1);
  }
#if ENABLE_TIMERS
  REQUIRE(t1->get_total() == Approx(1.0));
  REQUIRE(t1->get_num_calls() == 1);
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
  tm.set_timer_threshold(timer_level_fine);
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

#ifdef ENABLE_TIMERS
  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  REQUIRE(p.timeList[idx1] == Approx(5.9));
  REQUIRE(p.timeList[idx2] == Approx(3.6));

  REQUIRE(p.callList.size() == 2);
  REQUIRE(p.callList[idx1] == 5);
  REQUIRE(p.callList[idx2] == 3);
#endif
}

TEST_CASE("test_timer_nested_profile", "[utilities]")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
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

#ifdef ENABLE_TIMERS
  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  REQUIRE(p.timeList[idx1] == Approx(3*fake_cpu_clock_increment));
  REQUIRE(p.timeList[idx2] == Approx(fake_cpu_clock_increment));
#endif

  TimerManagerClass::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);

#ifdef ENABLE_TIMERS
  REQUIRE(p2.nameList.size() == 2);
  idx1 = p2.nameList.at("timer1");
  idx2 = p2.nameList.at("timer1/timer2");
  REQUIRE(p2.timeList.size() == 2);
  REQUIRE(p2.timeExclList.size() == 2);
  REQUIRE(p2.timeList[idx1] == Approx(3*fake_cpu_clock_increment));
  REQUIRE(p2.timeList[idx2] == Approx(fake_cpu_clock_increment));

  // Time in t1 minus time inside t2
  REQUIRE(p2.timeExclList[idx1] == Approx(2*fake_cpu_clock_increment));
  REQUIRE(p2.timeExclList[idx2] == Approx(fake_cpu_clock_increment));
#endif

}

TEST_CASE("test_timer_nested_profile_two_children", "[utilities]")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
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

#ifdef ENABLE_TIMERS
  REQUIRE(p2.names.size() == 3);
  REQUIRE(p2.names[0] == "timer1");
  REQUIRE(p2.names[1] == "timer1/timer2");
  REQUIRE(p2.names[2] == "timer1/timer3");
#endif


  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp.xml");
  // To really test this, should read the file in and inspect the contents.
  // For now, it makes for quick iterations on writing the file.

}

TEST_CASE("test_timer_nested_profile_alt_routes", "[utilities]")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
  NewTimer t1("timer1");
  tm.addTimer(&t1);
  NewTimer t2("timer2");
  tm.addTimer(&t2);
  NewTimer t3("timer3");
  tm.addTimer(&t3);
  NewTimer t4("timer4");
  tm.addTimer(&t4);
  NewTimer t5("timer5");
  tm.addTimer(&t5);


  fake_cpu_clock_increment = 1.1;
  t1.start();
    t2.start();
      t3.start();
        t4.start();
        t4.stop();
        t5.start();
        t5.stop();
      t3.stop();
    t2.stop();
    t3.start();
      t4.start();
      t4.stop();
    t3.stop();
  t1.stop();

  TimerManagerClass::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);
  //tm.print_stack(NULL);
#ifdef ENABLE_TIMERS
  REQUIRE(p2.names.size() == 7);
  REQUIRE(p2.names[0] == "timer1");
  REQUIRE(p2.names[1] == "timer1/timer2");
  REQUIRE(p2.names[2] == "timer1/timer2/timer3");
  REQUIRE(p2.names[3] == "timer1/timer2/timer3/timer4");
  REQUIRE(p2.names[4] == "timer1/timer2/timer3/timer5");
  REQUIRE(p2.names[5] == "timer1/timer3");
  REQUIRE(p2.names[6] == "timer1/timer3/timer4");

  REQUIRE(p2.timeList[0] == Approx(14.3));
  REQUIRE(p2.timeExclList[0] == Approx(3.3));
  REQUIRE(p2.timeList[1] == Approx(7.7));
  REQUIRE(p2.timeExclList[1] == Approx(2.2));
#endif

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp2.xml");
}

TEST_CASE("test_timer_nested_profile_collate", "[utilities]")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
  NewTimer t1("timer1");
  tm.addTimer(&t1);
  NewTimer t2("timer2");
  tm.addTimer(&t2);
  NewTimer t2b("timer2");
  tm.addTimer(&t2b);
  NewTimer t3("timer3");
  tm.addTimer(&t3);


  fake_cpu_clock_increment = 1.1;
  t1.start();
    t2.start();
      t3.start();
      t3.stop();
    t2.stop();
    t2b.start();
      t3.start();
      t3.stop();
    t2b.stop();
    t2.start();
      t3.start();
      t3.stop();
    t2.stop();
    t2b.start();
      t3.start();
      t3.stop();
    t2b.stop();
  t1.stop();


  TimerManagerClass::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);
  //tm.print_stack(NULL);
#ifdef ENABLE_TIMERS
  REQUIRE(p2.names.size() == 3);
  REQUIRE(p2.names[0] == "timer1");
  REQUIRE(p2.names[1] == "timer1/timer2");
  REQUIRE(p2.names[2] == "timer1/timer2/timer3");
#endif

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp3.xml");
}

#ifdef ENABLE_TIMERS
TEST_CASE("test stack key")
{
  StackKey sk;
  REQUIRE(timer_max_level_exceeded == false);
  for (int i = 0; i < StackKey::max_level+1; i++)
  {
    sk.add_id(1);
  }
  REQUIRE(timer_max_level_exceeded == true);
}

TEST_CASE("test stack exceeded message")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
  std::vector<NewTimer *> timer_list;
  for (int i = 0; i < StackKey::max_level+1; i++)
  {
    std::ostringstream name;
    name << "timer" << i;
    NewTimer *t = new NewTimer(name.str());
    tm.addTimer(t);
    timer_list.push_back(t);
  }
  for (int i = 0; i < StackKey::max_level+1; i++)
  {
    timer_list[i]->start();
  }
  for (int i = 0; i < StackKey::max_level+1; i++)
  {
    timer_list[i]->stop();
  }
  REQUIRE(timer_max_level_exceeded == true);

  //tm.print_stack(NULL);

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp4.xml");
}

TEST_CASE("test max exceeded message")
{
  TimerManagerClass tm;
  tm.set_timer_threshold(timer_level_fine);
  std::vector<NewTimer *> timer_list;
  for (int i = 0; i < std::numeric_limits<timer_id_t>::max()+1; i++)
  {
    std::ostringstream name;
    name << "timer" << i;
    NewTimer *t = new NewTimer(name.str());
    tm.addTimer(t);
    timer_list.push_back(t);
  }

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp5.xml");
}
#endif

#if __cplusplus >=201103l
// Define a list of timers indexed by an enum
// First, define an enum with the timers
enum TestTimer
{
  MyTimer1,
  MyTimer2,
};

// Next define a structure mapping the enum to a string name
TimerNameList_t<TestTimer> TestTimerNames =
{
  {MyTimer1, "Timer name 1"},
  {MyTimer2, "Timer name 2"}
};

TEST_CASE("test setup timers","[utilities]")
{
  TimerManagerClass tm;
  // Create  a list of timers and initialize it
  TimerList_t Timers;
  setup_timers(Timers, TestTimerNames, timer_level_coarse);

  fake_cpu_clock_increment = 1.0;
  Timers[MyTimer1]->start();
  Timers[MyTimer1]->stop();

#ifdef ENABLE_TIMERS
  REQUIRE(Timers[MyTimer1]->get_total() == Approx(1.0));
  REQUIRE(Timers[MyTimer1]->get_num_calls() == 1);
#endif
}
#endif

}
