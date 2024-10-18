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

#include <string>
#include <vector>
#include "Utilities/TimerManager.h"

using namespace std::chrono_literals;

namespace qmcplusplus
{
using FakeTimerManager = TimerManager<FakeTimer>;

template<class CLOCK>
void set_total_time(TimerType<CLOCK>* timer, double total_time_input)
{
  timer->total_time = total_time_input;
}

template<class CLOCK>
void set_num_calls(TimerType<CLOCK>* timer, long num_calls_input)
{
  timer->num_calls = num_calls_input;
}

// Convert duration input type to nanosecond duration
template<typename T>
FakeChronoClock::duration convert_to_ns(T in)
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(in);
}

// Convert duration input type to seconds as double precision type
template<typename T>
double convert_to_s(T in)
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(in).count();
}

TEST_CASE("test_timer_stack", "[utilities]")
{
  // Use a local version rather than the global timer_manager, otherwise
  //  changes will persist from test to test.
  FakeTimerManager tm;
  FakeTimer* t1 = tm.createTimer("timer1", timer_level_coarse);
#if defined(ENABLE_TIMERS)
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
  FakeTimerManager tm;
  FakeTimer* t1 = tm.createTimer("timer1", timer_level_coarse);
  {
    ScopedFakeTimer st(*t1);
  }
#if defined(ENABLE_TIMERS)
  CHECK(t1->get_total() == Approx(1.0));
  REQUIRE(t1->get_num_calls() == 1);
#endif


  const std::string prefix_str("Prefix::");
  FakeTimer& t2(*(tm.createTimer(prefix_str + "timer2", timer_level_coarse)));
  {
    ScopedFakeTimer st(t2);
  }
#if defined(ENABLE_TIMERS)
  CHECK(t1->get_total() == Approx(1.0));
  REQUIRE(t1->get_num_calls() == 1);
#endif
}

TEST_CASE("test_timer_flat_profile", "[utilities]")
{
  FakeTimerManager tm;
  FakeTimer* t1 = tm.createTimer("timer1");
  set_total_time(t1, 1.1);
  set_num_calls(t1, 2);

  FakeTimerManager::FlatProfileData p;
  tm.collate_flat_profile(NULL, p);

  REQUIRE(p.nameList.size() == 1);
  REQUIRE(p.nameList.at("timer1") == 0);
  REQUIRE(p.timeList.size() == 1);
  CHECK(p.timeList[0] == Approx(1.1));
  REQUIRE(p.callList.size() == 1);
  REQUIRE(p.callList[0] == 2);
}

TEST_CASE("test_timer_flat_profile_same_name", "[utilities]")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  FakeTimer* t1 = tm.createTimer("timer1");
  FakeTimer* t2 = tm.createTimer("timer2");
  FakeTimer* t3 = tm.createTimer("timer1");

  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.1s);
  t1->start();
  t1->stop();
  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.2s);
  for (int i = 0; i < 3; i++)
  {
    t2->start();
    t2->stop();

    t3->start();
    t3->stop();
  }
  t3->start();
  t3->stop();

  FakeTimerManager::FlatProfileData p;

  tm.collate_flat_profile(NULL, p);

#ifdef ENABLE_TIMERS
  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  CHECK(p.timeList[idx1] == Approx(5.9));
  CHECK(p.timeList[idx2] == Approx(3.6));

  REQUIRE(p.callList.size() == 2);
  REQUIRE(p.callList[idx1] == 5);
  REQUIRE(p.callList[idx2] == 3);
#endif
}

TEST_CASE("test_timer_nested_profile", "[utilities]")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  FakeTimer* t1 = tm.createTimer("timer1");
  FakeTimer* t2 = tm.createTimer("timer2");

  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.1s);
  t1->start();
  t2->start();
  t2->stop();
  t1->stop();

  FakeTimerManager::FlatProfileData p;
  tm.collate_flat_profile(NULL, p);

#ifdef ENABLE_TIMERS
  REQUIRE(p.nameList.size() == 2);
  int idx1 = p.nameList.at("timer1");
  int idx2 = p.nameList.at("timer2");
  REQUIRE(p.timeList.size() == 2);
  CHECK(p.timeList[idx1] == Approx(3 * convert_to_s(FakeChronoClock::fake_chrono_clock_increment)));
  CHECK(p.timeList[idx2] == Approx(convert_to_s(FakeChronoClock::fake_chrono_clock_increment)));
#endif

  FakeTimerManager::StackProfileData p2;
  tm.collate_stack_profile(NULL, p2);

#ifdef ENABLE_TIMERS
  REQUIRE(p2.nameList.size() == 2);
  idx1 = p2.nameList.at("timer1");
  idx2 = p2.nameList.at("timer1/timer2");
  REQUIRE(p2.timeList.size() == 2);
  REQUIRE(p2.timeExclList.size() == 2);
  CHECK(p2.timeList[idx1] == Approx(convert_to_s(3 * FakeChronoClock::fake_chrono_clock_increment)));
  CHECK(p2.timeList[idx2] == Approx(convert_to_s(FakeChronoClock::fake_chrono_clock_increment)));

  // Time in t1 minus time inside t2
  CHECK(p2.timeExclList[idx1] == Approx(2 * convert_to_s(FakeChronoClock::fake_chrono_clock_increment)));
  CHECK(p2.timeExclList[idx2] == Approx(convert_to_s(FakeChronoClock::fake_chrono_clock_increment)));
#endif
}

TEST_CASE("test_timer_nested_profile_two_children", "[utilities]")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  FakeTimer* t1 = tm.createTimer("timer1");
  FakeTimer* t2 = tm.createTimer("timer2");
  FakeTimer* t3 = tm.createTimer("timer3");

  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.1s);
  t1->start();
  t2->start();
  t2->stop();
  t3->start();
  t3->stop();
  t1->stop();

  FakeTimerManager::StackProfileData p2;
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
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  FakeTimer* t1 = tm.createTimer("timer1");
  FakeTimer* t2 = tm.createTimer("timer2");
  FakeTimer* t3 = tm.createTimer("timer3");
  FakeTimer* t4 = tm.createTimer("timer4");
  FakeTimer* t5 = tm.createTimer("timer5");


  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.1s);
  t1->start();
  t2->start();
  t3->start();
  t4->start();
  t4->stop();
  t5->start();
  t5->stop();
  t3->stop();
  t2->stop();
  t3->start();
  t4->start();
  t4->stop();
  t3->stop();
  t1->stop();

  FakeTimerManager::StackProfileData p2;
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

  CHECK(p2.timeList[0] == Approx(14.3));
  CHECK(p2.timeExclList[0] == Approx(3.3));
  CHECK(p2.timeList[1] == Approx(7.7));
  CHECK(p2.timeExclList[1] == Approx(2.2));
#endif

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp2.xml");
}

TEST_CASE("test_timer_nested_profile_collate", "[utilities]")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  FakeTimer* t1  = tm.createTimer("timer1");
  FakeTimer* t2  = tm.createTimer("timer2");
  FakeTimer* t2b = tm.createTimer("timer2");
  FakeTimer* t3  = tm.createTimer("timer3");


  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.1s);
  t1->start();
  t2->start();
  t3->start();
  t3->stop();
  t2->stop();
  t2b->start();
  t3->start();
  t3->stop();
  t2b->stop();
  t2->start();
  t3->start();
  t3->stop();
  t2->stop();
  t2b->start();
  t3->start();
  t3->stop();
  t2b->stop();
  t1->stop();


  FakeTimerManager::StackProfileData p2;
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
  for (int i = 0; i < StackKey::max_level + 1; i++)
  {
    sk.add_id(1);
  }
  REQUIRE(timer_max_level_exceeded == true);
}

TEST_CASE("test stack exceeded message")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  std::vector<FakeTimer*> timer_list;
  for (int i = 0; i < StackKey::max_level + 1; i++)
  {
    std::ostringstream name;
    name << "timer" << i;
    FakeTimer* t = tm.createTimer(name.str());
    timer_list.push_back(t);
  }

  for (int i = 0; i < StackKey::max_level + 1; i++)
    timer_list[i]->start();

  for (int i = StackKey::max_level; i >= 0; i--)
    timer_list[i]->stop();

  REQUIRE(timer_max_level_exceeded == true);

  //tm.print_stack(NULL);

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp4.xml");
}

TEST_CASE("test max exceeded message")
{
  FakeTimerManager tm;
  tm.set_timer_threshold(timer_level_fine);
  std::vector<FakeTimer*> timer_list;
  for (int i = 0; i < std::numeric_limits<timer_id_t>::max() + 1; i++)
  {
    std::ostringstream name;
    name << "timer" << i;
    FakeTimer* t = tm.createTimer(name.str());
    timer_list.push_back(t);
  }

  Libxml2Document doc;
  doc.newDoc("resources");
  tm.output_timing(NULL, doc, doc.getRoot());
  doc.dump("tmp5.xml");
}
#endif

// Define a list of timers indexed by an enum
// First, define an enum with the timers
enum TestTimer
{
  MyTimer1,
  MyTimer2,
};

// Next define a structure mapping the enum to a string name
TimerNameList_t<TestTimer> TestTimerNames = {{MyTimer1, "Timer name 1"}, {MyTimer2, "Timer name 2"}};

TEST_CASE("test setup timers", "[utilities]")
{
  FakeTimerManager tm;
  // Create  a list of timers and initialize it
  TimerList Timers(tm, TestTimerNames, timer_level_coarse);

  FakeChronoClock::fake_chrono_clock_increment = convert_to_ns(1.0s);
  Timers[MyTimer1].get().start();
  Timers[MyTimer1].get().stop();

#ifdef ENABLE_TIMERS
  CHECK(Timers[MyTimer1].get().get_total() == Approx(1.0));
  REQUIRE(Timers[MyTimer1].get().get_num_calls() == 1);
#endif
}

} // namespace qmcplusplus
