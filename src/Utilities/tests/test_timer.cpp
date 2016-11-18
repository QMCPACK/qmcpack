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
  t1.start();
  REQUIRE(tm.current_timer() == &t1);
  t1.stop();
  REQUIRE(tm.current_timer() == NULL);
#endif
}

TEST_CASE("test_timer_flat_profile", "[utilities]")
{
  TimerManagerClass tm;
  FakeTimer t1("timer1");
  tm.addTimer(&t1);
  t1.set_total_time(1.1);
  t1.set_num_calls(2);

  TimerManagerClass::nameList_t nameList; 
  TimerManagerClass::timeList_t timeList;
  TimerManagerClass::callList_t callList;
  tm.compute_flat_profile(NULL, nameList, timeList, callList);

  REQUIRE(nameList.size() == 1);
  REQUIRE(nameList.at("timer1") == 0);
  REQUIRE(timeList.size() == 1);
  REQUIRE(timeList[0] == Approx(1.1));
  REQUIRE(callList.size() == 1);
  REQUIRE(callList[0] == 2);
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

  t1.set_total_time(1.1);
  t1.set_num_calls(1);
  t2.set_total_time(4.0);
  t2.set_num_calls(3);
  t3.set_total_time(2.1);
  t3.set_num_calls(4);

  TimerManagerClass::nameList_t nameList; 
  TimerManagerClass::timeList_t timeList;
  TimerManagerClass::callList_t callList;
  tm.compute_flat_profile(NULL, nameList, timeList, callList);

  REQUIRE(nameList.size() == 2);
  int idx1 = nameList.at("timer1");
  int idx2 = nameList.at("timer2");
  REQUIRE(timeList.size() == 2);
  REQUIRE(timeList[idx1] == Approx(3.2));
  REQUIRE(timeList[idx2] == Approx(4.0));
  
  REQUIRE(callList.size() == 2);
  REQUIRE(callList[idx1] == 5);
  REQUIRE(callList[idx2] == 3);
}

}
