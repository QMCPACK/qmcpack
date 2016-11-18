//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file NewTimer.h
 * @brief NewTimer class various high-resolution timers.
 */
#ifndef QMCPLUSPLUS_NEW_TIMER_H
#define QMCPLUSPLUS_NEW_TIMER_H

#include <Utilities/Clock.h>
#include <vector>
#include <string>
#include <algorithm>

class Communicate;

namespace qmcplusplus
{

class NewTimer;

class TimerManagerClass
{
protected:
  std::vector<NewTimer*> TimerList;
  std::vector<NewTimer*> CurrentTimerStack;
public:
  inline void addTimer (NewTimer* t)
  {
    #pragma omp critical
    {
      TimerList.push_back(t);
    }
  }

  void push_timer(NewTimer *t)
  {
    CurrentTimerStack.push_back(t);
  }

  void pop_timer()
  {
    CurrentTimerStack.pop_back();
  }

  NewTimer *current_timer()
  {
    if (CurrentTimerStack.size() > 0)
    {
      return CurrentTimerStack.back();
    }
    return NULL;
  }


  void reset();
  void print (Communicate* comm);
};

extern TimerManagerClass TimerManager;

/* Timer using omp_get_wtime  */
class NewTimer
{
protected:
  double start_time;
  double total_time;
  long num_calls;
  std::string name;
public:
#if not(ENABLE_TIMER)
  inline void start() {}
  inline void stop() {}
#else
  inline void start()
  {
    TimerManager.push_timer(this);
    start_time = cpu_clock();
  }

  inline void stop()
  {
    total_time += cpu_clock() - start_time;
    num_calls++;
    TimerManager.pop_timer();
  }
#endif

  inline double    get_total() const
  {
    return total_time;
  }

  inline long  get_num_calls() const
  {
    return num_calls;
  }

  inline std::string get_name() const
  {
    return name;
  }

  inline void reset()
  {
    num_calls = 0;
    total_time=0.0;
  }

  NewTimer(const std::string& myname) :
    total_time(0.0), num_calls(0), name(myname)
  { }

  void set_name(const std::string& myname)
  {
    name=myname;
  }
};

struct TimerComparator
{
  inline bool operator()(const NewTimer *a, const NewTimer *b)
  {
    return a->get_name() < b->get_name();
  }
};




}

#endif
