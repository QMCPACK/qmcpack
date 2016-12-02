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
#include <OhmmsData/Libxml2Doc.h>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <iostream>

#define USE_STACK_TIMERS

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
  void addTimer (NewTimer* t);

  void push_timer(NewTimer *t)
  {
    #pragma omp master
    {
      CurrentTimerStack.push_back(t);
    }
  }

  void pop_timer()
  {
    #pragma omp master
    {
      CurrentTimerStack.pop_back();
    }
  }

  NewTimer *current_timer()
  {
    NewTimer *current = NULL;
    # pragma omp critical
    if (CurrentTimerStack.size() > 0)
    {
       current = CurrentTimerStack.back();
    }
    return current;
  }


  void reset();
  void print (Communicate* comm);
  void print_flat (Communicate* comm);
  void print_stack (Communicate* comm);

  typedef std::map<std::string, int> nameList_t;
  typedef std::vector<double> timeList_t;
  typedef std::vector<long> callList_t;
  typedef std::vector<std::string> names_t;

  struct FlatProfileData {
    nameList_t nameList;
    timeList_t timeList;
    callList_t callList;
  };

  struct StackProfileData {
    names_t    names;
    nameList_t nameList;
    timeList_t timeList;
    timeList_t timeExclList;
    callList_t callList;
  };

  void collate_flat_profile(Communicate *comm, FlatProfileData &p);

  void collate_stack_profile(Communicate *comm, StackProfileData &p);

  void output_timing(Communicate *comm, Libxml2Document &doc, xmlNodePtr root);

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
#ifdef USE_STACK_TIMERS
  TimerManagerClass *manager;
  NewTimer *parent;

  std::vector<NewTimer *> children;
  std::map<std::string, double> per_stack_total_time;
  std::map<std::string, long> per_stack_num_calls;
#endif
public:
#if not(ENABLE_TIMER)
  inline void start() {}
  inline void stop() {}
#else
  void start()
  {
#ifdef USE_STACK_TIMERS
    #pragma omp master
    if (manager)
    {
      if (this == manager->current_timer())
      {
         std::cerr << "Timer loop: " << name << std::endl;
      }
      parent = manager->current_timer();
      if (parent)
      {
        parent->add_child(this);
      }
      manager->push_timer(this);
    }

#endif
    start_time = cpu_clock();
  }

  void stop()
  {
    double elapsed = cpu_clock() - start_time;
    total_time += elapsed;
    num_calls++;


#ifdef USE_STACK_TIMERS
    #pragma omp master
    {
      std::string stack_name = get_stack_name();
      per_stack_total_time[stack_name] += elapsed;
      per_stack_num_calls[stack_name] += 1;

      if (manager)
      {
        manager->pop_timer();
      }
    }
#endif
  }
#endif

#ifdef USE_STACK_TIMERS
  std::string get_stack_name();
#endif

#ifdef USE_STACK_TIMERS
  std::map<std::string, double>& get_per_stack_total_time()
  {
    return per_stack_total_time;
  }
#endif

  inline double    get_total() const
  {
    return total_time;
  }

#ifdef USE_STACK_TIMERS
  inline double get_total(const std::string &stack_name)
  {
    return per_stack_total_time[stack_name];
  }
#endif

  inline long  get_num_calls() const
  {
    return num_calls;
  }

#ifdef USE_STACK_TIMERS
  inline long  get_num_calls(const std::string &stack_name)
  {
    return per_stack_num_calls[stack_name];
  }
#endif

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
#ifdef USE_STACK_TIMERS
  ,manager(NULL), parent(NULL)
#endif
  { }

  void set_name(const std::string& myname)
  {
    name=myname;
  }

  void set_manager(TimerManagerClass *mymanager)
  {
#ifdef USE_STACK_TIMERS
    manager = mymanager;
#endif
  }

#ifdef USE_STACK_TIMERS
  NewTimer *get_parent()
  {
    return parent;
  }
#endif

#ifdef USE_STACK_TIMERS
  std::vector<NewTimer *> &get_children()
  {
    return children;
  }
#endif


#ifdef USE_STACK_TIMERS
  void add_child(NewTimer *t)
  {
    bool found = false;
    for (int i = 0; i < children.size(); i++)
    {
       if (t == children[i])
       {
         found = true;
         break;
       }
    }
    if (!found)
    {
      children.push_back(t);
    }
  }
#endif
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
