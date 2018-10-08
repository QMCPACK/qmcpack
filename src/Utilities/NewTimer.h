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
#ifdef USE_VTUNE_TASKS
#include <ittnotify.h>
#endif

#define USE_STACK_TIMERS

class Communicate;

namespace qmcplusplus
{

class NewTimer;

enum timer_levels {
  timer_level_none,   // The 'none' settting is not for individual timers.
                      // It is for setting a threshold to turn all timers off.
  timer_level_coarse,
  timer_level_medium,
  timer_level_fine
};

const char TIMER_STACK_SEPARATOR='/';

// Unsigned char gives 254 timers (0 is reserved).
// Use a longer type (eg. unsigned short) to increase the limit.
typedef unsigned char timer_id_t;

extern bool timer_max_level_exceeded;

// Key for tracking time per stack.  Parametered by size.

template<int N> class StackKeyParam
{
public:
    // The union is for a performance hack
    // Use the array of small types to store the stack of timer id's.
    // Use the larger types for fast comparison for storage in a map.

    // If timer_id_t is char, there can be up to 254 timers.
    // N is the number of long ints to store timer nesting levels.
    // Each N gives (64 bits/long int) / (8 bits/char) = 8 levels
    union {
        long int long_buckets[N];
        timer_id_t short_buckets[sizeof(long int)*N/sizeof(timer_id_t)];
    };

    static const int max_level = sizeof(long int)*N;

    StackKeyParam() : level(0) {
      for (int j = 0; j < N; j++) {
        long_buckets[j] = 0;
      }
    }

    int level;

    void add_id(timer_id_t c1)
    {
      short_buckets[level] = c1;
      if (level >= max_level-1) {
        timer_max_level_exceeded = true;
      } else {
        level++;
      }
    }

    void put_id(timer_id_t c1)
    {
        short_buckets[level] = c1;
    }

    timer_id_t get_id(int idx) const
    {
      return short_buckets[idx];
    }

    bool operator==(const StackKeyParam &rhs)
    {
        bool same = false;
        for (int j = 0; j < N; j++)
        {
            same &= this->long_buckets[j] ==  rhs.long_buckets[j];
        }
        return same;
    }

    bool operator<(const StackKeyParam &rhs) const
    {
        for (int j = 0; j < N; j++)
        {
            if (!(this->long_buckets[j] == rhs.long_buckets[j]))
            {
              return this->long_buckets[j] < rhs.long_buckets[j];
            }
        }
        return this->long_buckets[N-1] < rhs.long_buckets[N-1];
    }
};

// N = 2 gives 16 nesting levels
typedef StackKeyParam<2> StackKey;

class TimerManagerClass
{
protected:
  std::vector<NewTimer*> TimerList;
  std::vector<NewTimer*> CurrentTimerStack;
  timer_levels timer_threshold;
  timer_id_t max_timer_id;
  bool max_timers_exceeded;
  std::map<timer_id_t, std::string> timer_id_name;
  std::map<std::string, timer_id_t> timer_name_to_id;
public:
#ifdef USE_VTUNE_TASKS
  __itt_domain *task_domain;
#endif

  TimerManagerClass():timer_threshold(timer_level_coarse),max_timer_id(1),
    max_timers_exceeded(false) {
#ifdef USE_VTUNE_TASKS
      task_domain = __itt_domain_create("QMCPACK");
#endif
    }
  void addTimer (NewTimer* t);
  NewTimer *createTimer(const std::string& myname, timer_levels mytimer = timer_level_fine);

  void push_timer(NewTimer *t)
  {
    {
      CurrentTimerStack.push_back(t);
    }
  }

  void pop_timer()
  {
    {
      CurrentTimerStack.pop_back();
    }
  }

  NewTimer *current_timer()
  {
    NewTimer *current = NULL;
    if (CurrentTimerStack.size() > 0)
    {
       current = CurrentTimerStack.back();
    }
    return current;
  }

  void set_timer_threshold(const timer_levels threshold);

  bool maximum_number_of_timers_exceeded() const
  {
    return max_timers_exceeded;
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

  void get_stack_name_from_id(const StackKey &key, std::string &name);

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
  bool active;
  timer_levels timer_level;
  timer_id_t timer_id;
#ifdef USE_STACK_TIMERS
  TimerManagerClass *manager;
  NewTimer *parent;
  StackKey current_stack_key;

  std::map<StackKey, double> per_stack_total_time;
  std::map<StackKey, long> per_stack_num_calls;
#endif

#ifdef USE_VTUNE_TASKS
  __itt_string_handle *task_name;
#endif
public:
#if not(ENABLE_TIMERS)
  inline void start() {}
  inline void stop() {}
#else
  void start()
  {
    if (active)
    {
#ifdef USE_STACK_TIMERS

#ifdef USE_VTUNE_TASKS
      __itt_id parent_task = __itt_null;
      __itt_task_begin(manager->task_domain, __itt_null, parent_task, task_name);
#endif

      #pragma omp master
      {
        if (manager)
        {
          if (this == manager->current_timer())
          {
             std::cerr << "Timer loop: " << name << std::endl;
          }
          if (parent != manager->current_timer())
          {
            parent = manager->current_timer();
            if (parent)
            {
              current_stack_key = parent->get_stack_key();
              current_stack_key.add_id(timer_id);
            }
          }
          if (parent == NULL)
          {
            current_stack_key = StackKey();
            current_stack_key.add_id(timer_id);
          }
          manager->push_timer(this);
        }
        start_time = cpu_clock();
      }
#else
      start_time = cpu_clock();
#endif
    }
  }

  void stop()
  {
    if (active)
    {
#ifdef USE_STACK_TIMERS

#ifdef USE_VTUNE_TASKS
      __itt_task_end(manager->task_domain);
#endif

      #pragma omp master
#endif
      {
        double elapsed = cpu_clock() - start_time;
        total_time += elapsed;
        num_calls++;

#ifdef USE_STACK_TIMERS
        per_stack_total_time[current_stack_key] += elapsed;
        per_stack_num_calls[current_stack_key] += 1;

        if (manager)
        {
          manager->current_timer()->set_parent(NULL);
          manager->pop_timer();
        }
#endif
      }
    }
  }
#endif

#ifdef USE_STACK_TIMERS
  std::map<StackKey, double>& get_per_stack_total_time()
  {
    return per_stack_total_time;
  }

  StackKey &get_stack_key()
  {
    return current_stack_key;
  }
#endif



  inline double    get_total() const
  {
    return total_time;
  }

#ifdef USE_STACK_TIMERS
  inline double get_total(const StackKey &key)
  {
    return per_stack_total_time[key];
  }
#endif

  inline long  get_num_calls() const
  {
    return num_calls;
  }

#ifdef USE_STACK_TIMERS
  inline long  get_num_calls(const StackKey &key)
  {
    return per_stack_num_calls[key];
  }
#endif

  timer_id_t get_id() const
  {
    return timer_id;
  }

  void set_id(timer_id_t id)
  {
    timer_id = id;
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

  NewTimer(const std::string& myname, timer_levels mytimer = timer_level_fine) :
    total_time(0.0), num_calls(0), name(myname), active(true), timer_level(mytimer)
    ,timer_id(0)
#ifdef USE_STACK_TIMERS
  ,manager(NULL), parent(NULL)
#endif
  {
#ifdef USE_VTUNE_TASKS
    task_name = __itt_string_handle_create(myname.c_str());
#endif
  }



  void set_name(const std::string& myname)
  {
    name=myname;
  }

  void set_active(const bool &is_active)
  {
    active = is_active;
  }

  void set_active_by_timer_threshold(const timer_levels threshold);

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

  void set_parent(NewTimer *new_parent)
  {
    parent = new_parent;
  }
#endif
};

// Wrapper for timer that starts on construction and stops on destruction
class ScopedTimer
{
public:
  ScopedTimer(NewTimer *t) : timer(t)
  {
    if (timer) timer->start();
  }

  ~ScopedTimer()
  {
    if (timer) timer->stop();
  }
private:
  NewTimer *timer;
};

// Helpers to make it easier to define a set of timers
// See tests/test_timer.cpp for an example


typedef std::vector<NewTimer *> TimerList_t;

template <class T> struct TimerIDName_t
{
  T id;
  const std::string name;
};

// C++ 11 type aliasing
#if __cplusplus >= 201103L
template <class T> using TimerNameList_t = std::vector<TimerIDName_t<T>>;

template <class T> void setup_timers(TimerList_t &timers, TimerNameList_t<T> timer_list,
                                     timer_levels timer_level = timer_level_fine)
{
  timers.resize(timer_list.size());
  for (int i = 0; i < timer_list.size(); i++)
  {
    timers[timer_list[i].id] = TimerManager.createTimer(timer_list[i].name, timer_level);
  }
}
#endif


struct TimerComparator
{
  inline bool operator()(const NewTimer *a, const NewTimer *b)
  {
    return a->get_name() < b->get_name();
  }
};


}

#endif
