//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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

#include <string>
#include <algorithm>
#include <map>
#include "config.h"
#include "Clock.h"

#ifdef USE_VTUNE_TASKS
#include <ittnotify.h>
#endif

#ifdef USE_NVTX_API
#include <nvToolsExt.h>
#endif

#define USE_STACK_TIMERS

namespace qmcplusplus
{
template<class TIMER>
class TimerManager;

enum timer_levels
{
  timer_level_none, // The 'none' settting is not for individual timers.
                    // It is for setting a threshold to turn all timers off.
  timer_level_coarse,
  timer_level_medium,
  timer_level_fine,
  num_timer_levels // this is not a timer level but to count the elements in this enum
};

extern bool timer_max_level_exceeded;

// Unsigned char gives 254 timers (0 is reserved).
// Use a longer type (eg. unsigned short) to increase the limit.
using timer_id_t = unsigned char;

// Key for tracking time per stack.  Parametered by size.
template<int N>
class StackKeyParam
{
public:
  // The union is for a performance hack
  // Use the array of small types to store the stack of timer id's.
  // Use the larger types for fast comparison for storage in a map.

  // If timer_id_t is char, there can be up to 254 timers.
  // N is the number of long ints to store timer nesting levels.
  // Each N gives (64 bits/long int) / (8 bits/char) = 8 levels
  union
  {
    long int long_buckets[N];
    timer_id_t short_buckets[sizeof(long int) * N / sizeof(timer_id_t)];
  };

  static const int max_level = sizeof(long int) * N;

  StackKeyParam() : level(0)
  {
    for (int j = 0; j < N; j++)
    {
      long_buckets[j] = 0;
    }
  }

  int level;

  void add_id(timer_id_t c1)
  {
    short_buckets[level] = c1;
    if (level >= max_level - 1)
    {
      timer_max_level_exceeded = true;
    }
    else
    {
      level++;
    }
  }

  void put_id(timer_id_t c1) { short_buckets[level] = c1; }

  timer_id_t get_id(int idx) const { return short_buckets[idx]; }

  bool operator==(const StackKeyParam& rhs)
  {
    bool same = false;
    for (int j = 0; j < N; j++)
    {
      same &= this->long_buckets[j] == rhs.long_buckets[j];
    }
    return same;
  }

  bool operator<(const StackKeyParam& rhs) const
  {
    for (int j = 0; j < N; j++)
    {
      if (!(this->long_buckets[j] == rhs.long_buckets[j]))
      {
        return this->long_buckets[j] < rhs.long_buckets[j];
      }
    }
    return this->long_buckets[N - 1] < rhs.long_buckets[N - 1];
  }
};

// N = 2 gives 16 nesting levels
using StackKey = StackKeyParam<2>;

/** Timer accumulates time and call counts
 * @tparam CLOCK can be CPUClock or FakeCPUClock
 */
template<class CLOCK>
class TimerType
{
protected:
  /// start time of the current measurement
  double start_time;
  /// total time accumulated of all the calls
  double total_time;
  /// total call counts
  long num_calls;
  /// name of this timer
  std::string name;
  /// if false, timer start/stop becomes no-op.
  bool active;
  /// timer level
  timer_levels timer_level;
  /// timer id in registered in the manager
  timer_id_t timer_id;
  /// timer manager which allocated this timer object. nullptr if USE_STACK_TIMERS is not used.
  TimerManager<TimerType<CLOCK>>* manager;
#ifdef USE_STACK_TIMERS
  /// stack key of the current measurement
  StackKey current_stack_key;

  /// total time accumulated per stack key
  std::map<StackKey, double> per_stack_total_time;
  /// total call counts per stack key
  std::map<StackKey, long> per_stack_num_calls;
#endif

#ifdef USE_VTUNE_TASKS
  __itt_string_handle* task_name;
#endif
public:
  void start();
  void stop();

#ifdef USE_STACK_TIMERS
  std::map<StackKey, double>& get_per_stack_total_time() { return per_stack_total_time; }

  StackKey& get_stack_key() { return current_stack_key; }
#endif


  inline double get_total() const { return total_time; }
  inline long get_num_calls() const { return num_calls; }

#ifdef USE_STACK_TIMERS
  inline double get_total(const StackKey& key) { return per_stack_total_time[key]; }
  inline long get_num_calls(const StackKey& key) { return per_stack_num_calls[key]; }
#endif

  timer_id_t get_id() const { return timer_id; }

  void set_id(timer_id_t id) { timer_id = id; }

  inline std::string get_name() const { return name; }

  inline void reset()
  {
    num_calls  = 0;
    total_time = 0.0;
  }

  TimerType(const std::string& myname,
            TimerManager<TimerType<CLOCK>>* mymanager,
            timer_levels mytimer = timer_level_fine)
      : total_time(0.0),
        num_calls(0),
        name(myname),
        active(true),
        timer_level(mytimer),
        timer_id(0),
#ifdef USE_STACK_TIMERS
        manager(mymanager)
#else
        manager(nullptr)
#endif
  {
#ifdef USE_VTUNE_TASKS
    task_name = __itt_string_handle_create(myname.c_str());
#endif
  }

  TimerType(const TimerType& o) = delete;

  void set_name(const std::string& myname) { name = myname; }

  void set_active(const bool& is_active) { active = is_active; }

  void set_active_by_timer_threshold(const timer_levels threshold);

  // Functions for unit testing
  template<class CLOCK1>
  friend void set_total_time(TimerType<CLOCK1>* timer, double total_time_input);

  template<class CLOCK1>
  friend void set_num_calls(TimerType<CLOCK1>* timer, long num_calls_input);
};

using NewTimer  = TimerType<CPUClock>;
using FakeTimer = TimerType<FakeCPUClock>;
extern template class TimerType<CPUClock>;
extern template class TimerType<FakeCPUClock>;

// Wrapper for timer that starts on construction and stops on destruction
template<class TIMER = NewTimer>
class ScopeGuard
{
public:
  ScopeGuard(TIMER& t) : timer(t) { timer.start(); }

  ScopeGuard(const ScopeGuard&)                = delete;
  ScopeGuard& operator=(const ScopeGuard&)     = delete;
  ScopeGuard(ScopeGuard&&) noexcept            = default;
  ScopeGuard& operator=(ScopeGuard&&) noexcept = default;

  ~ScopeGuard() { timer.stop(); }

private:
  TIMER& timer;
};

using ScopedTimer     = ScopeGuard<NewTimer>;
using ScopedFakeTimer = ScopeGuard<FakeTimer>;

} // namespace qmcplusplus

#endif
