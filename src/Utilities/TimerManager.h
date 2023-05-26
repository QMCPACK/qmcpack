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


/** @file TimerManager.h
 * @brief timer_manager class.
 */
#ifndef QMCPLUSPLUS_TIMER_MANAGER_H
#define QMCPLUSPLUS_TIMER_MANAGER_H

#include <vector>
#include <string>
#include <mutex>
#include <map>
#include <memory>
#include <type_traits>
#include "NewTimer.h"
#include "config.h"
#include "OhmmsData/Libxml2Doc.h"

#ifdef USE_VTUNE_TASKS
#include <ittnotify.h>
#endif

class Communicate;

namespace qmcplusplus
{
/** Manager creates timers and handle reports
 * @tparam TIMER regular or fake timer
 *
 * TimerManager is generally not thread-safe.
 * Thread-safe functions are noted below.
 */
template<class TIMER>
class TimerManager
{
private:
  /// All the timers created by this manager
  std::vector<std::unique_ptr<TIMER>> timer_storage_;
  /// mutex for TimerList
  std::mutex timer_list_lock_;
  /// The stack of nested active timers
  std::vector<TIMER*> CurrentTimerStack;
  /// The threshold for active timers
  timer_levels timer_threshold;
  /// The current maximal timer id
  timer_id_t max_timer_id;
  /// status of maxmal timer id reached
  bool max_timers_exceeded;
  /// timer id to name mapping
  std::map<timer_id_t, std::string> timer_id_name;
  /// name to timer id mapping
  std::map<std::string, timer_id_t> timer_name_to_id;

  void initializeTimer(TIMER& t);

  void print_flat(Communicate* comm);
  void print_stack(Communicate* comm);

public:
#ifdef USE_VTUNE_TASKS
  __itt_domain* task_domain;
#endif

  TimerManager() : timer_threshold(timer_level_coarse), max_timer_id(1), max_timers_exceeded(false)
  {
#ifdef USE_VTUNE_TASKS
    task_domain = __itt_domain_create("QMCPACK");
#endif
  }

  /// Create a new timer object registred in this manager. This call is thread-safe.
  TIMER* createTimer(const std::string& myname, timer_levels mytimer = timer_level_fine);

  void push_timer(TIMER* t);

  void pop_timer(TIMER* t);

  TIMER* current_timer()
  {
    TIMER* current = nullptr;
    if (CurrentTimerStack.size() > 0)
      current = CurrentTimerStack.back();

    return current;
  }

  void set_timer_threshold(const timer_levels threshold);
  void set_timer_threshold(const std::string& threshold);
  std::string get_timer_threshold_string() const;

  bool maximum_number_of_timers_exceeded() const { return max_timers_exceeded; }

  void reset();
  void print(Communicate* comm);

  using nameList_t = std::map<std::string, int>;
  using timeList_t = std::vector<double>;
  using callList_t = std::vector<long>;
  using names_t    = std::vector<std::string>;

  struct FlatProfileData
  {
    nameList_t nameList;
    timeList_t timeList;
    callList_t callList;
  };

  struct StackProfileData
  {
    names_t names;
    nameList_t nameList;
    timeList_t timeList;
    timeList_t timeExclList;
    callList_t callList;
  };

  void collate_flat_profile(Communicate* comm, FlatProfileData& p);

  void collate_stack_profile(Communicate* comm, StackProfileData& p);

  void output_timing(Communicate* comm, Libxml2Document& doc, xmlNodePtr root);

  void get_stack_name_from_id(const StackKey& key, std::string& name);
};

extern template class TimerManager<NewTimer>;
extern template class TimerManager<FakeTimer>;

TimerManager<NewTimer>& getGlobalTimerManager();

NewTimer& createGlobalTimer(const std::string& myname, timer_levels mylevel = timer_level_fine);

// Helpers to make it easier to define a set of timers
// See tests/test_timer.cpp for an example

template<class T>
struct TimerIDName_t
{
  T id;
  const std::string name;
};

template<class T>
using TimerNameList_t = std::vector<TimerIDName_t<T>>;

template<class TIMER>
class TimerList : public std::vector<std::reference_wrapper<TIMER>>
{
public:
  template<class T>
  TimerList(TimerManager<TIMER>& manager,
            const TimerNameList_t<T>& timer_list,
            timer_levels timer_level = timer_level_fine)
  {
    this->reserve(timer_list.size());
    for (std::size_t i = 0; i < timer_list.size(); i++)
    {
      if (i != static_cast<std::underlying_type_t<T>>(timer_list[i].id))
        throw std::runtime_error("Mismatch between index and enumeration");
      this->push_back(*manager.createTimer(timer_list[i].name, timer_level));
    }
  }
};

using TimerList_t = TimerList<NewTimer>;

} // namespace qmcplusplus
#endif
