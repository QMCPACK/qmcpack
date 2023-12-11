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


/** @file NewTimer.cpp
 * @brief Implements NewTimer member functions
 */
#include "NewTimer.h"
#include <iostream>
#include "Concurrency/OpenMP.h"
#include "config.h"
#include "TimerManager.h"

namespace qmcplusplus
{
bool timer_max_level_exceeded = false;

#ifndef ENABLE_TIMERS
template<class CLOCK>
void TimerType<CLOCK>::start()
{}
template<class CLOCK>
void TimerType<CLOCK>::stop()
{}
#else
template<class CLOCK>
void TimerType<CLOCK>::start()
{
  if (active)
  {
#ifdef USE_STACK_TIMERS

#ifdef USE_VTUNE_TASKS
    __itt_id parent_task = __itt_null;
    __itt_task_begin(manager->task_domain, __itt_null, parent_task, task_name);
#endif

#ifdef USE_NVTX_API
    nvtxRangePushA(name.c_str());
#endif

    bool is_true_master(true);
    for (int level = omp_get_level(); level > 0; level--)
      if (omp_get_ancestor_thread_num(level) != 0)
        is_true_master = false;
    if (is_true_master)
    {
      if (manager)
      {
        // compute current_stack_key
        TimerType* parent = manager->current_timer();
        if (parent)
        {
          current_stack_key = parent->get_stack_key();
          current_stack_key.add_id(timer_id);
        }
        else
        {
          current_stack_key = StackKey();
          current_stack_key.add_id(timer_id);
        }

        manager->push_timer(this);
      }
      start_time = CLOCK::now();
    }
#else
    start_time                            = CLOCK::now();
#endif
  }
}

template<class CLOCK>
void TimerType<CLOCK>::stop()
{
  if (active)
  {
#ifdef USE_STACK_TIMERS

#ifdef USE_VTUNE_TASKS
    __itt_task_end(manager->task_domain);
#endif

#ifdef USE_NVTX_API
    nvtxRangePop();
#endif

    bool is_true_master(true);
    for (int level = omp_get_level(); level > 0; level--)
      if (omp_get_ancestor_thread_num(level) != 0)
        is_true_master = false;
    if (is_true_master)
    {
      std::chrono::duration<double> elapsed = CLOCK::now() - start_time;
      total_time += elapsed.count();
      num_calls++;

      per_stack_total_time[current_stack_key] += elapsed.count();
      per_stack_num_calls[current_stack_key] += 1;

      if (manager)
        manager->pop_timer(this);
    }
#else
    std::chrono::duration<double> elapsed = CLOCK::now() - start_time;
    total_time += elapsed.count();
    num_calls++;
#endif
  }
}
#endif

template<class CLOCK>
void TimerType<CLOCK>::set_active_by_timer_threshold(const timer_levels threshold)
{
  if (timer_level <= threshold)
    active = true;
  else
    active = false;
}

template class TimerType<ChronoClock>;
template class TimerType<FakeChronoClock>;

} // namespace qmcplusplus
