////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by:
// Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by:
// Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_THREAD_HPP
#define QMCPLUSPLUS_THREAD_HPP
/** @file
 *  @brief Abstractions threads performing identical independent tasks
 *
 *  We assume they will come to a common synchronization point
 *  TaskWrapper code from
 *  https://stackoverflow.com/questions/48268322/wrap-stdthread-call-function
 *
 *  Currently reference arguments to the "task" need to be wrapped in std::ref
 *  to survive std::forward used in the STD specialization.  This doesn't cause trouble
 *  for the OPENMP variant except for the TaskBlockBarrier which breaks with std::ref.
 */
#include <iostream>
#include <functional>
#include <type_traits>
#include <thread>
#include <utility>

#include <omp.h>

#include "Concurrency/Info.hpp"

namespace qmcplusplus
{

template<typename F>
struct TaskWrapper
{
  F f;

  template<typename... T>
  void operator()(T&&... args)
  {
    f(std::forward<T>(args)...);
  }
};
    
/** Abstraction for simple 1 task to 1 thread concurrency
 *
 *  Functor takes num_threads to construct
 *  Call with F and args...
 *  F takes task_id, args..
 */
template<Threading TT>
class TasksOneToOne
{
public:
  TasksOneToOne(unsigned int num_threads) : num_threads_(num_threads) {}
  template<typename F, typename... Args>
  void operator()(F&& f, Args&&... args);

private:
  unsigned int num_threads_;
};

template<>
template<typename F, typename... Args>
void TasksOneToOne<Threading::OPENMP>::operator()(F&& f, Args&&... args)
{
#pragma omp parallel num_threads(num_threads_)
  {
      f(omp_get_thread_num(), std::forward<Args>(args)...);
  }
}

template<>
template<typename F, typename... Args>
void TasksOneToOne<Threading::STD>::operator()(F&& f, Args&&... args)

{
  std::vector<std::thread> threads(num_threads_);

  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
      threads[task_id] = std::thread(TaskWrapper<F>{std::forward<F>(f)}, task_id, std::forward<Args>(args)...);
  }

  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
    threads[task_id].join();
  }
}


} // namespace qmcplusplus

#endif
