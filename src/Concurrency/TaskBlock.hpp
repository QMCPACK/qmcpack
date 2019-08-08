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
#include <condition_variable>
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

template<Threading TT>
class TaskBlockBarrier
{
public:
  TaskBlockBarrier(unsigned int num_threads) : num_threads_(num_threads) {}
  void wait();
  
private:
  unsigned int num_threads_;
};

template<>
class TaskBlockBarrier<Threading::OPENMP>
{
public:
  TaskBlockBarrier(unsigned int num_threads) : num_threads_(num_threads)
  {
    omp_init_lock(&mutex_);
    omp_init_lock(&barrier_);
  }
  void wait()
  {
  omp_set_lock(&mutex_);
  threads_at_barrier_++;
  int local_tab = threads_at_barrier_;
  omp_unset_lock(&mutex_);
  omp_set_lock(&barrier_);
  while(local_tab < num_threads_)
  {
    omp_set_lock(&mutex_);
    local_tab = threads_at_barrier_;
    omp_unset_lock(&mutex_);
  }
  omp_unset_lock(&barrier_);
  }
  
private:
  omp_lock_t mutex_;
  omp_lock_t barrier_;
  unsigned int num_threads_;
  int threads_at_barrier_;
};

    
template<>
class TaskBlockBarrier<Threading::STD>
{
public:
  TaskBlockBarrier(unsigned int num_threads) : num_threads_(num_threads),
						 threads_at_barrier_(0)  {}
  void wait();


private:
  std::condition_variable barrier_;
  int num_threads_;
  std::mutex mutex_;
  int threads_at_barrier_;
};
    

void TaskBlockBarrier<Threading::STD>::wait()
{
    {
	std::lock_guard<std::mutex> lock(mutex_);
	threads_at_barrier_++;
    }
    if (threads_at_barrier_ == num_threads_)
    {
	barrier_.notify_all();
    }
    else
    {
	std::unique_lock<std::mutex> lock(mutex_);
	barrier_.wait(lock);
    }
}

template<Threading TT>
class TaskBlock
{
public:
  TaskBlock(unsigned int num_threads) : num_threads_(num_threads) {}
  template<typename F, typename... Args>
  void operator()(F&& f, Args&&... args);
  template<typename F, typename... Args>
  void operator()(F&& f, TaskBlockBarrier<TT>& barrier, Args&&... args);

private:
  unsigned int num_threads_;
};

template<Threading TT>
template<typename F, typename... Args>
void TaskBlock<TT>::operator()(F&& f, Args&&... args)

{
#pragma omp parallel for
  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
    f(task_id, std::forward<Args>(args)...);
  }
}

template<>
template<typename F, typename... Args>
void TaskBlock<Threading::OPENMP>::operator()(F&& f, TaskBlockBarrier<Threading::OPENMP>& barrier, Args&&... args)
{
  omp_set_num_threads(num_threads_);
#pragma omp parallel for
  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
    f(task_id, barrier, std::forward<Args>(args)...);
  }
}

template<>
template<typename F, typename... Args>
void TaskBlock<Threading::STD>::operator()(F&& f, Args&&... args)

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

template<>
template<typename F, typename... Args>
void TaskBlock<Threading::STD>::operator()(F&& f, TaskBlockBarrier<Threading::STD>& barrier, Args&&... args)

{
  std::vector<std::thread> threads(num_threads_);

  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
      threads[task_id] = std::thread(TaskWrapper<F>{std::forward<F>(f)}, task_id, barrier, std::forward<Args>(args)...);
  }

  for (int task_id = 0; task_id < num_threads_; ++task_id)
  {
    threads[task_id].join();
  }
}

} // namespace qmcplusplus

#endif
