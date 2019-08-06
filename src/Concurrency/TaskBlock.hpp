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
 */

#include <functional>
#include <thread>

namespace qmcplusplus
{
template<typename F>
class TaskWrapper
{
  F f;
  template<typename... T>
  void operator()(T&&... args)
  {
    f(std::forward<T>(args)...);
  }
};

enum class Threading
{
  OPENMP,
  STD
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

template<Threading TT>
inline void TaskBlockBarrier<TT>::wait()
{
#pragma omp barrier
}

template<>
inline void TaskBlockBarrier<Threading::STD>::wait()
{
  std::cout << "We'll need to implement this for std::thread\n";
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

template<Threading TT>
template<typename F, typename... Args>
void TaskBlock<TT>::operator()(F&& f, TaskBlockBarrier<TT>& barrier, Args&&... args)

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
} // namespace qmcplusplus

#endif
