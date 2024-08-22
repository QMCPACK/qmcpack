////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
////////////////////////////////////////////////////////////////////////////////


/** @file
 *  @brief implementation of std::thread specialization of ParallelExecutor
 */
#ifndef QMCPLUSPLUS_PARALLELEXECUTOR_STDTHREADS_HPP
#define QMCPLUSPLUS_PARALLELEXECUTOR_STDTHREADS_HPP

#include <vector>
#include <functional>
#include <utility>
#include <thread>

#include "Concurrency/ParallelExecutor.hpp"

namespace qmcplusplus
{
/** implements parallel tasks executed by STD threads. One task one thread mapping.
 */
template<>
template<typename F, typename... Args>
void ParallelExecutor<Executor::STD_THREADS>::operator()(int num_tasks, F&& f, Args&&... args)
{
  std::vector<std::thread> threads(num_tasks);

  for (int task_id = 0; task_id < num_tasks; ++task_id)
  {
    threads[task_id] = std::thread(f, task_id, std::ref(std::forward<Args>(args))...);
  }

  for (int task_id = 0; task_id < num_tasks; ++task_id)
  {
    threads[task_id].join();
  }
}

} // namespace qmcplusplus

#endif
