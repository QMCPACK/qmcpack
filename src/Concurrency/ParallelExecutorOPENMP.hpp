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
 *  @brief implementation of openmp specialization of ParallelExecutor
 */
#ifndef QMCPLUSPLUS_PARALLELEXECUTOR_OPENMP_HPP
#define QMCPLUSPLUS_PARALLELEXECUTOR_OPENMP_HPP

#include <atomic>
#include <stdexcept>
#include <string>
#include <exception>
#include <iostream>

#include "Concurrency/ParallelExecutor.hpp"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
/** Implements parallel tasks executed by threads in an OpenMP thread pool.
 *
 *  Nested OpenMP execution is not allowed; calling this when omp_get_level() > 0
 *  throws a std::runtime_error.
 *
 *  Exceptions thrown by tasks, including those not derived from std::exception, are caught
 *  at the thread level to avoid std::terminate. The first caught exception is captured, prepended with
 *  the task ID, and rethrown after the parallel region completes. Once an exception
 *  occurs, subsequent loop iterations are skipped. Any later exceptions caught from
 *  concurrent threads are output to std::cerr.
 */
template<>
template<typename F, typename... Args>
void ParallelExecutor<Executor::OPENMP>::operator()(int num_tasks, F&& f, Args&&... args)
{
  const std::string nesting_error{"ParallelExecutor should not be used for nested openmp threading\n"};
  if (omp_get_level() > 0)
    throw std::runtime_error(nesting_error);
  std::exception_ptr first_exception;
  // lock-free flag so threads can check for failure without racing on first_exception
  std::atomic<bool> exception_occurred{false};

  // Record the first failure to be rethrown later, report any further ones to std::cerr.
  auto record_exception = [&](const std::string& msg) {
#pragma omp critical
    if (!first_exception)
    {
      first_exception = std::make_exception_ptr(std::runtime_error(msg));
      exception_occurred.store(true, std::memory_order_release);
    }
    else
    {
      std::cerr << msg;
      if (msg.empty() || msg.back() != '\n')
        std::cerr << '\n';
    }
  };

#pragma omp parallel for
  for (int task_id = 0; task_id < num_tasks; ++task_id)
  {
    if (exception_occurred.load(std::memory_order_acquire))
      continue;
    const std::string error_prefix = "Unexpected exception thrown in threaded section: task_id " + std::to_string(task_id);
    try
    {
      f(task_id, std::forward<Args>(args)...);
    }
    catch (const std::exception& e)
    {
      record_exception(error_prefix + " " + e.what());
    }
    catch (...)
    {
      record_exception(error_prefix + " unknown exception type not derived from std::exception");
    }
  }
  if (exception_occurred.load(std::memory_order_acquire))
    std::rethrow_exception(first_exception);
}

} // namespace qmcplusplus

#endif
