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
 *  @brief For threads performing identical independent tasks
 *
 *
 */
#include "Concurrency/Info.hpp"

namespace qmcplusplus
{
/** Abstraction for simple 1 task to 1 thread concurrency
 *
 *  Construct with num_threads to to run on
 *  then call operator(F, args...) 
 *  F is lambda or function with form
 *      void F(int task_id, args...)
 *      [](int task_id, args...){...}
 *
 */
template<Threading TT = Threading::OPENMP>
class TasksOneToOne
{
public:
  TasksOneToOne(unsigned int num_threads) : num_threads_(num_threads) {}
  template<typename F, typename... Args>
  void operator()(F&& f, Args&&... args);

private:
  unsigned int num_threads_;
};

} // namespace qmcplusplus

// Implementation includes must follow functor declaration
#include "Concurrency/TasksOPENMP.hpp"
#include "Concurrency/TasksSTD.hpp"
// Additional implementations enabled by cmake options would go here
#endif
