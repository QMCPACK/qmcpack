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

#ifndef QMCPLUSPLUS_UTILITYFUNCTIONS_HPP
#define QMCPLUSPLUS_UTILITYFUNCTIONS_HPP
/** @file
 *  @brief utility functions for executors
 */

#include "Concurrency/Info.hpp"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
namespace Concurrency
{
/** A service class to restore active avaiable threads upon destruction as the thread count recorded during construction.
 */
template<Executor TT = Executor::OPENMP>
class ThreadCountProtector
{
  ThreadCountProtector() {}
};

template<>
class ThreadCountProtector<Executor::OPENMP>
{
private:
  const int original_max_threads_;

public:
  ThreadCountProtector() : original_max_threads_(omp_get_max_threads()) {}

  ~ThreadCountProtector() { omp_set_num_threads(original_max_threads_); }
};

/** A service class to override active avaiable threads upon construction.
 * Restore active avaiable threads upon destruction as the thread count recorded during construction.
 */
template<Executor TT = Executor::OPENMP>
class OverrideMaxCapacity : private ThreadCountProtector<TT>
{
  OverrideMaxCapacity(int max_threads) {}
};

template<>
class OverrideMaxCapacity<Executor::OPENMP> : private ThreadCountProtector<Executor::OPENMP>
{
public:
  OverrideMaxCapacity(int max_threads) { omp_set_num_threads(max_threads); }
};

} // namespace Concurrency
} // namespace qmcplusplus

#endif
