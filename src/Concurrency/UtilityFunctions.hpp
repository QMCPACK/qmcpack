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
template<Executor TT = Executor::OPENMP>
class OverrideMaxCapacity;

template<>
class OverrideMaxCapacity<Executor::OPENMP>
{
private:
  int original_max_threads_;

public:
  OverrideMaxCapacity(int max_threads)
  {
    original_max_threads_ = omp_get_max_threads();
    omp_set_num_threads(max_threads);
  }

  ~OverrideMaxCapacity() { omp_set_num_threads(original_max_threads_); }
};

template<Executor TT>
class OverrideMaxCapacity
{
  OverrideMaxCapacity(int max_threads) {}
};

} // namespace Concurrency
} // namespace qmcplusplus

#endif
