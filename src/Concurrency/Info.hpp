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


#ifndef QMCPLUSPLUS_CONCURRENCY_INFO_HPP
#define QMCPLUSPLUS_CONCURRENCY_INFO_HPP
/** @file
 *  @brief Abstraction of information on executor environments
 */

#include "OpenMP.h"

#ifdef QMC_EXP_THREADING
#include <thread>
#endif

namespace qmcplusplus
{
enum class Executor
{
  OPENMP,
#ifdef QMC_EXP_THREADING
  STD_THREADS
#endif
};

namespace Concurrency
{
using qmcplusplus::Executor;

template<Executor TT = Executor::OPENMP>
unsigned int maxCapacity();

template<>
inline unsigned int maxCapacity<Executor::OPENMP>()
{
  return omp_get_max_threads();
}

template<Executor TT = Executor::OPENMP>
unsigned int getWorkerId();

template<>
inline unsigned int getWorkerId<Executor::OPENMP>()
{
  return omp_get_thread_num();
}

#ifdef QMC_EXP_THREADING
template<>
inline unsigned int maxCapacity<Executor::STD_THREADS>()
{
  // Does taskset fix what this reports?  i.e. deal with binding to socket properly
  return std::thread::hardware_concurrency();
}
#endif

} // namespace Concurrency
} // namespace qmcplusplus
#endif
