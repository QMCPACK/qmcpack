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

#ifndef QMCPLUSPLUS_CONCURRENCY_INFO_HPP
#define QMCPLUSPLUS_CONCURRENCY_INFO_HPP
/** @file
 *  @brief Abstraction of information on threading environment
 */

#include <omp.h>
#include <thread>

namespace qmcplusplus
{
enum class Threading
{
  OPENMP,
  STD
};

namespace Concurrency
{
using Threading = typename qmcplusplus::Threading;
template<Threading TT = Threading::OPENMP>
unsigned int maxThreads();

template<>
inline unsigned int maxThreads<Threading::OPENMP>()
{
  return omp_get_max_threads();
}

template<>
inline unsigned int maxThreads<Threading::STD>()
{
  // Does taskset fix what this reports?  i.e. deal with binding to socket properly
  return std::thread::hardware_concurrency();
}

} // namespace Concurrency
} // namespace qmcplusplus
#endif
