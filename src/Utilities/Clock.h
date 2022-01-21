//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CLOCK_H
#define QMCPLUSPLUS_CLOCK_H

#include "Concurrency/OpenMP.h"

#include <sys/time.h>
#include <stddef.h>

namespace qmcplusplus
{
/** functor for fake clock
 * calling FakeCPUClock()() returns the clock value
 */
class FakeCPUClock
{
public:
  static double fake_cpu_clock_value;
  static double fake_cpu_clock_increment;

  double operator()()
  {
    fake_cpu_clock_value += fake_cpu_clock_increment;
    return fake_cpu_clock_value;
  }
};

/** functor for high precision clock
 * calling CPUClock()() returns the clock value
 */
class CPUClock
{
public:
  double operator()()
  {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (1.e-6) * tv.tv_usec;
#endif
  }
};

} // namespace qmcplusplus
#endif
