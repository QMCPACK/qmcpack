//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CLOCK_H
#define QMCPLUSPLUS_CLOCK_H

#include "Message/OpenMP.h"

#include <sys/time.h>
#include <stddef.h>

namespace qmcplusplus
{
#if defined(USE_FAKE_CLOCK)
extern double fake_cpu_clock_value;
extern double fake_cpu_clock_increment;
inline double fake_cpu_clock()
{
  fake_cpu_clock_value += fake_cpu_clock_increment;
  return fake_cpu_clock_value;
}
#define cpu_clock fake_cpu_clock
#elif defined(ENABLE_OPENMP)
inline double cpu_clock() { return omp_get_wtime(); }
#else
inline double cpu_clock()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (1.e-6) * tv.tv_usec;
}
#endif
} // namespace qmcplusplus
#endif
