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

#if defined(__bgp__)
#include <dcmf.h>
#endif

namespace qmcplusplus
{
#if defined(USE_FAKE_CLOCK)
extern double fake_cpu_clock_value;
extern double fake_cpu_clock_increment;
double fake_cpu_clock()
{
    fake_cpu_clock_value += fake_cpu_clock_increment;
    return fake_cpu_clock_value;
}
#define cpu_clock fake_cpu_clock
#else
#if defined(__bgp__)

inline double cpu_clock()
{
  return DCMF_Timer();
}
//  __inline__ unsigned long long getticks(void)
//  {
//    unsigned long long int result=0;
//    unsigned long int upper, lower,tmp;
//    __asm__ volatile(
//        "0:                  \n"
//        "\tmftbu   %0           \n"
//        "\tmftb    %1           \n"
//        "\tmftbu   %2           \n"
//        "\tcmpw    %2,%0        \n"
//        "\tbne     0b         \n"
//        : "=r"(upper),"=r"(lower),"=r"(tmp)
//        );
//    result = upper;
//    result = result<<32;
//    result = result|lower;
//
//    return(result);
//  }
//
//  inline double cpu_clock()
//  {
//    //using 512e6 ticks per second
//    const double SEC_PER_TICKS=1.953125e-9;
//    return static_cast<double>(getticks())*SEC_PER_TICKS;
//  }

#else
#if defined(ENABLE_OPENMP)
inline double cpu_clock()
{
  return omp_get_wtime();
}
#else
inline double cpu_clock()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec+(1.e-6)*tv.tv_usec;
}
#endif //
#endif
#endif
}
#endif
