#ifndef QMCPLUSPLUS_CLOCK_H
#define QMCPLUSPLUS_CLOCK_H

#include "Message/OpenMP.h"

#if defined(__bgp__)
#include <dcmf.h>
#endif

namespace qmcplusplus
{
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
}
#endif
