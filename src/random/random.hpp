//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file random.h
 * @brief Declare a global Random Number Generator with OHMMS_PRECISION
 *
 * qmcplusplus::Random() returns a random number [0,1) For OpenMP is enabled,
 * it is important to use thread-safe boost::random. Each thread uses its own
 * random number generator with a distinct seed. This prevents a use of global
 * lock which will slow down the applications significantly.
 * c++11 is required
 */
#ifndef QMCPLUSPLUS_RANDOM_MASTER_H
#define QMCPLUSPLUS_RANDOM_MASTER_H
#include "config.h"
#include <cmath>
#include <ctime>
#include <omp.h>

#include <stdint.h>
#include <iostream>

struct BoxMuller2
{ 
  template<typename RNG>
  static inline void generate(RNG& rng, double* restrict a, int n)
  {
    for (int i=0; i+1<n; i+=2) 
    {
      double temp1=1.0-0.9999999999*rng(), temp2=rng();
      a[i]  =sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
      a[i+1]=sqrt(-2.0*log(temp1))*sin(6.283185306*temp2);
    }
    if (n%2==1) {
      double temp1=1-0.9999999999*rng(), temp2=rng();
      a[n-1]=sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
    }
  }

  template<typename RNG>
  static inline void generate(RNG& rng, float* restrict a, int n)
  {
    for (int i=0; i+1<n; i+=2) 
    {
      float temp1=1.0f-0.9999999999f*rng(), temp2=rng();
      a[i]  =sqrtf(-2.0*logf(temp1))*cosf(6.283185306f*temp2);
      a[i+1]=sqrtf(-2.0*logf(temp1))*sinf(6.283185306f*temp2);
    }
    if (n%2==1) {
      float temp1=1.0f-0.9999999999f*rng(), temp2=rng();
      a[n-1]=sqrtf(-2.0f*logf(temp1))*cosf(6.283185306f*temp2);
    }
  }
};


inline uint32_t MakeSeed(int i, int n)
{
  const uint32_t u=1<<10;
  return static_cast<uint32_t>(std::time(nullptr))%u+(i+1)*n+i;
}


#if (__cplusplus< 201103L)
#include "Utilities/RandomGenerator.h"
namespace qmcplusplus
{
  typedef BoostRandom<OHMMS_PRECISION> RandomGenerator_t;
  extern RandomGenerator_t Random;
}
#else
#include "random/BoostRandom.h"
namespace qmcplusplus
{
  template<class T> using RandomGenerator=BoostRandom2<T>;
  typedef RandomGenerator<OHMMS_PRECISION> RandomGenerator_t;
  extern RandomGenerator_t Random;
}
#endif /* if (__cplusplus < 201103L) */

#endif

