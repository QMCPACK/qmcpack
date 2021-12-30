//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file RandomGenerator.h
 * @brief Declare a global Random Number Generator
 *
 * Selected among
 * - boost::random
 * - sprng
 * - math::random
 * qmcplusplus::Random() returns a random number [0,1)
 * For OpenMP is enabled, it is important to use thread-safe boost::random. Each
 * thread uses its own random number generator with a distinct seed. This prevents
 * a use of global lock which will slow down the applications significantly.
 */
#ifndef OHMMS_RANDOMGENERATOR
#define OHMMS_RANDOMGENERATOR
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <cmath>
#include <ctime>

#include <stdint.h>

uint32_t make_seed(int i, int n);

// The definition of the fake RNG should always be available for unit testing
#include "Utilities/FakeRandom.h"
#ifdef USE_FAKE_RNG
namespace qmcplusplus
{
typedef FakeRandom RandomGenerator_t;

} // namespace qmcplusplus
#else

#ifdef HAVE_LIBBOOST

#include "Utilities/BoostRandom.h"
namespace qmcplusplus
{
template<class T>
using RandomGenerator   = BoostRandom<T>;
using RandomGenerator_t = BoostRandom<OHMMS_PRECISION_FULL>;
} // namespace qmcplusplus
#else

#error -DHAVE_LIBBOOST is missing in the compile line. A cmake dependency fix is needed.

#endif
#endif


namespace qmcplusplus
{
class RNGThreadSafe : public RandomGenerator_t
{
public:
  inline RandomGenerator_t::result_type rand()
  {
    result_type result;
// This should be a named section but at least clang 9 doesn't seem to support
// and warns of extra tokens.
#pragma omp critical
    {
      result = RandomGenerator_t::rand();
    }
    return result;
  }

  /** return a random number [0,1)
   */
  inline RandomGenerator_t::result_type operator()()
  {
    result_type result;
#pragma omp critical
    {
      result = RandomGenerator_t::rand();
    }
    return result;
  }

  /** generate a series of random numbers */
  template<typename T1>
  inline void generate_uniform(T1* restrict d, int n)
  {
#pragma omp critical
    {
      for (int i = 0; i < n; ++i)
        d[i] = RandomGenerator_t::rand();
    }
  }
};

extern RNGThreadSafe Random;

} // namespace qmcplusplus

#endif
