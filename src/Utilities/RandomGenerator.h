//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file RandomGenerator.h
 * @brief Declare a global Random Number Generator
 *
 * Selected among
 * - std::mt19937
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
#include <cstdint>
// The definition of the fake RNG should always be available for unit testing
#include "FakeRandom.h"
#include "StdRandom.h"

uint32_t make_seed(int i, int n);

namespace qmcplusplus
{
template<class RNG>
class RNGThreadSafe : public RNG
{
public:
  using result_type = typename RNG::result_type;

  /** return a random number [0,1)
   */
  result_type operator()() override;
};

extern template class RNGThreadSafe<FakeRandom<OHMMS_PRECISION_FULL>>;
extern template class RNGThreadSafe<StdRandom<OHMMS_PRECISION_FULL>>;

using RandomGenerator = StdRandom<OHMMS_PRECISION_FULL>;
extern RNGThreadSafe<RandomGenerator> random_global;
#define Random random_global
} // namespace qmcplusplus

#endif
