//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "RandomGenerator.h"
#include <ctime>

uint32_t make_seed(int i, int n) { return static_cast<uint32_t>(std::time(0)) % 10474949 + (i + 1) * n + i; }

namespace qmcplusplus
{

template<class RNG>
typename RNG::result_type RNGThreadSafe<RNG>::rand()
{
  result_type result;
// This should be a named section but at least clang 9 doesn't seem to support
// and warns of extra tokens.
#pragma omp critical
  {
    result = RNG::rand();
  }
  return result;
}

template<class RNG>
typename RNG::result_type RNGThreadSafe<RNG>::operator()()
{
  result_type result;
#pragma omp critical
  {
    result = RNG::rand();
  }
  return result;
}

template class RNGThreadSafe<FakeRandom>;
template class RNGThreadSafe<StdRandom<float>>;
template class RNGThreadSafe<StdRandom<double>>;

RNGThreadSafe<StdRandom<OHMMS_PRECISION_FULL>> boost_random_global;
RNGThreadSafe<FakeRandom> fake_random_global;
} // namespace qmcplusplus
