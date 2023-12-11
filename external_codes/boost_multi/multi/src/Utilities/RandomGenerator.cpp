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
typename RNG::result_type RNGThreadSafe<RNG>::operator()()
{
  result_type result;
#pragma omp critical
  {
    result = RNG::operator()();
  }
  return result;
}

template class RNGThreadSafe<FakeRandom<OHMMS_PRECISION_FULL>>;
template class RNGThreadSafe<RandomGenerator>;

RNGThreadSafe<FakeRandom<OHMMS_PRECISION_FULL>> fake_random_global;
RNGThreadSafe<RandomGenerator> random_global;
} // namespace qmcplusplus
