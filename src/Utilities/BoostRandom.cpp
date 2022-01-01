//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "BoostRandom.h"
#include <cmath>

uint32_t make_seed(int i, int n);

template<typename T, typename RNG>
void BoostRandom<T, RNG>::init(int iseed_in)
{
  uint_type baseSeed = iseed_in;
  uni.engine().seed(baseSeed);
}

template class BoostRandom<float>;
template class BoostRandom<double>;
