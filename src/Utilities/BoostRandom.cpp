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

struct BoxMuller2
{
  template<typename RNG>
  static inline void generate(RNG& rng, double* restrict a, int n)
  {
    for (int i = 0; i + 1 < n; i += 2)
    {
      double temp1 = 1.0 - 0.9999999999 * rng(), temp2 = rng();
      a[i]     = sqrt(-2.0 * log(temp1)) * cos(6.283185306 * temp2);
      a[i + 1] = sqrt(-2.0 * log(temp1)) * sin(6.283185306 * temp2);
    }
    if (n % 2 == 1)
    {
      double temp1 = 1 - 0.9999999999 * rng(), temp2 = rng();
      a[n - 1] = sqrt(-2.0 * log(temp1)) * cos(6.283185306 * temp2);
    }
  }

  template<typename RNG>
  static inline void generate(RNG& rng, float* restrict a, int n)
  {
    for (int i = 0; i + 1 < n; i += 2)
    {
      float temp1 = 1.0f - 0.9999999999f * rng(), temp2 = rng();
      a[i]     = sqrtf(-2.0f * logf(temp1)) * cosf(6.283185306f * temp2);
      a[i + 1] = sqrtf(-2.0f * logf(temp1)) * sinf(6.283185306f * temp2);
    }
    if (n % 2 == 1)
    {
      float temp1 = 1.0f - 0.9999999999f * rng(), temp2 = rng();
      a[n - 1] = sqrtf(-2.0f * logf(temp1)) * cosf(6.283185306f * temp2);
    }
  }
};

template<typename T, typename RNG>
void BoostRandom<T, RNG>::init(int i, int nstr, int iseed_in, uint_type offset)
{
  uint_type baseSeed = iseed_in;
  myContext          = i;
  nContexts          = nstr;
  if (iseed_in <= 0)
    baseSeed = make_seed(i, nstr);
  baseOffset = offset;
  uni.engine().seed(baseSeed);
}

template<typename T, typename RNG>
void BoostRandom<T, RNG>::generate_normal(T* restrict d, int n)
{
  BoxMuller2::generate(*this, d, n);
}

template class BoostRandom<float>;
template class BoostRandom<double>;
