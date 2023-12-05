//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "RandomForTest.h"
#include <algorithm>

namespace qmcplusplus
{
namespace testing
{
template<typename VT>
RandomForTest<VT>::RandomForTest()
{
  rng.init(111);
}

template<typename VT>
std::vector<VT> RandomForTest<VT>::getRngVec(int ncount)
{
  std::vector<VT> rng_reals;
  rng_reals.reserve(ncount);
  std::generate_n(std::back_inserter(rng_reals), ncount, rng);
  return rng_reals;
}

template<typename VT>
std::vector<std::complex<VT>> RandomForTest<VT>::getRngVecComplex(int ncount)
{
  std::vector<std::complex<VT>> rngs_cplx(ncount);
  for( auto& rng_cplx : rngs_cplx)
    rng_cplx = {rng(), rng()};
  return rngs_cplx;
}


template<typename VT>
void RandomForTest<VT>::fillVecRng(std::vector<VT>& rngReals)
{
  // until c++ std = 17
  //std::generate(rng_reals.begin(), rng_reals.end(), rng());
  for (auto& rng_real : rngReals)
    rng_real = rng();
}

template<typename VT>
void RandomForTest<VT>::fillVecRng(std::vector<std::complex<VT>>& cplx_nums)
{
  // until c++ std = 17
  //std::generate(rng_reals.begin(), rng_reals.end(), rng());
  for (auto& cplx_num : cplx_nums)
    cplx_num = std::complex<VT>{rng(), rng()};
}

template<typename VT>
void RandomForTest<VT>::fillBufferRng(VT* rng_reals, size_t count)
{
  for (size_t ir = 0; ir < count; ++ir)
  {
    *rng_reals = rng();
    ++rng_reals;
  }
}

template<typename VT>
void RandomForTest<VT>::fillBufferRng(std::complex<VT>* cplx_nums, size_t count)
{
  for (size_t i = 0; i < count; ++i)
  {
    cplx_nums->real(rng());
    cplx_nums->imag(rng());
    ++cplx_nums;
  }
}

template<typename VT>
VT RandomForTest<VT>::operator()()
{
  return rng();
}

template class RandomForTest<double>;
template class RandomForTest<float>;
} // namespace testing
} // namespace qmcplusplus
