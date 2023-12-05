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

#ifndef QMCPLUSPLUS_RANDOMFORTEST_H
#define QMCPLUSPLUS_RANDOMFORTEST_H

#include <complex>
#include <vector>
#include "Utilities/StdRandom.h"

namespace qmcplusplus
{
namespace testing
{

/** Get a known sequence of random numbers for testing.
 *  VT is the floating point precision 
 *  While inelegant to have separate named calls for the cplx types in the same class
 *  separate class templates for RandomForTest<double> and RandomForTest<std::complex<double>>
 *  turned out to be surprisingly difficult. Someone is welcome to try when we required > c++14
 */
template<typename VT>
class RandomForTest
{
public:
  RandomForTest();
  std::vector<VT> getRngVec(int ncount);
  std::vector<std::complex<VT>> getRngVecComplex(int ncount);
  void fillVecRng(std::vector<VT>& rng_reals);
  void fillVecRng(std::vector<std::complex<VT>>& rng_reals);
  void fillBufferRng(VT* rng_reals, size_t number);
  void fillBufferRng(std::complex<VT>* rng_reals, size_t number);
  VT operator()();
private:
  StdRandom<VT> rng;
};
  
  extern template class RandomForTest<double>;
  extern template class RandomForTest<float>;
} // namespace testing
} // namespace qmcplusplus
#endif
