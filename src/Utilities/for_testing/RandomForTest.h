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

#include <vector>
#include "Utilities/StdRandom.h"

namespace qmcplusplus
{
namespace testing
{
template<typename REAL>
class RandomForTest
{
public:
  RandomForTest();
  std::vector<REAL> getRealRandoms(int ncount);
  void makeRngReals(std::vector<REAL>& rng_reals);

private:
  StdRandom<REAL> rng;
};

extern template class RandomForTest<double>;
extern template class RandomForTest<float>;
} // namespace testing
} // namespace qmcplusplus
#endif
