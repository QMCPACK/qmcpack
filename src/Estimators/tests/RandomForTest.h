//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RANDOMFORTEST_H
#define QMCPLUSPLUS_RANDOMFORTEST_H

#include <vector>
#include "Configuration.h"
#include "Utilities/StdRandom.h"

namespace qmcplusplus
{
namespace testing
{
class RandomForTest
{
public:
  using Real = QMCTraits::RealType;
  RandomForTest();
  std::vector<Real> getRealRandoms(int ncount);
  void makeRngReals(std::vector<Real>& rng_reals);

private:
  StdRandom<Real> rng;
};
} // namespace testing
} // namespace qmcplusplus
#endif
