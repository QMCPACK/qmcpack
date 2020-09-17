//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include <complex>
#include "catch.hpp"
#include "type_traits/QMCTypes.h"


namespace qmcplusplus
{
template<typename P>
class TestQMCTypes
{
  using TA = QMCTypes<P, 3>;
  typename TA::RealType testReal;
  typename TA::ComplexType testComplex;
};

TEST_CASE("QMCTypes", "[type_traits]")
{
  TestQMCTypes<float> float_test;
  TestQMCTypes<double> double_test;

  // This should cause compiler error
  // Realtype and ValueType precision do not match
  //TestDeviceCUDA<float, double> pv_mismatch_test;
  //REQUIRE(floatTest);
}

} // namespace qmcplusplus
