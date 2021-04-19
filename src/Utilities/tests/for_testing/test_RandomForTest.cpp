//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "RandomForTest.cpp"
#include "../StlPrettyPrint.hpp"
#include <iostream>
#include <string>

/** \file
 */
namespace qmcplusplus
{
// When linking containers does cause asan failure this should return
TEST_CASE("RandomForTest_generateRngReals", "[utilities][for_testing]")
{
  testing::RandomForTest<double> rng_for_test;
  std::vector<double> rng_vector(5, 0.0);
  rng_for_test.generateRngReals(rng_vector.data(), 5);
  std::vector<double> reference{0.120758, 0.301789, 0.853906, 0.297716, 0.377862};
  auto checkArray = [](auto A, auto B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == Approx(B[i]));
    }
  };
  checkArray(rng_vector.data(), reference.data(),5);
}

TEST_CASE("RandomForTest_makeRngReals", "[utilities][for_testing]")
{
  testing::RandomForTest<double> rng_for_test;
  std::vector<double> rng_vector(5);
  rng_for_test.makeRngReals(rng_vector);
  std::vector<double> reference{0.120758, 0.301789, 0.853906, 0.297716, 0.377862};
  auto checkArray = [](auto A, auto B, int n) {
    for (int i = 0; i < n; ++i)
    {
      CHECK(A[i] == Approx(B[i]));
    }
  };
  checkArray(rng_vector.data(), reference.data(),5);
}

TEST_CASE("RandomForTest_call_operator", "[utilities][for_testing]")
{
  testing::RandomForTest<double> rng_for_test;
  double test = rng_for_test();
  CHECK(test == Approx(0.120758));
}


} // namespace qmcplusplus
