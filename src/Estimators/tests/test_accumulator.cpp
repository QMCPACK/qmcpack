//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewin@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Message/catch_mpi_main.hpp"

#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Estimators/accumulators.h"


#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

template<typename T>
void test_real_accumulator_basic()
{
  //std::cout << "int eps = " << std::numeric_limits<T>::epsilon() << std::endl;
  //std::cout << "int max = " << std::numeric_limits<T::max() << std::endl;
  accumulator_set<T> a1;
  REQUIRE(a1.count() == 0);
  REQUIRE(a1.good() == false);
  REQUIRE(a1.bad() == true);
  REQUIRE(a1.mean() == Approx(0.0));
  REQUIRE(a1.mean2() == Approx(0.0));
  REQUIRE(a1.variance() == Approx(std::numeric_limits<T>::max()));

  a1(2.0);
  REQUIRE(a1.count() == 1);
  REQUIRE(a1.good() == true);
  REQUIRE(a1.bad() == false);
  REQUIRE(a1.result() == Approx(2.0));
  REQUIRE(a1.result2() == Approx(4.0));
  REQUIRE(a1.mean() == Approx(2.0));
  REQUIRE(a1.mean2() == Approx(4.0));
  REQUIRE(a1.variance() == Approx(0.0));

  a1.clear();
  REQUIRE(a1.count() == 0);
  REQUIRE(a1.result() == Approx(0.0));
  REQUIRE(a1.result2() == Approx(0.0));
}

TEST_CASE("accumulator basic float", "[estimators]")
{
  test_real_accumulator_basic<double>();
}

TEST_CASE("accumulator basic double", "[estimators]")
{
  test_real_accumulator_basic<double>();
}

template <typename T>
void test_real_accumulator()
{
  accumulator_set<T> a;
  a(2.0);
  a(3.1);
  a(-1.1);
  a(4.8);
  REQUIRE(a.count() == Approx(4.0));

  REQUIRE(a.result() == Approx(8.8));
  REQUIRE(a.result2() == Approx(37.86));
  REQUIRE(a.mean() == Approx(2.2));
  REQUIRE(a.mean2() == Approx(9.465));
  REQUIRE(a.variance() == Approx(4.625));

  std::pair<T, T> mv = a.mean_and_variance();
  REQUIRE(mv.first == Approx(a.mean()));
  REQUIRE(mv.second == Approx(a.variance()));

  REQUIRE(mean(a) == Approx(a.mean()));

  // check that this doesn't crash
  std::stringstream o;
  o << a;
}

TEST_CASE("accumulator some values float", "[estimators]")
{
  test_real_accumulator<float>();
}

TEST_CASE("accumulator some values double", "[estimators]")
{
  test_real_accumulator<double>();
}

template <typename T>
void test_real_accumulator_weights()
{
  accumulator_set<T> a;
  a(2.0, 1.0);
  a(3.1, 2.4);
  a(-1.1, 2.1);
  a(4.8, 3.3);
  REQUIRE(a.count() == Approx(8.8));

  REQUIRE(a.result() == Approx(22.97));
  REQUIRE(a.result2() == Approx(105.637));
  REQUIRE(a.mean() == Approx(2.6102272727));
  REQUIRE(a.mean2() == Approx(12.004204545));
  REQUIRE(a.variance() == Approx(5.1909181302));

  std::pair<T, T> mv = a.mean_and_variance();
  REQUIRE(mv.first == Approx(a.mean()));
  REQUIRE(mv.second == Approx(a.variance()));
}

TEST_CASE("accumulator with weights float", "[estimators]")
{
  test_real_accumulator_weights<float>();
}

TEST_CASE("accumulator with weights double", "[estimators]")
{
  test_real_accumulator_weights<double>();
}

}
