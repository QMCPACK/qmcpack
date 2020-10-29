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
#include "type_traits/CUDATypes.h"


namespace qmcplusplus
{
template<typename P, typename V>
class TestDeviceCUDA
{
  using CTS = CUDATypes<P, V, 3>;
  typename CTS::RealType testReal;
  typename CTS::ComplexType testComplex;
};

TEST_CASE("CUDA_Type_Aliases", "[type_traits][CUDA]")
{
  TestDeviceCUDA<float, float> float_test;
  TestDeviceCUDA<double, double> double_test;
  TestDeviceCUDA<float, std::complex<float>> complex_float_test;
  TestDeviceCUDA<double, std::complex<double>> complex_double_test;

  // This should cause compiler error
  // Realtype and ValueType precision do not match
  //TestDeviceCUDA<float, double> pv_mismatch_test;
  //REQUIRE(floatTest);
}

} // namespace qmcplusplus
