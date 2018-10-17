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

#include "catch.hpp"
#include "CUDA/CudaTypeAliases.h"

template<typename P>
class TestDeviceCUDA
{
  using CTA = CudaTypeAliases<P>;
  CTA::RealType testReal;
  CTA::ComplexType testComplex;
}

TEST_CASE("CUDA_Type_Aliases", "[CUDA]")
{
  TestDeviceCUDA<float> floatTest;
  TestDeviceCUDA<double> doubleTest;
}
