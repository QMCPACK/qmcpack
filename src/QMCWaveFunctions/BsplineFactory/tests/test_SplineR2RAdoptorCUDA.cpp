//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "QMCWaveFunctions/BsplineFactory/BsplineDeviceCUDA.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorVectorized.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2RAdoptorCUDA.h"
#include <iostream>

namespace qmcplusplus
{

TEST_CASE("SplineAdoptorVectorized_Instantiation", "[wavefunction]")
{
  SplineAdoptorVectorized<BsplineDeviceCUDA, double, 3> testAdoptor;
}
  
TEST_CASE("SplineR2RAdoptorCUDA_Instantiation", "[wavefunction]")
{
  SplineR2RAdoptorCUDA<double, double> testAdoptor;
}

}

