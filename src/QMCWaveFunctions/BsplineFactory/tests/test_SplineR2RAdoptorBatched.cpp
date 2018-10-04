//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "QMCWaveFunctions/BsplineFactory/BsplineDeviceCUDA.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorBatched.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2RAdoptorBatched.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"

#include "QMCWaveFunctions/BsplineFactory/mocks/MockEinsplineSetBuilder.h"
#include <iostream>

#include "Batching.h"

namespace qmcplusplus
{

TEST_CASE("SplineAdoptorBatched_Instantiation", "[wavefunction]")
{
  SplineAdoptorBatched<BsplineDeviceCUDA, double, 3> testAdoptor;
}
  
TEST_CASE("SplineR2RAdoptorBatched_Instantiation", "[wavefunction]")
{
  SplineR2RAdoptorBatched<double, double> testAdoptor;
}

TEST_CASE("SplineAdoptorReader<SplineR2RAdoptorBatched>", "[wavefunction]")
{
  MockEinsplineSetBuilder* e = new MockEinsplineSetBuilder();
  SplineAdoptorReader<SplineR2RAdoptorBatched<double, double>, Batching::BATCHED>(dynamic_cast<SplineBuilder*>(e));
  
}

}
