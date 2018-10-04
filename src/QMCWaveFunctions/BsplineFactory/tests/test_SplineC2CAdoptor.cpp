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

#include "QMCWaveFunctions/BsplineFactory/SplineC2CAdoptor.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"

#include "QMCWaveFunctions/BsplineFactory/mocks/MockEinsplineSetBuilder.h"
#include <iostream>

#include "Batching.h"

#include "Configuration.h"

namespace qmcplusplus
{

using QMCT = QMCTraits;
TEST_CASE("SplineAdoptorBatched_Instantiation", "[wavefunction]")
{
  SplineAdoptor<double, 3> testAdoptor;
}
  
TEST_CASE("SplineR2RAdoptorBatched_Instantiation", "[wavefunction]")
{
  SplineC2CAdoptor<double, QMCT::RealType> testAdoptor;
}

TEST_CASE("SplineAdoptorReader<SplineC2CAdoptorBatched>", "[wavefunction]")
{
  MockEinsplineSetBuilder* e = new MockEinsplineSetBuilder();
  BsplineReaderBase* aReader=nullptr;
  aReader = new SplineAdoptorReader<SplineC2CAdoptor<double, QMCT::RealType>,Batching::SINGLE>(dynamic_cast<SplineBuilder*>(e));
  
}

}
