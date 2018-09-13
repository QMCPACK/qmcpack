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

#include "QMCWaveFunctions/BsplineFactory/SplineR2RAdoptor.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"

#include "QMCWaveFunctions/BsplineFactory/mocks/MockEinsplineSetBuilder.h"
#include <iostream>

#include "QMCWaveFunctions/Batching.h"

#include "Configuration.h"

namespace qmcplusplus
{

using QMCT = QMCTraits;
TEST_CASE("SplineAdoptorBatched_Instantiation", "[wavefunction]")
{
  SplineAdoptor<double, 3> testAdoptor;
}
  
TEST_CASE("SplineR2RAdoptor", "[wavefunction]")
{
  SplineR2RAdoptor<double, QMCT::RealType> testAdoptor;
}

TEST_CASE("SplineAdoptorReader<SplineR2RAdoptor>", "[wavefunction]")
{
  MockEinsplineSetBuilder* e = new MockEinsplineSetBuilder();
  BsplineReaderBase* aReader=nullptr;
  aReader = new SplineAdoptorReader<SplineR2RAdoptor<double, QMCT::RealType>,Batching::SINGLE>(dynamic_cast<SplineBuilder*>(e));
  
}

}
