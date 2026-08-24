//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"
#include "EstimatorTesting.h"
#include "Message/UniformCommunicateError.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("SpinDensityInput::readXML", "[estimators]")
{
  using input = testing::ValidSpinDensityInput;
  for (auto input_xml : input::xml)
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input_xml));
    xmlNodePtr node = doc.getRoot();

    SpinDensityInput sdi(node);

    Lattice lattice;
    if (sdi.get_cell().explicitly_defined == true)
      lattice = sdi.get_cell();
    else
      lattice = testing::makeTestLattice();

    SpinDensityInput::DerivedParameters dev_par = sdi.calculateDerivedParameters(lattice);

    CHECK(dev_par.npoints == 1000);
    TinyVector<int, SpinDensityInput::DIM> grid(10, 10, 10);
    CHECK(dev_par.grid == grid);
    TinyVector<int, SpinDensityInput::DIM> gdims(100, 10, 1);
    CHECK(dev_par.gdims == gdims);
  }
}

TEST_CASE("SpinDensityInput invalid input", "[estimators]")
{
  constexpr std::array<std::string_view, 4> invalid_xml{
      R"XML(<estimator type="spindensity"><parameter name="dr">1 1 1</parameter><parameter name="grid">1 1 1</parameter></estimator>)XML",
      R"XML(<estimator type="spindensity"><parameter name="grid">1 1 1</parameter><parameter name="corner">0 0 0</parameter><parameter name="center">0 0 0</parameter></estimator>)XML",
      R"XML(<estimator type="spindensity"><parameter name="grid">1 1 1</parameter><parameter name="cell">1 0 0 0 1 0 0 0 1</parameter></estimator>)XML",
      R"XML(<estimator type="spindensity"><parameter name="grid">1 1 1</parameter><parameter name="center">0 0 0</parameter><parameter name="cell">1 0 0 0 1 0 0 0</parameter></estimator>)XML"};

  for (const auto input_xml : invalid_xml)
  {
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(input_xml));
    CHECK_THROWS_AS(SpinDensityInput(doc.getRoot()), UniformCommunicateError);
  }
}

} // namespace qmcplusplus
