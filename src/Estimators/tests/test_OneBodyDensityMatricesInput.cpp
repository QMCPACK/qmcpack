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

#include "OneBodyDensityMatricesInput.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <iostream>

namespace qmcplusplus
{
TEST_CASE("OneBodyDensityMatricesInput::from_xml", "[estimators]")
{
  using POLT    = PtclOnLatticeTraits;
  using Lattice = POLT::ParticleLayout;
  using namespace testing::onebodydensitymatrices;
  for (auto input_xml : valid_one_body_density_matrices_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);
  }

  for (auto input_xml : testing::invalid_one_body_density_matrices_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    CHECK_THROWS_AS(OneBodyDensityMatricesInput(node), UniformCommunicateError);
  }
}

} // namespace qmcplusplus
