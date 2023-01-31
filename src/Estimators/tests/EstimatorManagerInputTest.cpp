//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File refactored from: Refactored from test_manager.cpp
//////////////////////////////////////////////////////////////////////////////////////
#include "EstimatorManagerInputTest.h"

#include "catch.hpp"

#include "ValidOneBodyDensityMatricesInput.h"
#include "ValidSpinDensityInput.h"
#include "ValidMomentumDistributionInput.h"
#include "ValidScalarEstimatorInput.h"

namespace qmcplusplus
{
namespace testing
{

Libxml2Document createEstimatorManagerNewGlobalInputXML()
{
  const int max_node_recurse = 3;
  Libxml2Document estimators_doc;
  estimators_doc.newDoc("Estimators");
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }

  return estimators_doc;
}

Libxml2Document createEstimatorManagerNewInputXML()
{
  const int max_node_recurse = 3;
  Libxml2Document estimators_doc;
  estimators_doc.newDoc("Estimators");
  {
    using namespace testing::onebodydensitymatrices;
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_momentum_distribution_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }
  for (auto& input_xml : valid_scalar_estimator_input_sections)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }

  return estimators_doc;
}

Libxml2Document createEstimatorManagerNewVMCInputXML()
{
  const int max_node_recurse = 3;
  Libxml2Document estimators_doc;
  estimators_doc.newDoc("Estimators");
  {
    using namespace testing::onebodydensitymatrices;
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_momentum_distribution_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_scalar_estimator_input_sections[local_energy_input]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  }

  return estimators_doc;
}

  
} // namespace testing
} // namespace qmcplusplus
