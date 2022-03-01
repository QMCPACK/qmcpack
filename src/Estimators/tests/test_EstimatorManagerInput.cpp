//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "EstimatorManagerInput.h"
#include "test_EstimatorManagerInput.h"
#include "SpinDensityInput.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"

namespace qmcplusplus
{

namespace testing
{
class EstimatorManagerInputTests
{
public:
  /** @ingroup testing "private" methods of estimator class
   *  useful for this new implementation
   *  @{
   */
  /// The simplest insertion of a input class
  void testAppendMinimal(EstimatorManagerInput& emi) { emi.appendEstimatorInput<OneBodyDensityMatricesInput>(); }
  /// Actually from valid xml.
  template<class T>
  void testAppendFromXML(EstimatorManagerInput& emi, xmlNodePtr node)
  {
    emi.appendEstimatorInput<T>(node);
  }
  /** @} */
};

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
    bool okay = doc.parseFromString(valid_spin_density_input_sections[0]);
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
  return estimators_doc;
}

} // namespace testing

TEST_CASE("EstimatorManagerInput::testInserts", "[estimators]")
{
  using namespace testing;
  EstimatorManagerInputTests emit;
  EstimatorManagerInput emi;
  emit.testAppendMinimal(emi);

  {
    using namespace testing::onebodydensitymatrices;
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_one_body_density_matrices_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<OneBodyDensityMatricesInput>(emi, node);
  }
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<SpinDensityInput>(emi, node);
  }
}

TEST_CASE("EstimatorManagerInput::readXML", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());

  CHECK(emi.get_estimator_inputs().size() == 3);

  // CHECK EMI throws if unparsable estimators are in input.
  Libxml2Document doc;
  std::string bad_estimator = R"XML(
<estimator type="NeutrinoDensity" name="bad_estimator"/>
)XML";
  bool okay                 = doc.parseFromString(bad_estimator);
  REQUIRE(okay);
  xmlNodePtr node      = doc.getRoot();
  int max_node_recurse = 1;
  estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  CHECK_THROWS_AS(EstimatorManagerInput(estimators_doc.getRoot()), UniformCommunicateError);
}

} // namespace qmcplusplus
