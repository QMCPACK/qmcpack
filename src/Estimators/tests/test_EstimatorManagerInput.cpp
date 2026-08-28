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
#include "EstimatorManagerInput.h"
#include "EstimatorManagerInputTest.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "ValidSpinDensityInput.h"
#include "ValidStructureFactorInput.h"
#include "test_EstimatorManagerInput.h"
#include "PairCorrelationInput.h"

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
  /// Actually from valid xml.
  template<class T>
  void testAppendFromXML(EstimatorManagerInput& emi, xmlNodePtr node)
  {
    emi.appendEstimatorInput<T>(node);
  }
  /** @} */
};

} // namespace testing

TEST_CASE("EstimatorManagerInput::testInserts", "[estimators]")
{
  using namespace testing;
  EstimatorManagerInputTests emit;
  EstimatorManagerInput emi;

  {
    using Input = testing::ValidOneBodyDensityMatricesInput;
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(Input::getXml(Input::valid::VANILLA)));
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<OneBodyDensityMatricesInput>(emi, node);
  }
  {
    Libxml2Document doc;
    using spin_input = testing::SpinDensityInputs;
    REQUIRE(doc.parseFromString(spin_input::getXml(spin_input::valid::GRID)));
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<SpinDensityInput>(emi, node);
  }
}

TEST_CASE("EstimatorManagerInput::readXML", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());

  CHECK(emi.get_estimator_inputs().size() == n_opest_new_input_xml);
  CHECK(emi.get_scalar_estimator_inputs().size() == 4);

  // CHECK EMI throws if unparsable estimators are in input.
  Libxml2Document doc;
  std::string bad_estimator = R"XML(
<estimator type="NeutrinoDensity" name="bad_estimator"/>
)XML";
  REQUIRE(doc.parseFromString(bad_estimator));
  xmlNodePtr node      = doc.getRoot();
  int max_node_recurse = 1;
  estimators_doc.addChild(xmlCopyNode(node, max_node_recurse));
  CHECK_THROWS_AS(EstimatorManagerInput(estimators_doc.getRoot()), UniformCommunicateError);
}

template<class INPUT>
class TakesAMovedInput
{
public:
  TakesAMovedInput(INPUT&& obdmi) : input_(obdmi) {}

private:
  const INPUT input_;
};

TEST_CASE("EstimatorManagerInput::moveFromEstimatorInputs", "[estimators]")
{
  using namespace testing;
  EstimatorManagerInputTests emit;
  EstimatorManagerInput emi;

  {
    using Input = testing::ValidOneBodyDensityMatricesInput;
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(Input::getXml(Input::valid::VANILLA)));
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<OneBodyDensityMatricesInput>(emi, node);
  }
  {
    Libxml2Document doc;
    using spin_input = testing::SpinDensityInputs;
    REQUIRE(doc.parseFromString(spin_input::getXml(spin_input::valid::GRID)));
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<SpinDensityInput>(emi, node);
  }
  TakesAMovedInput<SpinDensityInput> takes_sdi(
      std::move(std::get<SpinDensityInput>(emi.get_estimator_inputs().back())));
  emi.get_estimator_inputs().pop_back();
  TakesAMovedInput<OneBodyDensityMatricesInput> takes_obdmi(
      std::move(std::get<OneBodyDensityMatricesInput>(emi.get_estimator_inputs().back())));
  emi.get_estimator_inputs().pop_back();
}

TEST_CASE("EstimatorManagerInput::moveConstructor", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi(estimators_doc.getRoot());

  CHECK(emi.get_estimator_inputs().size() == n_opest_new_input_xml);
  CHECK(emi.get_scalar_estimator_inputs().size() == 4);

  EstimatorManagerInput emi_moved_to(std::move(emi));

  CHECK(emi_moved_to.get_estimator_inputs().size() == n_opest_new_input_xml);
  CHECK(emi_moved_to.get_scalar_estimator_inputs().size() == 4);
}

TEST_CASE("EstimatorManagerInput::MergeConstructor", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc        = createEstimatorManagerNewInputXML();
  Libxml2Document global_estimators_doc = createEstimatorManagerNewGlobalInputXML();
  EstimatorManagerInput emi_global(global_estimators_doc.getRoot());
  EstimatorManagerInput emi_local(estimators_doc.getRoot());
  EstimatorManagerInput emi_merged{emi_global, emi_local};

  CHECK(emi_merged.get_estimator_inputs().size() == n_opest_new_input_xml);
  CHECK(emi_merged.get_scalar_estimator_inputs().size() == 5);
}

TEST_CASE("EstimatorManagerInput::Name")
{
  // This test case covers name initialization when it is missing from input.
  // It also covers breakage in the canned estimator manager input mocking.
  using namespace testing;
  Libxml2Document estimators_doc = createEstimatorManagerNewInputXML();
  EstimatorManagerInput emi_local(estimators_doc.getRoot());

  auto checkName = [](EstimatorManagerInput& emi, auto ent) {
    using EstimatorInput = typename decltype(ent)::Type;
    auto indexes         = emi.getEstimatorTypeIndexes<EstimatorInput>();
    for (auto index : indexes)
    {
      auto& eden_input = std::get<EstimatorInput>(emi.get_estimator_inputs()[index]);
      CHECK(eden_input.get_name() == ent.name);
    }
  };
  checkName(emi_local, ExpectedEstimatorInputName<EnergyDensityInput>{});
  checkName(emi_local, ExpectedEstimatorInputName<OneBodyDensityMatricesInput>{});
  checkName(emi_local, ExpectedEstimatorInputName<MomentumDistributionInput>{});
  checkName(emi_local, ExpectedEstimatorInputName<StructureFactorInput>{});
}

TEST_CASE("PairCorrelationInput retired debug attribute", "[estimators]")
{
  Libxml2Document valid_doc;
  REQUIRE(valid_doc.parseFromString(
      R"(<estimator type="PairCorrelation" name="pc" num_bin="8" rmax="4.0" dr="0.5"/>)"));
  PairCorrelationInput valid_input(valid_doc.getRoot());
  CHECK(valid_input.get_name() == "pc");
  CHECK(valid_input.get_nbins() == 8);
  CHECK(valid_input.get_rmax() == 4.0);
  CHECK(valid_input.get_delta() == 0.5);

  Libxml2Document retired_doc;
  REQUIRE(retired_doc.parseFromString(R"(<estimator type="PairCorrelation" debug="true"/>)"));
  CHECK_THROWS_AS(PairCorrelationInput(retired_doc.getRoot()), UniformCommunicateError);
}

} // namespace qmcplusplus
