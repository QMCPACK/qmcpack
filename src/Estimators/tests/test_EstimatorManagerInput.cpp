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

#include "catch.hpp"
#include "EstimatorManagerInput.h"
#include "EstimatorManagerInputTest.h"
#include "EstimatorInputDelegates.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "ValidSpinDensityInput.h"

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
    Input input;
    Libxml2Document doc;
    bool okay = doc.parseFromString(input[Input::valid::VANILLA]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<OneBodyDensityMatricesInput>(emi, node);
  }
  {
    Libxml2Document doc;
    using spin_input = testing::ValidSpinDensityInput;
    bool okay = doc.parseFromString(spin_input::xml[spin_input::GRID]);
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
  CHECK(emi.get_scalar_estimator_inputs().size() == 5);

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
    Input input;
    Libxml2Document doc;
    bool okay = doc.parseFromString(input[Input::valid::VANILLA]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    emit.testAppendFromXML<OneBodyDensityMatricesInput>(emi, node);
  }
  {
    Libxml2Document doc;
    using spin_input = testing::ValidSpinDensityInput;
    bool okay = doc.parseFromString(spin_input::xml[spin_input::GRID]);
    REQUIRE(okay);
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

  CHECK(emi.get_estimator_inputs().size() == 3);
  CHECK(emi.get_scalar_estimator_inputs().size() == 5);

  EstimatorManagerInput emi_moved_to(std::move(emi));

  CHECK(emi_moved_to.get_estimator_inputs().size() == 3);
  CHECK(emi_moved_to.get_scalar_estimator_inputs().size() == 5);
}

TEST_CASE("EstimatorManagerInput::MergeConstructor", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc        = createEstimatorManagerNewInputXML();
  Libxml2Document global_estimators_doc = createEstimatorManagerNewGlobalInputXML();
  EstimatorManagerInput emi_global(global_estimators_doc.getRoot());
  EstimatorManagerInput emi_local(estimators_doc.getRoot());
  EstimatorManagerInput emi_merged{emi_global, emi_local};

  CHECK(emi_merged.get_estimator_inputs().size() == 4);
  CHECK(emi_merged.get_scalar_estimator_inputs().size() == 5);
}

TEST_CASE("EstimatorManagerInput::AddAnInputDirectly", "[estimators]")
{
  using namespace testing;
  Libxml2Document estimators_doc        = createEstimatorManagerNewInputXML();
  Libxml2Document global_estimators_doc = createEstimatorManagerNewGlobalInputXML();
  EstimatorManagerInput emi_global(global_estimators_doc.getRoot());
  EstimatorManagerInput emi_local(estimators_doc.getRoot());
  EstimatorManagerInput emi_merged{emi_global, emi_local};

  CHECK(emi_merged.get_estimator_inputs().size() == 4);
  CHECK(emi_merged.get_scalar_estimator_inputs().size() == 5);
}

} // namespace qmcplusplus
