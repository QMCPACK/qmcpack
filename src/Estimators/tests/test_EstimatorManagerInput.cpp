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

#include "EstimatorManagerInput.h"
#include "SpinDensityInput.h"
#include "MomentumDistributionInput.h"
#include "OneBodyDensityMatricesInput.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "ValidSpinDensityInput.h"
#include "OhmmsData/Libxml2Doc.h"

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
  void testAppendMinimal(EstimatorManagerInput& emi)
  {
    emi.appendEstimatorInput<OneBodyDensityMatricesInput>();
  }
  /// Actually from valid xml.
  template<class T>
  void testAppendFromXML(EstimatorManagerInput& emi, xmlNodePtr node)
  {
    emi.appendEstimatorInput<T>(node);
  }
  /** @} */
};
}
  
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
}
