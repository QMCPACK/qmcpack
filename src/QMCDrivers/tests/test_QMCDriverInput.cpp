//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "QMCDrivers/QMCDriverInput.h"
#include "EstimatorInputDelegates.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
TEST_CASE("QMCDriverInput Instantiation", "[drivers]") { QMCDriverInput driver_input; }

TEST_CASE("QMCDriverInput readXML", "[drivers]")
{
  auto xml_test = [](std::string_view driver_xml) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(driver_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(node);
    REQUIRE(qmcdriver_input.get_qmc_method().size() > 0);
  };

  using VMCInputs = testing::VmcInputs;
  std::for_each(VMCInputs::begin(), VMCInputs::end(), xml_test);
  using DMCInputs = testing::DmcInputs;
  std::for_each(DMCInputs::begin(), DMCInputs::end(), xml_test);
}
} // namespace qmcplusplus
