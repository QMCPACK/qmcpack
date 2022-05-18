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
  auto xml_test = [](const char* driver_xml) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(driver_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    QMCDriverInput qmcdriver_input;
    qmcdriver_input.readXML(node);
    REQUIRE(qmcdriver_input.get_qmc_method().size() > 0);
  };

  std::for_each(testing::valid_vmc_input_sections.begin() + testing::valid_vmc_input_vmc_batch_index,
                testing::valid_vmc_input_sections.end(), xml_test);

  std::for_each(testing::valid_dmc_input_sections.begin() + testing::valid_dmc_input_dmc_batch_index,
                testing::valid_dmc_input_sections.end(), xml_test);
}
} // namespace qmcplusplus
