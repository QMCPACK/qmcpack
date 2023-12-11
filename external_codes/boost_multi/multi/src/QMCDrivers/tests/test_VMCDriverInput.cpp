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

#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
TEST_CASE("VMCDriverInput readXML", "[drivers]")
{
  auto xml_test = [](const char* driver_xml) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(driver_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    VMCDriverInput vmcdriver_input;
    vmcdriver_input.readXML(node);
    REQUIRE(vmcdriver_input.get_use_drift() == false);
  };

  std::for_each(testing::valid_vmc_input_sections.begin() + testing::valid_vmc_input_vmc_batch_index,
                testing::valid_vmc_input_sections.end(), xml_test);
}


} // namespace qmcplusplus
