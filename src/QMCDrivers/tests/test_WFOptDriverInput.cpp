//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "QMCDrivers/WFOpt/WFOptDriverInput.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
TEST_CASE("WFOptDriverInput Instantiation", "[drivers]") { WFOptDriverInput driver_input; }

TEST_CASE("WFOptDriverInput readXML", "[drivers]")
{
  auto xml_test = [](const char* driver_xml) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(driver_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    WFOptDriverInput wfoptdriver_input;
    wfoptdriver_input.readXML(node);

    if (wfoptdriver_input.get_opt_method() == "")
    {
      REQUIRE(wfoptdriver_input.get_opt_num_crowds() == 0);
      REQUIRE(wfoptdriver_input.get_opt_crowd_size() == 1);
    }
    else if (wfoptdriver_input.get_opt_method() == "test")
    {
      REQUIRE(wfoptdriver_input.get_opt_num_crowds() == 4);
      REQUIRE(wfoptdriver_input.get_opt_crowd_size() == 8);
    }
    else
    {
      std::cout << "Unknown opt method: " << wfoptdriver_input.get_opt_method() << std::endl;
      REQUIRE(false); // optimizer method name not one of the two options
    }
  };

  std::for_each(testing::valid_opt_input_sections.begin(), testing::valid_opt_input_sections.end(), xml_test);
}
} // namespace qmcplusplus
