//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"

#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("SpinDensityInput::readXML", "[estimators]")
{
  auto xml_test = [](const char* input_xml) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    SpinDensityInput sdi;
    sdi.readXML(node);
    CHECK(sdi.get_npoints() == 1000);
    TinyVector<int, SpinDensityInput::DIM> grid(10, 10, 10);
    CHECK(sdi.get_grid() == grid);
    TinyVector<int, SpinDensityInput::DIM> gdims(100, 10, 1);
    CHECK(sdi.get_gdims() == gdims);
  };
  std::for_each(testing::valid_spin_density_input_sections.begin(),testing::valid_spin_density_input_sections.end(), xml_test);
}

} // namespace qmcplusplus
