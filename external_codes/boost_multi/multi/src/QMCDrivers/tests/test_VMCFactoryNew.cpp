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

#include "OhmmsData/Libxml2Doc.h"

#include "QMCDrivers/VMC/VMCFactoryNew.h"
#include "QMCDrivers/DriverTraits.h"
#include "QMCDrivers/tests/ValidQMCInputSections.h"

namespace qmcplusplus
{
TEST_CASE("VMCFactory Instantiation", "[drivers]")
{
  using namespace testing;
  Libxml2Document doc;

  bool okay = doc.parseFromString(valid_vmc_input_sections[valid_vmc_input_vmc_batch_index]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  std::bitset<QMC_MODE_MAX> vmc_mode;
  vmc_mode[UPDATE_MODE] = true;
  int qmc_counter       = 0;
  VMCFactoryNew vmc_factory(node, vmc_mode[UPDATE_MODE]);
}

} // namespace qmcplusplus
