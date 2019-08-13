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

namespace qmcplusplus
{
TEST_CASE("VMCFactory Instantiation", "[drivers]") {
  // clang-format off
  const char* driver_xml = R"(
  <qmc method="vmc" move="pbyp">
    <estimator name="LocalEnergy" hdf5="no" />
    <parameter name="walkers">                1 </parameter>
    <parameter name="stepsbetweensamples">    1 </parameter>
    <parameter name="warmupSteps">            5 </parameter>
    <parameter name="substeps">               5 </parameter>
    <parameter name="steps">                  1 </parameter>
    <parameter name="blocks">                 2 </parameter>
    <parameter name="timestep">             1.0 </parameter>
    <parameter name="usedrift">              no </parameter>
  </qmc>
)";
  // clang-format on

  Libxml2Document doc;
  bool okay = doc.parseFromString(driver_xml);
  REQUIRE(okay);
  xmlNodePtr node                           = doc.getRoot();
  std::bitset<QMC_MODE_MAX> vmc_mode;
  vmc_mode[SPACEWARP_MODE] = false;
  vmc_mode[MULTIPLE_MODE]  = false;
  vmc_mode[UPDATE_MODE] = true;
  vmc_mode[GPU_MODE] = false;
  int qmc_counter = 0;
  VMCFactoryNew vmc_factory(node,
                            vmc_mode.to_ulong(),
                            qmc_counter);
  
}

} // namespace qmcplusplus
