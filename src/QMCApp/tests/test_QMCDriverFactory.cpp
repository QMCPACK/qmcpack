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

#include "Utilities/RandomGenerator.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCApp/QMCDriverFactory.h"


#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>

#include "Message/Communicate.h"

namespace qmcplusplus
{

TEST_CASE("QMCDriverFactory::VMCBatchedDriver", "[qmcapp]")
{
  Communicate comm;
  QMCDriverFactory driver_factory(&comm);
  // clang-format off
  const char* driver_xml = R"(
  <qmc method="vmc_batch" move="pbyp">
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
  xmlNodePtr node = doc.getRoot();
  QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(0, node);
  REQUIRE(das.new_run_type == QMCDriverFactory::QMCRunType::VMC_BATCH);
  REQUIRE(das.append_run == false);
}
}
