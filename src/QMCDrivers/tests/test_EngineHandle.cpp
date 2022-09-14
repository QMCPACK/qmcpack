//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu, University of California Berkeley
//
// File created by: Leon Otis, leon_otis@berkeley.edu, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "QMCDrivers/WFOpt/EngineHandle.h"
#include "Configuration.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
using ValueType         = qmcplusplus::QMCTraits::ValueType;

///This provides a basic test of constructing an EngineHandle object and checking information in it
TEST_CASE("EngineHandle construction", "[drivers]")
{
  Communicate* c = OHMMS::Controller;


  const std::string engine_input("<tmp> </tmp>");

  Libxml2Document doc;
  bool okay = doc.parseFromString(engine_input);
  REQUIRE(okay);

  xmlNodePtr fakeXML = doc.getRoot();

  DescentEngine descentEngineObj = DescentEngine(c, fakeXML);

  descentEngineObj.processXML(fakeXML);

  app_log() << "Test of DescentEngineHandle construction" << std::endl;
  std::unique_ptr<DescentEngineHandle> handle = std::make_unique<DescentEngineHandle>(descentEngineObj);


  const int fake_num_params = 5;
  const int fake_sample_num = 100;
  handle->prepareSampling(fake_num_params, fake_sample_num);
  auto& test_der_rat_samp = handle->getVector();

  REQUIRE(test_der_rat_samp.size() == 6);
  REQUIRE(test_der_rat_samp[0] == 0.0);
}
} // namespace qmcplusplus
