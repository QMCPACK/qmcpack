//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//
// File created by: Leon Otis, leon_otis@berkeley.edu University, University of California Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "QMCDrivers/WFOpt/EngineData.h"
#include "Configuration.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

using FullPrecValueType = qmcplusplus::QMCTraits::FullPrecValueType;
using ValueType         = qmcplusplus::QMCTraits::ValueType;

///This provides a basic test of constructing an EngineData struct and checking information in it
TEST_CASE("EngineData construction", "[drivers]")
{
  Communicate* c = OHMMS::Controller;


  const std::string engine_input("<tmp> </tmp>");

  Libxml2Document doc;
  bool okay = doc.parseFromString(engine_input);
  REQUIRE(okay);

  xmlNodePtr fakeXML = doc.getRoot();

  DescentEngine descentEngineObj = DescentEngine(c,fakeXML);
  DescentEngine* descentEnginePtr = &descentEngineObj;

  descentEnginePtr->processXML(fakeXML);

  std::string name = "descent";

  app_log() << "Test of engineData construction" << std::endl;
  engineData testData;

  testData.descentEngine = descentEnginePtr;
  testData.method = name;

  bool test_val = testData.descentEngine->targetingExcited();
  REQUIRE(testData.method=="descent");
  REQUIRE(test_val==false);


}
} // namespace qmcplusplus
