//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/OutputManager.h"

namespace qmcplusplus {

// redirect streams to stringstream for testing
static std::ostringstream summary_out;
static std::ostringstream app_out;
static std::ostringstream err_out;
static std::ostringstream debug_out;

void reset_string_output()
{
  summary_out.str("");
  app_out.str("");
  err_out.str("");
  debug_out.str("");
}

void init_string_output()
{
  infoSummary.setStream(&summary_out);
  infoLog.setStream(&app_out);
  infoError.setStream(&err_out);
  infoDebug.setStream(&debug_out);
  reset_string_output();
}


TEST_CASE("OutputManager basic", "[utilities]")
{

  OutputManagerClass om;

  // Low verbosity
  om.setVerbosity(Verbosity::LOW);
  REQUIRE(om.isDebugActive() == false);
  REQUIRE(om.isHighActive() == false);

  // High verbosity
  om.setVerbosity(Verbosity::HIGH);
  REQUIRE(om.isDebugActive() == false);
  REQUIRE(om.isHighActive() == true);

  // Debug verbosity
  om.setVerbosity(Verbosity::DEBUG);
  REQUIRE(om.isDebugActive() == true);
  REQUIRE(om.isHighActive() == true);
}

TEST_CASE("OutputManager output", "[utilities]")
{

  init_string_output();

  // Low verbosity
  outputManager.setVerbosity(Verbosity::LOW);

  app_summary() << "test1";
  app_log() << "test2";
  app_debug() << "test3";
  app_error() << "test4";
  REQUIRE(summary_out.str() == "test1");
  REQUIRE(app_out.str() == "");
  REQUIRE(debug_out.str() == "");
  REQUIRE(err_out.str() == "ERROR test4");

  reset_string_output();

  // High verbosity
  outputManager.setVerbosity(Verbosity::HIGH);

  app_summary() << "test5";
  app_log() << "test6";
  app_debug() << "test7";
  app_error() << "test8";
  REQUIRE(summary_out.str() == "test5");
  REQUIRE(app_out.str() == "test6");
  REQUIRE(debug_out.str() == "");
  REQUIRE(err_out.str() == "ERROR test8");

  reset_string_output();

  // Debug verbosity
  outputManager.setVerbosity(Verbosity::DEBUG);

  app_summary() << "test9";
  app_log() << "testA";
  app_debug() << "testB";
  app_error() << "testC";
  REQUIRE(summary_out.str() == "test9");
  REQUIRE(app_out.str() == "testA");
  REQUIRE(debug_out.str() == "testB");
  REQUIRE(err_out.str() == "ERROR testC");
}

TEST_CASE("OutputManager pause", "[utilities]")
{
  init_string_output();

  outputManager.pause();

  app_summary() << "test1";
  app_log() << "test2";

  REQUIRE(summary_out.str() == "");
  REQUIRE(app_out.str() == "");

  reset_string_output();

  outputManager.resume();
  app_summary() << "test3";
  app_log() << "test4";

  REQUIRE(summary_out.str() == "test3");
  REQUIRE(app_out.str() == "test4");
}


TEST_CASE("OutputManager shutoff", "[utilities]")
{
  init_string_output();

  outputManager.shutOff();

  app_summary() << "test1";
  app_log() << "test2";
  app_error() << "test3";
  app_debug() << "test4";

  REQUIRE(summary_out.str() == "");
  REQUIRE(app_out.str() == "");
  REQUIRE(err_out.str() == "");
  REQUIRE(debug_out.str() == "");
}

}
