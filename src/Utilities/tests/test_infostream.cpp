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

#include "Utilities/InfoStream.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

TEST_CASE("InfoStream basic", "[utilities]")
{
  std::ostringstream out;
  InfoStream info(&out);
  info.getStream() << "test";
  REQUIRE(out.str() == "test");

  info.pause();
  info.getStream() << "empty";
  REQUIRE(out.str() == "test");

  info.resume();
  info.getStream() << "second";
  REQUIRE(out.str() == "testsecond");
}

TEST_CASE("InfoStream redirect to new stream", "[utilities]")
{
  std::ostringstream out1;
  std::ostringstream out2;
  InfoStream info1(&out1);
  InfoStream info2(&out2);

  info1.getStream() << "test1";
  info2.getStream() << "test2";

  info1.redirectToSameStream(info2);
  info1.getStream() << "test3";

  REQUIRE(out2.str() == "test2test3");
  REQUIRE(out1.str() == "test1");
}

TEST_CASE("InfoStream redirect to file", "[utilities]")
{
  std::ostringstream out;
  InfoStream *info = new InfoStream(&out);
  info->redirectToFile("test_infostream.log");
  info->getStream() << "test";
  delete(info);

  std::ifstream in("test_infostream.log");
  std::string s;
  in >> s;
  REQUIRE(s == "test");
}

TEST_CASE("InfoStream double pause", "[utilities]")
{
  std::ostringstream out;
  InfoStream info(&out);
  info.pause();
  info.pause();
  info.resume();
  info.getStream() << "test";
  REQUIRE(out.str() == "test");
}

TEST_CASE("InfoStream shutOff", "[utilities]")
{
  std::ostringstream out;
  InfoStream info(&out);
  info.shutOff();
  info.pause();
  info.resume();
  info.getStream() << "test";
  REQUIRE(out.str() == "");
}

TEST_CASE("InfoStream operator", "[utilities]")
{
  std::ostringstream out;
  InfoStream info(&out);
  info << "test";
  info << 1;
  REQUIRE(out.str() == "test1");

}
