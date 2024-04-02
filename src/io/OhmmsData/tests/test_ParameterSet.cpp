//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsData/ParameterSet.h"

TEST_CASE("ParameterSet", "[xml]")
{
  const char* content = R"(
<simulation>
   <parameter name="p1">1</parameter>
   <parameter name="p4"> -2 -3 </parameter>
   <p2>2</p2>
</simulation>)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(content);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();
  ParameterSet param;
  int p1_val = 0;
  int p2_val = 0;
  int p3_val = 0;
  qmcplusplus::astring p4_str;
  param.add(p1_val, "p1");
  param.add(p2_val, "p2");
  param.add(p3_val, "p3");
  param.add(p4_str, "p4");
  param.put(root);

  CHECK(p1_val == 1);
  CHECK(p2_val == 2);
  CHECK(p3_val == 0);

  auto p4_vals = qmcplusplus::convertStrToVec<int>(p4_str.s);
  REQUIRE(p4_vals.size() == 2);
  CHECK(p4_vals[0] == -2);
  CHECK(p4_vals[1] == -3);
}

TEST_CASE("ParameterSet_bool", "[xml]")
{
  {
    const char* content = R"(
      <simulation>
       <parameter name="p1"> yes </parameter>
       <parameter name="p2"> no </parameter>
      </simulation>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    ParameterSet param;
    bool p1_val = false;
    bool p2_val = true;
    bool p3_val = false;
    qmcplusplus::astring p4_str;
    param.add(p1_val, "p1");
    param.add(p2_val, "p2");
    param.add(p3_val, "p3", {true, false});
    param.put(doc.getRoot());

    CHECK(p1_val == true);
    CHECK(p2_val == false);
    CHECK(p3_val == true);
  }

  {
    const char* content = R"(
      <simulation>
        <parameter name="p1">  </parameter>
      </simulation>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    ParameterSet param;
    bool p1_val = false;
    qmcplusplus::astring p4_str;
    param.add(p1_val, "p1");
    CHECK_THROWS_WITH(param.put(doc.getRoot()), "p1 requires a single value input.");
  }

  {
    const char* content = R"(
      <simulation>
        <parameter name="p1"> yes no </parameter>
      </simulation>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    ParameterSet param;
    bool p1_val = false;
    qmcplusplus::astring p4_str;
    param.add(p1_val, "p1");
    CHECK_THROWS_WITH(param.put(doc.getRoot()), "p1 only accepts a single value input.");
  }

  {
    const char* content = R"(
      <simulation>
        <parameter name="p1"> here </parameter>
      </simulation>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    ParameterSet param;
    bool p1_val = false;
    qmcplusplus::astring p4_str;
    param.add(p1_val, "p1");
    CHECK_THROWS_WITH(param.put(doc.getRoot()), "p1 only accepts 'yes'/'no'/'true'/'false' but the input value is 'here'.");
  }
}
