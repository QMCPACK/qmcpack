//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "InputSection.h"
#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>


namespace qmcplusplus
{
using RealType = QMCTraits::RealType;

// clang-format: off
class TestInput : public InputSection
{
public:
  TestInput()
  {
    section_name   = "Test";
    attributes     = {"name", "samples", "kmax", "full"};
    parameters     = {"label", "count", "width", "rational"};
    required       = {"count", "full"};
    strings        = {"name", "label"};
    integers       = {"samples", "count"};
    reals          = {"kmax", "width"};
    bools          = {"full", "rational"};
    default_values = {{"name", std::string("demo")},
                      {"samples", int(20)},
                      {"width", RealType(1.0)},
                      {"rational", bool(false)}};
  };
};
// clang-format: on


TEST_CASE("InputSection::InputSection", "[estimators]")
{
  // empty input section
  TestInput ti;
  // no values assigned yet
  CHECK(!ti.has("name"));
  CHECK(!ti.has("samples"));
  CHECK(!ti.has("kmax"));
  CHECK(!ti.has("full"));
  CHECK(!ti.has("label"));
  CHECK(!ti.has("count"));
  CHECK(!ti.has("width"));
  CHECK(!ti.has("rational"));
}


TEST_CASE("InputSection::readXML", "[estimators]")
{
  SECTION("minimum required attributes and parameters")
  {
    // clang-format: off
    const char* xml = R"(
<test full="no">
  <parameter name="count"> 15 </parameter>
</test>
)";
    // clang-format: on

    // parse xml doc
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml);
    REQUIRE(okay);
    xmlNodePtr cur = doc.getRoot();

    // read xml
    TestInput ti;
    ti.readXML(cur);

    // assigned from xml
    CHECK(ti.has("full"));
    CHECK(ti.has("count"));
    // assigned from defaults
    CHECK(ti.has("name"));
    CHECK(ti.has("samples"));
    CHECK(ti.has("width"));
    CHECK(ti.has("rational"));
    // unassigned
    CHECK(!ti.has("kmax"));
    CHECK(!ti.has("label"));
    // check value correctness
    CHECK(ti.get<bool>("full") == false);
    CHECK(ti.get<int>("count") == 15);
    CHECK(ti.get<std::string>("name") == "demo");
    CHECK(ti.get<int>("samples") == 20);
    CHECK(ti.get<RealType>("width") == Approx(1.0));
    CHECK(ti.get<bool>("rational") == false);
  }


  SECTION("complete attributes and parameters")
  {
    // clang-format: off
    const char* xml = R"(
<test name="alice" samples="10" kmax="3.0" full="no">
  <parameter name="label"   >  relative  </parameter>
  <parameter name="count"   >  15        </parameter>
  <parameter name="width"   >  2.5       </parameter>
  <parameter name="rational">  yes       </parameter>
</test>
)";
    // clang-format: on

    // parse xml doc
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml);
    REQUIRE(okay);
    xmlNodePtr cur = doc.getRoot();

    // read xml
    TestInput ti;
    ti.readXML(cur);

    // assigned from xml
    CHECK(ti.has("name"));
    CHECK(ti.has("samples"));
    CHECK(ti.has("kmax"));
    CHECK(ti.has("full"));
    CHECK(ti.has("label"));
    CHECK(ti.has("count"));
    CHECK(ti.has("width"));
    CHECK(ti.has("rational"));
    // check value correctness
    CHECK(ti.get<std::string>("name") == "alice");
    CHECK(ti.get<int>("samples") == 10);
    CHECK(ti.get<RealType>("kmax") == Approx(3.0));
    CHECK(ti.get<bool>("full") == false);
    CHECK(ti.get<std::string>("label") == "relative");
    CHECK(ti.get<int>("count") == 15);
    CHECK(ti.get<RealType>("width") == Approx(2.5));
    CHECK(ti.get<bool>("rational") == true);
  }


  SECTION("invalid inputs")
  {
    //   cases that should result in thrown exception:
    //     required attribute/parameter is missing
    //     unrecognized attribute/parameter encountered

    std::unordered_map<std::string, const char*> invalid_inputs = {
        {"missing_attribute", R"(
<test>
  <parameter name="count"> 15 </parameter>
</test>
)"},
        {"missing_parameter", R"(
<test full="no"/>
)"},
        {"foreign_attribute", R"(
<test full="no" area="51">
  <parameter name="count"> 15 </parameter>
</test>
)"},
        {"foreign_parameter", R"(
<test full="no">
  <parameter name="count"> 15 </parameter>
  <parameter name="area" > 51 </parameter>
</test>
)"},
    };

    for (auto& [label, xml] : invalid_inputs)
    {
      // parse xml doc
      Libxml2Document doc;
      bool okay = doc.parseFromString(xml);
      REQUIRE(okay);
      xmlNodePtr cur = doc.getRoot();

      // read xml
      TestInput ti;
      CHECK_THROWS_AS(ti.readXML(cur), std::runtime_error);
    }
  }
}


TEST_CASE("InputSection::init", "[estimators]")
{
  SECTION("minimum required attributes and parameters")
  {
    // initialize
    TestInput ti;
    ti.init({{"full", bool(false)}, {"count", int(15)}});

    // assigned from initializer-list
    CHECK(ti.has("full"));
    CHECK(ti.has("count"));
    // assigned from defaults
    CHECK(ti.has("name"));
    CHECK(ti.has("samples"));
    CHECK(ti.has("width"));
    CHECK(ti.has("rational"));
    // unassigned
    CHECK(!ti.has("kmax"));
    CHECK(!ti.has("label"));
    // check value correctness
    CHECK(ti.get<bool>("full") == false);
    CHECK(ti.get<int>("count") == 15);
    CHECK(ti.get<std::string>("name") == "demo");
    CHECK(ti.get<int>("samples") == 20);
    CHECK(ti.get<RealType>("width") == Approx(1.0));
    CHECK(ti.get<bool>("rational") == false);
  }


  SECTION("complete attributes and parameters")
  {
    // initialize
    TestInput ti;
    ti.init({
        {"name", std::string("alice")},
        {"samples", int(10)},
        {"kmax", RealType(3.0)},
        {"full", bool(false)},
        {"label", std::string("relative")},
        {"count", int(15)},
        {"width", RealType(2.5)},
        {"rational", bool(true)},
    });

    // assigned from initializer-list
    CHECK(ti.has("name"));
    CHECK(ti.has("samples"));
    CHECK(ti.has("kmax"));
    CHECK(ti.has("full"));
    CHECK(ti.has("label"));
    CHECK(ti.has("count"));
    CHECK(ti.has("width"));
    CHECK(ti.has("rational"));
    // check value correctness
    CHECK(ti.get<std::string>("name") == "alice");
    CHECK(ti.get<int>("samples") == 10);
    CHECK(ti.get<RealType>("kmax") == Approx(3.0));
    CHECK(ti.get<bool>("full") == false);
    CHECK(ti.get<std::string>("label") == "relative");
    CHECK(ti.get<int>("count") == 15);
    CHECK(ti.get<RealType>("width") == Approx(2.5));
    CHECK(ti.get<bool>("rational") == true);
  }


  SECTION("invalid type assignment")
  {
    TestInput ti;
    CHECK_THROWS_AS(ti.init({{"full", bool(false)}, {"count", RealType(15.)}}), std::bad_cast);
  }
}


TEST_CASE("InputSection::get", "[estimators]")
{
  TestInput ti;
  ti.init({
      {"name", std::string("alice")},
      {"samples", int(10)},
      {"kmax", RealType(3.0)},
      {"full", bool(false)},
      {"label", std::string("relative")},
      {"count", int(15)},
      {"width", RealType(2.5)},
      {"rational", bool(true)},
  });

  // invalid type access results in thrown exception
  CHECK_THROWS_AS(ti.get<int>("name"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<RealType>("samples"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<bool>("kmax"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<std::string>("full"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<RealType>("label"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<bool>("count"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<std::string>("width"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<int>("rational"), std::bad_cast);
}


} // namespace qmcplusplus
