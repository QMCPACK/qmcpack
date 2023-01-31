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
#include "OhmmsData/AttributeSet.h"
#include <string>

using std::string;

TEST_CASE("AttributeSet", "[xml]")
{
  const char* content = R"(
    <simulation name="here" deprecated_tag="lmn">
    </simulation>)";
  Libxml2Document doc;
  bool okay = doc.parseFromString(content);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();
  OhmmsAttributeSet pattrib;
  string name  = "default_name";
  string other = "default";
  string deprecated_tag;
  pattrib.add(name, "name");
  pattrib.add(other, "other");
  pattrib.add(deprecated_tag, "deprecated_tag", {"abc", "def"}, TagStatus::DEPRECATED);
  CHECK_THROWS_WITH(pattrib.put(doc.getRoot()), Catch::Matchers::Contains("is not valid"));

  REQUIRE(name == "here");
  REQUIRE(other == "default");
  REQUIRE(deprecated_tag == "lmn");
}

TEST_CASE("AttributeSet_bool", "[xml]")
{
  {
    const char* content = R"(<simulation/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == false);
  }

  {
    const char* content = R"(<simulation use_feature=" yes  "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == true);
  }

  {
    const char* content = R"(<simulation use_feature=" no  "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == false);
  }

  {
    const char* content = R"(<simulation use_feature=" YES  "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == true);
  }

  {
    const char* content = R"(<simulation use_feature=" No"/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == false);
  }

  {
    const char* content = R"(<simulation use_feature=" true  "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == true);
  }

  {
    const char* content = R"(<simulation use_feature="false"/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    pattrib.put(doc.getRoot());
    CHECK(use_feature == false);
  }

  {
    const char* content = R"(<simulation use_feature=" "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    CHECK_THROWS_WITH(pattrib.put(doc.getRoot()), "use_feature requires a single value input.");
  }

  {
    const char* content = R"(<simulation use_feature=" no a "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    CHECK_THROWS_WITH(pattrib.put(doc.getRoot()), "use_feature only accepts a single value input.");
  }

  {
    const char* content = R"(<simulation use_feature=" here  "/>)";
    Libxml2Document doc;
    bool okay = doc.parseFromString(content);
    REQUIRE(okay == true);

    OhmmsAttributeSet pattrib;
    bool use_feature = true;
    pattrib.add(use_feature, "use_feature", {false, true});
    CHECK_THROWS_WITH(pattrib.put(doc.getRoot()),
                      "use_feature only accepts 'yes'/'no'/'true'/'false' but the input value is 'here'.");
  }
}
