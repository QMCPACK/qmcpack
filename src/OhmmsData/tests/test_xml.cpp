//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

TEST_CASE("read_file", "[xml]")
{
  Libxml2Document doc;
  bool okay = doc.parse("bad.xml");
  REQUIRE(okay == false);

  okay = doc.parse("simple.xml");
  REQUIRE(okay == true);
}

TEST_CASE("parseString", "[xml]")
{
    string s1(
"<?xml version=\"1.0\"?> \
<simulation> \
</simulation> "
);

  Libxml2Document doc;
  bool okay = doc.parseFromString(s1);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();
  REQUIRE(root != NULL);

  REQUIRE((char *)root->name == string("simulation"));

  string root_name;
  getNodeName(root_name, root);

  REQUIRE(root_name == "simulation");

}

TEST_CASE("putContent", "[xml]")
{
  string s1(
"<?xml version=\"1.0\"?> \
<simulation> \
  <item1>2</item1> \
  <item2>3.5</item2> \
  <item3>4.0</item3> \
  <item4>4.0 5.1</item4> \
  <item5>4.0 5.1,</item5> \
</simulation>");

  Libxml2Document doc;
  bool okay = doc.parseFromString(s1);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();

  int a;
  xmlNodePtr item = xmlFirstElementChild(root);
  REQUIRE(string((char *)item->name) == "item1");
  putContent(a, item);
  REQUIRE(a == 2);

  double b;
  item = xmlNextElementSibling(item);
  REQUIRE(string((char *)item->name) == "item2");
  putContent(b, item);
  REQUIRE(b == 3.5);

  float c;
  putContent(c, item);
  REQUIRE(c == Approx(3.5));

  vector<double> d;
  item = xmlNextElementSibling(item);
  REQUIRE(string((char *)item->name) == "item3");
  putContent(d, item);
  REQUIRE(d.size() == 1);
  vector<double> e;
  item = xmlNextElementSibling(item);
  REQUIRE(string((char *)item->name) == "item4");
  putContent(e, item);
  REQUIRE(e.size() == 2);

  vector<double> f;
  item = xmlNextElementSibling(item);
  REQUIRE(string((char *)item->name) == "item5");
  // Will hang, don't test for now
  //putContent(f, item);
  //REQUIRE(f.size() == 2);
}

TEST_CASE("AttributeSet", "[xml]")
{
  const char *content = " \
<simulation name=\"here\"> \
</simulation>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(content);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();
  OhmmsAttributeSet pattrib;
  string name="default_name";
  string other="default";
  pattrib.add(name, "name");
  pattrib.add(other, "other");
  pattrib.put(root);

  REQUIRE(name == "here");
  REQUIRE(other == "default");
}


TEST_CASE("ParameterSet", "[xml]")
{
  const char *content = " \
<simulation> \
   <parameter name=\"p1\">1</parameter> \
   <p2>2</p2> \
</simulation>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(content);
  REQUIRE(okay == true);

  xmlNodePtr root = doc.getRoot();
  ParameterSet param;
  int p1_val = 0;
  int p2_val = 0;
  int p3_val = 0;
  param.add(p1_val, "p1", "none");
  param.add(p2_val, "p2", "none");
  param.add(p3_val, "p3", "none");
  param.put(root);

  REQUIRE(p1_val == 1);
  REQUIRE(p2_val == 2);
  REQUIRE(p3_val == 0);
}

TEST_CASE("write_file", "[xml]")
{
  Libxml2Document doc;
  doc.newDoc("root");

  xmlNodePtr node1 = doc.addChild(doc.getRoot(), "node1");
  doc.addChild(node1, "value1", 1);
  doc.addChild(node1, "value2", 3.2);
  doc.dump("tmp.out.xml");
}

