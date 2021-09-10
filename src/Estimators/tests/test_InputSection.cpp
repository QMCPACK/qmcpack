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

class TestInput : public InputSection
{
public:
  TestInput()
  {
    section_name   = "Test";
    attributes     = {"name","samples","kmax","full"};
    parameters     = {"label","count","width","rational"};
    required       = {"count","full"};
    strings        = {"name","label"};
    integers       = {"samples","count"};
    reals          = {"kmax","width"};
    bools          = {"full","rational"};
    default_values = { {"name"    , std::string("demo") },
                       {"samples" , 20     },
                       {"width"   , 1.0    },
                       {"rational", false  }};
  };
};


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
  {// test read: minimum required attributes and parameters

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
    CHECK(ti.get<bool>("full")==false);
    CHECK(ti.get<int>("count")==15);
    CHECK(ti.get<std::string>("name")=="demo");
    CHECK(ti.get<int>("samples")==20);
    CHECK(ti.get<RealType>("width")==Approx(1.0));
    CHECK(ti.get<bool>("rational")==false);
  }


  {// test read: complete attributes and parameters

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
    CHECK(ti.get<std::string>("name")=="alice");
    CHECK(ti.get<int>("samples")==10);
    CHECK(ti.get<RealType>("kmax")==Approx(3.0));
    CHECK(ti.get<bool>("full")==false);
    CHECK(ti.get<std::string>("label")=="relative");
    CHECK(ti.get<int>("count")==15);
    CHECK(ti.get<RealType>("width")==Approx(2.5));
    CHECK(ti.get<bool>("rational")==true);
  }


  {// test read: invalid inputs
   //   cases that should result in thrown exception:
   //     required attribute/parameter is missing
   //     unrecognized attribute/parameter encountered

    std::unordered_map<std::string,const char*> invalid_inputs = {
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

    for (auto &[label, xml] : invalid_inputs)
    {
      // parse xml doc
      Libxml2Document doc;
      bool okay = doc.parseFromString(xml);
      REQUIRE(okay);
      xmlNodePtr cur = doc.getRoot();
      
      // read xml
      TestInput ti;
      //app_log()<<"trying xml read for "<<label<<"\n";
      try
      {
        ti.readXML(cur);
      }
      catch(std::runtime_error& exception)
      {
        //app_log()<<label<<" is caught correctly\n";
      }
    }

  }

}


TEST_CASE("InputSection::init", "[estimators]")
{
  {// test init: minimum required attributes and parameters
    // initialize
    TestInput ti;
    ti.init({{"full",false},{"count",15}});

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
    CHECK(ti.get<bool>("full")==false);
    CHECK(ti.get<int>("count")==15);
    CHECK(ti.get<std::string>("name")=="demo");
    CHECK(ti.get<int>("samples")==20);
    CHECK(ti.get<RealType>("width")==Approx(1.0));
    CHECK(ti.get<bool>("rational")==false);
  }


  {// test init: complete attributes and parameters
    // initialize
    TestInput ti;
    ti.init({
        {"name"    ,std::string("alice")},
        {"samples" ,10   },
        {"kmax"    ,3.0  },
        {"full"    ,false},
        {"label"   ,std::string("relative")},
        {"count"   ,15   },
        {"width"   ,2.5  },
        {"rational",true },
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
    CHECK(ti.get<std::string>("name")=="alice");
    CHECK(ti.get<int>("samples")==10);
    CHECK(ti.get<RealType>("kmax")==Approx(3.0));
    CHECK(ti.get<bool>("full")==false);
    CHECK(ti.get<std::string>("label")=="relative");
    CHECK(ti.get<int>("count")==15);
    CHECK(ti.get<RealType>("width")==Approx(2.5));
    CHECK(ti.get<bool>("rational")==true);
  }

}




TEST_CASE("InputSection::get", "[estimators]")
{
  TestInput ti;
  ti.init({
      {"name"    ,std::string("alice")},
      {"samples" ,10   },
      {"kmax"    ,3.0  },
      {"full"    ,false},
      {"label"   ,std::string("relative")},
      {"count"   ,15   },
      {"width"   ,2.5  },
      {"rational",true },
      });
  
  // invalid type access results in thrown exception
  try
  {
    auto v = ti.get<int>("name");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<RealType>("samples");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<bool>("kmax");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<std::string>("full");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<RealType>("label");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<bool>("count");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<std::string>("width");
  }
  catch(std::bad_cast& exception)
  {
  }

  try
  {
    auto v = ti.get<int>("rational");
  }
  catch(std::bad_cast& exception)
  {
  }

}


} // namespace qmcplusplus
