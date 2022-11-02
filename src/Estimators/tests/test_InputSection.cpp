//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//                    Peter W.  Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "InputSection.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Utilities/string_utils.h"
#include "Utilities/ModernStringUtils.hpp"
#include "EstimatorInput.h"

#include <stdio.h>
#include <sstream>


namespace qmcplusplus
{
using RealType = QMCTraits::RealType;

enum class TestEnum1
{
  VALUE1,
  VALUE2
};

enum class TestEnum2
{
  VALUE1,
  VALUE2
};

std::unordered_map<std::string, std::any> lookup_input_enum_value{{"testenum1-value1", TestEnum1::VALUE1},
                                                                  {"testenum1-value2", TestEnum1::VALUE2},
                                                                  {"testenum2-value1", TestEnum2::VALUE1},
                                                                  {"testenum2-value2", TestEnum2::VALUE2}};

// clang-format: off
class TestInputSection : public InputSection
{
public:
  TestInputSection()
  {
    section_name   = "Test";
    attributes     = {"name", "samples", "kmax", "full"};
    parameters     = {"label", "count", "width", "rational", "testenum1", "testenum2", "sposets", "center"};
    required       = {"count", "full"};
    strings        = {"name", "label"};
    multi_strings  = {"sposets"};
    integers       = {"samples", "count"};
    reals          = {"kmax", "width"};
    positions      = {"center"};
    bools          = {"full", "rational"};
    enums          = {"testenum1", "testenum2"};
    default_values = {{"name", std::string("demo")},
                      {"samples", int(20)},
                      {"width", RealType(1.0)},
                      {"rational", bool(false)}};
  };
  // clang-format: on

  std::any assignAnyEnum(const std::string& name) const override
  {
    return lookupAnyEnum(name, get<std::string>(name), lookup_input_enum_value);
  }
};

class TestAssignAnyEnum : public InputSection
{
public:
  TestAssignAnyEnum() { enums = {"testenum"}; }
};


TEST_CASE("InputSection::InputSection", "[estimators]")
{
  // empty input section
  TestInputSection ti;
  // no values assigned yet
  CHECK(!ti.has("name"));
  CHECK(!ti.has("samples"));
  CHECK(!ti.has("kmax"));
  CHECK(!ti.has("full"));
  CHECK(!ti.has("label"));
  CHECK(!ti.has("count"));
  CHECK(!ti.has("width"));
  CHECK(!ti.has("rational"));
  CHECK(!ti.has("sposets"));
}

TEST_CASE("InputSection::assignAnyEnum", "[estimators]")
{
  enum class SomeEnum
  {
    VALUE1,
    VALUE2
  };

  TestAssignAnyEnum taae;
  CHECK_THROWS_AS(taae.get<SomeEnum>("testenum"), std::runtime_error);
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
    TestInputSection ti;
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
    CHECK(!ti.has("sposets"));
    CHECK(!ti.has("center"));
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
  <parameter name="testenum1"> Value1 </parameter>
  <parameter name="testenum2"> Value2 </parameter>
  <parameter name="sposets"> spo1 spo2 </parameter>
  <parameter name="center"> 0.0 0.0 0.1 </parameter>
</test>
)";
    // clang-format: on

    // parse xml doc
    Libxml2Document doc;
    bool okay = doc.parseFromString(xml);
    REQUIRE(okay);
    xmlNodePtr cur = doc.getRoot();

    // read xml
    TestInputSection ti;
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
    CHECK(ti.has("sposets"));
    // check value correctness
    CHECK(ti.get<std::string>("name") == "alice");
    CHECK(ti.get<int>("samples") == 10);
    CHECK(ti.get<RealType>("kmax") == Approx(3.0));
    CHECK(ti.get<bool>("full") == false);
    CHECK(ti.get<std::string>("label") == "relative");
    CHECK(ti.get<int>("count") == 15);
    CHECK(ti.get<RealType>("width") == Approx(2.5));
    CHECK(ti.get<bool>("rational") == true);
    CHECK(ti.get<TestEnum1>("testenum1") == TestEnum1::VALUE1);
    CHECK(ti.get<TestEnum2>("testenum2") == TestEnum2::VALUE2);
    CHECK(ti.get<std::vector<std::string>>("sposets") == std::vector<std::string>{"spo1", "spo2"});
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
      TestInputSection ti;
      CHECK_THROWS_AS(ti.readXML(cur), std::runtime_error);
    }
  }
}


TEST_CASE("InputSection::init", "[estimators]")
{
  SECTION("minimum required attributes and parameters")
  {
    // initialize
    TestInputSection ti;
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
    TestInputSection ti;
    ti.init({{"name", std::string("alice")},
             {"samples", int(10)},
             {"kmax", RealType(3.0)},
             {"full", bool(false)},
             {"label", std::string("relative")},
             {"count", int(15)},
             {"width", RealType(2.5)},
             {"rational", bool(true)},
             {"sposets", std::vector<std::string>{"spo1", "spo2"}},
             {"center", InputSection::Position(0.0, 0.0, 0.1)}});

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
    CHECK(ti.get<std::vector<std::string>>("sposets") == std::vector<std::string>{"spo1", "spo2"});
    CHECK(ti.get<InputSection::Position>("center") == InputSection::Position(0.0, 0.0, 0.1));
  }


  SECTION("invalid type assignment")
  {
    TestInputSection ti;
    CHECK_THROWS_AS(ti.init({{"full", bool(false)}, {"count", RealType(15.)}}), std::bad_cast);
  }
}


TEST_CASE("InputSection::get", "[estimators]")
{
  TestInputSection ti;
  ti.init({{"name", std::string("alice")},
           {"samples", int(10)},
           {"kmax", RealType(3.0)},
           {"full", bool(false)},
           {"label", std::string("relative")},
           {"count", int(15)},
           {"width", RealType(2.5)},
           {"rational", bool(true)},
           {"sposets", std::vector<std::string>{"spo1", "spo2"}},
           {"center", InputSection::Position(0.0, 0.0, 0.1)}});

  // invalid type access results in thrown exception
  CHECK_THROWS_AS(ti.get<int>("name"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<RealType>("samples"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<bool>("kmax"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<std::string>("full"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<RealType>("label"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<bool>("count"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<std::string>("width"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<int>("rational"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<std::string>("sposets"), std::bad_cast);
  CHECK_THROWS_AS(ti.get<RealType>("center"), std::bad_cast);
}

class CustomTestInput : public InputSection
{
public:
  struct WeirdStuff {
    std::string letters;
    std::array<int,3> numbers;
  };
  using Repeater = std::vector<std::pair<std::string, std::string>>;
  CustomTestInput()
  {
    section_name = "Test";
    attributes   = {"name", "samples", "kmax", "full"};
    parameters   = {"label", "count", "width"};
    strings      = {"name", "label"};
    reals        = {"kmax"};
    integers     = {"samples", "count"};
    bools        = {"full"};
    custom       = {"weird_stuff", "repeater"};
  }
  void setFromStreamCustom(const std::string& ename, const std::string& name, std::istringstream& svalue) override
  {
    if (ename == "weird_stuff")
    {
      WeirdStuff ws;
      svalue >> ws.letters;
      svalue >> ws.numbers[0];
      svalue >> ws.numbers[1];
      svalue >> ws.numbers[2];
      values_[name] = ws;
    }
    else if (ename == "repeater")
    {
      std::string compound;
      svalue >> compound;
      auto split_vstrv = split(compound, ":"); 
      if (has(name))
	std::any_cast<Repeater>(&(values_[name]))->emplace_back(split_vstrv[0],split_vstrv[1]);
      else
	values_[name] = Repeater{{split_vstrv[0],split_vstrv[1]}};
    }
    else
      throw std::runtime_error("bad name passed or custom setFromStream not implemented in derived class.");
  }
};

class FailCustomTestInput : public InputSection
{
public:
  struct WeirdStuff {
    std::string letters;
    std::array<int,3> numbers;
  };
  using Repeater = std::vector<std::pair<std::string, std::string>>;
  FailCustomTestInput()
  {
    section_name = "Test";
    attributes   = {"name", "samples", "kmax", "full"};
    parameters   = {"label", "count", "width"};
    strings      = {"name", "label"};
    reals        = {"kmax"};
    integers     = {"samples", "count"};
    bools        = {"full"};
    custom       = {"weird_stuff", "repeater"};
  }
};

TEST_CASE("InputSection::custom", "[estimators]")
{
  std::string_view xml = R"XML(
<test name="alice" samples="10" kmax="3.0" full="no">
  <parameter name="label"   >  relative  </parameter>
  <parameter name="count"   >  15        </parameter>
  <weird_stuff name="weird"> XQ 10 20 10 </weird_stuff>
  <Repeater> first:something </Repeater>
  <Repeater> second:else </Repeater>
</test>
)XML";

  Libxml2Document doc;
  bool okay = doc.parseFromString(xml);
  xmlNodePtr cur = doc.getRoot();
  CustomTestInput cti;
  cti.readXML(cur);

  auto ws = cti.get<decltype(cti)::WeirdStuff>("weird");
  CHECK(ws.letters == "XQ");
  std::array<int,3> exp_numbers{10, 20, 10};
  CHECK(ws.numbers == exp_numbers);

  auto repeater = cti.get<decltype(cti)::Repeater>("repeater");
  decltype(cti)::Repeater exp_repeater{{"first","something"},{"second","else"}};
  CHECK(repeater == exp_repeater);

  FailCustomTestInput fcti;
  CHECK_THROWS_AS(fcti.readXML(cur), std::runtime_error);
}

/** Generic delgate input class.
 *  This feature is targeted at supporting things like ReferencePointsInput, SpaceGridInput etc.
 *
 *  The delegate input requires:
 *  A constructor that takes an xmlNodePtr
 *  A factoryFunction that matches the signature discussed in InputSection::DelegateHandler
 */

class AnotherInput
{
public:
  class AnotherInputSection : public InputSection
  {
  public:
    AnotherInputSection()
    {
      section_name = "AnotherInputSection";
      attributes   = {"name", "optional"};
      strings      = {"name", "optional"};
    }
  };

  AnotherInput(xmlNodePtr cur)
  {
    input_section_.readXML(cur);
    auto setIfInInput = LAMBDA_setIfInInput;
    setIfInInput(name_, "name");
    setIfInInput(optional_, "optional");
    input_str_ = XMLNodeString{cur};
  }

  std::vector<std::string_view> getTokens() const { return modernstrutil::split(input_str_, " "); }
  const std::string& get_name() const { return name_; };

private:
  std::string name_ = "AnotherInputSection";
  AnotherInputSection input_section_;
  std::string input_str_;
  std::string optional_;
  std::vector<std::string_view> tokens_;
};

std::any makeAnotherInput(xmlNodePtr cur, std::string& value_label)
{
  AnotherInput another_input{cur};
  value_label = another_input.get_name();
  return another_input;
}

/** Generic input section delegating some inputs to other input classes.
 *  an example would be EnergyDensityInput
 */
class DelegatingInput
{
  class DelegatingInputSection : public InputSection
  {
  public:
    DelegatingInputSection()
    {
      section_name = "DelegateTest";
      attributes   = {"name", "full"};
      parameters   = {"label", "count", "width"};
      strings      = {"name", "label"};
      integers     = {"count"};
      bools        = {"full"};
      delegates    = {"anotherinput"};
      InputSection::registerDelegate("anotherinput", makeAnotherInput);
    }
  };

public:
  DelegatingInput(xmlNodePtr cur) { dins_.readXML(cur); }

  // rather than write a real input class we just pass through InputSection for this test.
  bool has(const std::string& name) const { return dins_.has(name); }
  template<typename T>
  T get(const std::string& name) const
  {
    return dins_.get<T>(name);
  }

private:
  std::string name_ = "DelegatingInput";
  DelegatingInputSection dins_;
};

TEST_CASE("InputSection::Delegate", "[estimators]")
{
  std::string_view xml = R"XML(
<test name="alice" full="no">
  <parameter name="label"   >  relative  </parameter>
  <parameter name="count"   >  15        </parameter>
  <AnotherInput name="ainput"> XQ 10 20 10 </AnotherInput>
</test>
)XML";

  Libxml2Document doc;
  bool okay = doc.parseFromString(xml);
  REQUIRE(okay);
  xmlNodePtr cur = doc.getRoot();

  DelegatingInput dti(cur);

  CHECK(dti.has("ainput"));
  AnotherInput ai(dti.get<AnotherInput>("ainput"));
  std::vector<std::string_view> ref_tokens{"XQ", "10", "20", "10"};
  CHECK(ai.getTokens() == ref_tokens);


  std::string_view xml_duplicate_delegate_name = R"XML(
<test name="alice" full="no">
  <parameter name="label"   >  relative  </parameter>
  <parameter name="count"   >  15        </parameter>
  <AnotherInput name="ainput"> XQ 10 20 10 </AnotherInput>
  <AnotherInput name="ainput"> XQ 10 20 10 </AnotherInput>
</test>
)XML";

  okay = doc.parseFromString(xml_duplicate_delegate_name);
  REQUIRE(okay);
  cur = doc.getRoot();

  CHECK_THROWS_AS(dti = DelegatingInput(cur), UniformCommunicateError);
}

} // namespace qmcplusplus
